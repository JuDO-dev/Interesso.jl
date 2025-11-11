function eval_funcs(optimizer::MOI.ModelLike, funcs::Vector{<:MOI.AbstractFunction})
    var_idxs = MOI.get(optimizer, MOI.ListOfVariableIndices())
    x_vals   = MOI.get(optimizer, MOI.VariablePrimal(), var_idxs)

    nl = MOI.Nonlinear.Model()
    backend = MOI.Nonlinear.SparseReverseMode()

    for f in funcs
        MOI.Nonlinear.add_constraint(nl, f, MOI.EqualTo(0.0))
    end
    
    evaluator = MOI.Nonlinear.Evaluator(nl, backend, var_idxs)
    MOI.initialize(evaluator, Symbol[])
    eval = fill(NaN, length(funcs))
    MOI.eval_constraint(evaluator, eval, x_vals)

    return eval
end

function eval_accuracy(model::Optimizer; q::Integer=10)
    if q < 1
        throw(DomainError(q, "Please ensure q ≥ 1."))
    end

    τ_nodes, τ_weights = FGQ.gausslegendre(q)

    total_res = 0.0

    for phase in model.phases
        mesh = get(model.meshes, phase, nothing)
        if mesh === nothing
            continue
        end

        dif_cons = collect(values(model.dif_cons[phase]))
        alg_cons = collect(values(model.alg_cons[phase]))

        if isempty(dif_cons) && isempty(alg_cons)
            continue
        end

        for mesh_i in get_points_meshes(mesh)
            Δt = 0.5 * (mesh_i.t_b - mesh_i.t_a)
            Σt = 0.5 * (mesh_i.t_b + mesh_i.t_a)

            for (τ, ω) in zip(τ_nodes, τ_weights)
                t = Σt + Δt * τ
                weight = Δt * ω

                for (dif_fun, _) in dif_cons
                    residual = _evaluate_differential_residual(model, dif_fun, t)
                    total_res += weight * abs(residual)
                end

                for (alg_fun, _) in alg_cons
                    residual = _evaluate_dynamic_function(model, alg_fun, t)
                    total_res += weight * abs(residual)
                end
            end
        end
    end

    return total_res
end

function _evaluate_value(
    sol::PiecewiseInterpolant{LagrangeInterpolant},
    t::Float64,
)
    idx = searchsortedlast(sol.pieces_initials, t)
    idx = clamp(idx, 1, length(sol.pieces))
    t_eval = t

    if idx < length(sol.pieces)
        tol = eps(Float64) * max(1.0, maximum(abs, sol.pieces_initials))
        next_initial = sol.pieces_initials[idx + 1]
        if abs(t - next_initial) ≤ tol
            idx += 1
            t_eval = max(t, next_initial)
        end
    end

    return sol.pieces[idx](t_eval)
end

function _get_dyn_var_solution(model::Optimizer, phase::PHS, dyn_var::DYN_VAR)
    sol_by_phase = get(model.sol_dyn_vars, phase, nothing)
    if sol_by_phase === nothing || !haskey(sol_by_phase, dyn_var)
        throw(ArgumentError("No solution stored for dynamic variable $(dyn_var)."))
    end
    return sol_by_phase[dyn_var]
end

function _get_derivative_solution(model::Optimizer, phase::PHS, dyn_var::DYN_VAR)
    derivative_solutions = get(model.sol_derivatives, phase, nothing)
    if derivative_solutions === nothing
        throw(ArgumentError("No derivative solutions stored for phase $(phase)."))
    end
    derivative_index = DOI.Derivative(dyn_var)
    if !haskey(derivative_solutions, derivative_index)
        throw(ArgumentError("No derivative solution stored for dynamic variable $(dyn_var)."))
    end
    return derivative_solutions[derivative_index]
end

function _evaluate_differential_residual(
    model::Optimizer,
    dif_fun::DIF_FUN,
    t::Float64,
)
    derivative_value = _evaluate_dynamic_function(model, DOI.Derivative(dif_fun.dyn_var), t)
    rhs = _evaluate_dynamic_function(model, dif_fun.dyn_fun, t)
    return derivative_value - rhs
end

function _evaluate_dynamic_function(
    model::Optimizer,
    fun::DOI.AbstractDynamicFunction,
    t::Float64,
)
    if fun isa DOI.PhaseIndex
        return t
    elseif fun isa DYN_VAR
        phase = DOI.phase_index(fun)
        sol = _get_dyn_var_solution(model, phase, fun)
        return _evaluate_value(sol, t)
    elseif fun isa DOI.Derivative{DYN_VAR}
        phase = DOI.phase_index(fun)
        derivative_sol = _get_derivative_solution(model, phase, fun.dyn_fun)
        derivative_value = derivative_sol(t)
        scale = _get_phase_duration(model, phase)
        if scale == 0.0
            throw(DomainError(scale, "Phase duration must be nonzero."))
        end
        return derivative_value / scale
    elseif fun isa DOI.LinearDynamicFunction
        total = zero(Float64)
        for term in fun.terms
            total += term.coefficient * _evaluate_dynamic_function(model, term.dyn_var, t)
        end
        return total
    elseif fun isa DOI.PureQuadraticDynamicFunction
        total = zero(Float64)
        for term in fun.terms
            total += term.coefficient *
                _evaluate_dynamic_function(model, term.dyn_var_1, t) *
                _evaluate_dynamic_function(model, term.dyn_var_2, t)
        end
        return total
    elseif fun isa DOI.NonlinearDynamicFunction
        args = [_evaluate_dynamic_argument(model, arg, t) for arg in fun.args]
        op = _resolve_operator(fun.head)
        return op(args...)
    elseif fun isa DOI.ExplicitDifferentialFunction
        return _evaluate_differential_residual(model, fun, t)
    elseif fun isa Real
        return fun
    else
        throw(ArgumentError("Unsupported dynamic function $(typeof(fun))."))
    end
end

function _get_phase_duration(model::Optimizer, phase::PHS)
    time_var = model.time_vars[phase]
    if time_var isa Float64
        return time_var
    end
    return MOI.get(model.inner, MOI.VariablePrimal(), time_var)
end

function _evaluate_dynamic_argument(
    model::Optimizer,
    arg,
    t::Float64,
)
    if arg isa DOI.AbstractDynamicFunction
        return _evaluate_dynamic_function(model, arg, t)
    elseif arg isa MOI.AbstractScalarFunction
        return eval_funcs(model.inner, [arg])[1]
    elseif arg isa MOI.VariableIndex
        return MOI.get(model.inner, MOI.VariablePrimal(), arg)
    elseif arg isa Bool || arg isa Real
        return arg
    else
        throw(ArgumentError("Unsupported dynamic argument $(typeof(arg))."))
    end
end

function _resolve_operator(head::Symbol)
    if head === :ifelse
        return ifelse
    end
    if isdefined(Base, head)
        op = getfield(Base, head)
        if op isa Function
            return op
        end
    end
    if isdefined(MOI.Nonlinear, head)
        op = getfield(MOI.Nonlinear, head)
        if op isa Function
            return op
        end
    end
    throw(ArgumentError("Unsupported nonlinear operator $(head)."))
end