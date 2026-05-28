struct IntervalResidual
    nodes::Vector{Float64}
    residual::Vector{Float64}
    label::String
    phase::PHS
    interval_residual::Vector{Float64}
end

function assess_solution(
    model::Optimizer;
    q::Integer=10,
    dump::Union{Nothing,AbstractString}=nothing,
)

    residuals = eval_accuracy(model; q)
    solution_error = sum(sum(r.interval_residual) for r in residuals)
    node_residuals = eval_node_residuals(model)
    residual_error = sum(sum(r.interval_residual) for r in node_residuals)
    quad_error = abs(solution_error - residual_error)

    println("Solution Error: ", solution_error)
    println("Residual Error: ", residual_error)
    println("Quadrature Error: ", quad_error)

    if dump !== nothing
        _dump_residuals(dump, node_residuals, residuals)
        println("Wrote per-node residuals to ", dump)
    end

    return (residuals, solution_error, residual_error, quad_error)
end

function _dump_residuals(
    filename::AbstractString,
    node_residuals::Vector{IntervalResidual},
    fine_residuals::Vector{IntervalResidual},
)
    open(filename, "w") do io
        println(io, "source,phase,label,node_idx,t,r2")
        for r in node_residuals
            for (k, (t, r2)) in enumerate(zip(r.nodes, r.residual))
                println(io, "node,", r.phase, ",", r.label, ",", k, ",", t, ",", r2)
            end
        end
        for r in fine_residuals
            for (k, (t, r2)) in enumerate(zip(r.nodes, r.residual))
                println(io, "fine,", r.phase, ",", r.label, ",", k, ",", t, ",", r2)
            end
        end
    end
    return nothing
end

function eval_funcs(
    optimizer::MOI.ModelLike,
    funcs::Vector{<:MOI.AbstractFunction},
    x_vals::Vector{Float64},
)
    var_idxs = MOI.get(optimizer, MOI.ListOfVariableIndices())

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

function eval_funcs(optimizer::MOI.ModelLike, funcs::Vector{<:MOI.AbstractFunction})
    var_idxs = MOI.get(optimizer, MOI.ListOfVariableIndices())
    x_vals   = MOI.get(optimizer, MOI.VariablePrimal(), var_idxs)
    return eval_funcs(optimizer, funcs, Float64.(x_vals))
end

function eval_node_residuals(model::Optimizer)
    results = IntervalResidual[]
    for phase in model.phases
        points = get(model.phase_points, phase, model.default_points)
        ω = points.quad_weights_τ
        τ_alg = points.points_alg_τ
        if length(points.points_dif_τ) == length(points.points_alg_τ)
            τ_dif = points.points_dif_τ
        else
            τ_dif = points.points_dif_τ[1:end-1]
        end

        for (_, (dif_fun, _, _)) in model.dif_cons[phase]
            push!(results, _eval_residual_function(model, phase, dif_fun, τ_dif, ω))
        end
        for (_, (alg_fun, set, _)) in model.alg_cons[phase]
            push!(results, _eval_residual_function(model, phase, alg_fun, set, τ_alg, ω))
        end
    end
    return results
end

function eval_accuracy(model::Optimizer; q::Integer=10)
    if q < 1
        throw(DomainError(q, "Please ensure q ≥ 1."))
    end

    τ_nodes, τ_weights = FGQ.gausslegendre(q)

    results = IntervalResidual[]

    for phase in model.phases
        for (_, (dif_fun, _, _)) in model.dif_cons[phase]
            push!(
                results,
                _eval_residual_function(model, phase, dif_fun, τ_nodes, τ_weights),
            )
        end

        for (_, (alg_fun, set, _)) in model.alg_cons[phase]
            push!(
                results,
                _eval_residual_function(model, phase, alg_fun, set, τ_nodes, τ_weights),
            )
        end
    end

    return results
end

function _eval_residual_function(
    model::Optimizer,
    phase::PHS,
    dif_fun::DIF_FUN,
    τ_nodes::Vector{Float64},
    τ_weights::Vector{Float64},
)
    nodes = Float64[]
    residual = Float64[]
    interval_residual = Float64[]

    for mesh_i in get_points_meshes(model.meshes[phase])
        Δt = 0.5 * (mesh_i.t_b - mesh_i.t_a)
        Σt = 0.5 * (mesh_i.t_b + mesh_i.t_a)
        l1 = 0.0

        for (τ, ω) in zip(τ_nodes, τ_weights)
            t = Σt + Δt * τ
            r = (_evaluate_differential_residual(model, dif_fun, t))^2
            push!(nodes, t)
            push!(residual, r)
            l1 += Δt * ω * r
        end
        push!(interval_residual, l1)
    end

    return IntervalResidual(
        nodes,
        residual,
        _residual_label(model, dif_fun),
        phase,
        interval_residual,
    )
end

function _eval_residual_function(
    model::Optimizer,
    phase::PHS,
    alg_fun::NDF,
    set::MOI.EqualTo{Float64},
    τ_nodes::Vector{Float64},
    τ_weights::Vector{Float64},
)
    nodes = Float64[]
    residual = Float64[]
    interval_residual = Float64[]
    for mesh_i in get_points_meshes(model.meshes[phase])
        Δt = 0.5 * (mesh_i.t_b - mesh_i.t_a)
        Σt = 0.5 * (mesh_i.t_b + mesh_i.t_a)
        l1 = 0.0

        for (τ, ω) in zip(τ_nodes, τ_weights)
            t = Σt + Δt * τ
            value = _evaluate_dynamic_function(model, alg_fun, t)
            r = (_constraint_violation(value, set))^2
            push!(nodes, t)
            push!(residual, r)
            l1 += Δt * ω * r
        end
        push!(interval_residual, l1)
    end

    return IntervalResidual(
        nodes,
        residual,
        _residual_label(model, alg_fun),
        phase,
        interval_residual,
    )
end

function aggregate_residual(residuals::Vector{IntervalResidual}; label::String="aggregate")
    base = residuals[1]
    residual = zeros(length(base.residual))
    interval_residual = zeros(length(base.interval_residual))

    for res in residuals
        residual .+= res.residual
        interval_residual .+= res.interval_residual
    end

    return IntervalResidual(
        base.nodes,
        residual,
        label,
        base.phase,
        interval_residual,
    )
end

function aggregate_residuals(residuals::Vector{IntervalResidual}; label::String="aggregate")
    phases = PHS[]
    for res in residuals
        res.phase in phases || push!(phases, res.phase)
    end

    aggregates = IntervalResidual[]

    for phase in phases
        phase_residuals = [res for res in residuals if res.phase == phase]
        push!(aggregates, aggregate_residual(phase_residuals; label))
    end

    return aggregates
end

function summarize_residual(
    residuals::Vector{IntervalResidual};
    filename::AbstractString="residual_summary.csv",
)
    open(filename, "w") do io
        println(io, "label,phase,mean,variance")

        for res in residuals
            μ = _residual_mean(res.residual)
            σ2 = _residual_variance(res.residual, μ)

            println(io, "$(res.label),$(res.phase),$(μ),$(σ2)")
        end
    end

    return filename
end

function summarize_residual(
    model::Optimizer;
    q::Integer=10,
    filename::AbstractString="residual_summary.csv",
)
    return summarize_residual(eval_accuracy(model; q); filename)
end

function _residual_mean(values::Vector{Float64})
    return sum(values) / length(values)
end

function _residual_variance(values::Vector{Float64}, mean::Float64)
    return sum((value - mean)^2 for value in values) / length(values)
end

function _residual_label(model::Optimizer, dif_fun::DIF_FUN)
    return get(model.dyn_var_names, dif_fun.dyn_var, string(dif_fun.dyn_var))
end

function _residual_label(::Optimizer, ::NDF)
    return "algebraic"
end

function _constraint_violation(value::Float64, set::MOI.EqualTo{Float64})
    return value - set.value
end

function _evaluate_value(
    sol::PiecewiseInterpolant{LagrangeInterpolant},
    t::Float64,
)
    idx, t_eval = _select_piece(sol, t)
    return sol.pieces[idx](t_eval)
end

function _evaluate_derivative(
    sol::PiecewiseInterpolant{LagrangeInterpolant},
    t::Float64,
)
    idx, t_eval = _select_piece(sol, t)
    piece = sol.pieces[idx]
    n = length(piece.points)
    d_values = Vector{Float64}(undef, n)
    for j in 1:n
        wj = piece.weights[j]
        yj = piece.values[j]
        s = 0.0
        for k in 1:n
            k == j && continue
            s += (piece.weights[k] / wj) / (piece.points[j] - piece.points[k]) *
                (piece.values[k] - yj)
        end
        d_values[j] = s
    end
    return _interpolate(piece.points, piece.weights, d_values, t_eval)
end

function _select_piece(sol::PiecewiseInterpolant, t::Float64)
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

    return idx, t_eval
end

function _get_dyn_var_solution(model::Optimizer, phase::PHS, dyn_var::DYN_VAR)
    sol_by_phase = get(model.sol_dyn_vars, phase, nothing)
    if sol_by_phase === nothing || !haskey(sol_by_phase, dyn_var)
        throw(ArgumentError("No solution stored for dynamic variable $(dyn_var)."))
    end
    return sol_by_phase[dyn_var]
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
        sol = _get_dyn_var_solution(model, phase, fun.dyn_fun)
        derivative_value = _evaluate_derivative(sol, t)
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
    elseif arg isa AbstractInterpolant
        return arg(t)
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

function residual_map(model::Optimizer; q::Integer=10)
    return eval_accuracy(model; q)
end

function plot_residual!(args...; kwargs...)
    throw(ArgumentError(
        "Plotting support requires `using Plots` first."
    ))
end

function plot_residual(args...; kwargs...)
    throw(ArgumentError(
        "Plotting support requires `using Plots` first."
    ))
end