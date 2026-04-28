struct IntervalResidual
    phase::PHS
    interval::Int
    t_a::Float64
    t_b::Float64
    residual::Float64
    point_mesh::Vector{Float64}
    point_res::Vector{Float64}
end

function assess_solution(model::Optimizer; q::Integer=10)

    residuals = eval_accuracy(model; q)
    solution_error = sum(r.residual for r in residuals)
    residual_error = sum(abs, eval_funcs(model.inner, model.res_funcs))
    quad_error = abs(solution_error - residual_error)

    println("Solution Error: ", solution_error)
    println("Residual Error: ", residual_error)
    println("Quadrature Error: ", quad_error)

    return (residuals, solution_error, residual_error, quad_error)
end

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

    results = IntervalResidual[]

    for phase in model.phases
        mesh = get(model.meshes, phase, nothing)
        mesh === nothing && continue

        dif_cons = collect(values(model.dif_cons[phase]))
        alg_cons = collect(values(model.alg_cons[phase]))
        (isempty(dif_cons) && isempty(alg_cons)) && continue

        for (idx, mesh_i) in enumerate(get_points_meshes(mesh))
            Δt = 0.5 * (mesh_i.t_b - mesh_i.t_a)
            Σt = 0.5 * (mesh_i.t_b + mesh_i.t_a)

            interval_res = 0.0
            point_mesh = Float64[]
            point_res = Float64[]

            for (τ, ω) in zip(τ_nodes, τ_weights)
                t = Σt + Δt * τ
                weight = Δt * ω

                R = 0.0

                for (dif_fun, _) in dif_cons
                    residual = _evaluate_differential_residual(model, dif_fun, t)
                    R += abs(residual)
                end

                for (alg_fun, set) in alg_cons
                    value = _evaluate_dynamic_function(model, alg_fun, t)
                    violation = _constraint_violation(value, set)
                    R += abs(violation)
                end

                interval_res += weight * R
                push!(point_mesh, t)
                push!(point_res, R)
            end

            push!(results, IntervalResidual(phase, idx, mesh_i.t_a, mesh_i.t_b, interval_res, point_mesh, point_res))
        end
    end

    return results
end

function _constraint_violation(value::Float64, set::MOI.EqualTo{Float64})
    return value - set.value
end

function _constraint_violation(value::Float64, set::MOI.LessThan{Float64})
    return max(0.0, value - set.upper)
end

function _constraint_violation(value::Float64, set::MOI.GreaterThan{Float64})
    return max(0.0, set.lower - value)
end

function _constraint_violation(value::Float64, set::MOI.Interval{Float64})
    if value > set.upper
        return value - set.upper
    elseif value < set.lower
        return set.lower - value
    else
        return 0.0
    end
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
    results = eval_accuracy(model; q)
    return (
        domain = reduce(vcat, r.point_mesh for r in results),
        l1_residual = reduce(vcat, r.point_res  for r in results),
        interval = reduce(vcat, fill(r.interval, length(r.point_mesh)) for r in results),
        q = q,
    )
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

function list_constraint_violations(
    model::Optimizer;
    q::Integer=10,
    path::AbstractString="constraint_violations.csv",
    map_path::AbstractString="constraint_violation_map.csv",
    min_magnitude::Float64=0.0,
    top_n::Union{Nothing,Int}=nothing,
)
    if q < 1
        throw(DomainError(q, "Please ensure q ≥ 1."))
    end
    if min_magnitude < 0.0
        throw(DomainError(min_magnitude, "Please ensure min_magnitude ≥ 0.0."))
    end
    if top_n !== nothing && top_n < 1
        throw(DomainError(top_n, "Please ensure top_n ≥ 1 or nothing."))
    end

    τ_nodes, τ_weights = FGQ.gausslegendre(q)

    rows = NamedTuple[]
    map_rows = NamedTuple[]

    for phase in model.phases
        mesh = get(model.meshes, phase, nothing)
        mesh === nothing && continue

        dif_cons = collect(pairs(model.dif_cons[phase]))
        alg_cons = collect(pairs(model.alg_cons[phase]))

        dif_ids = Dict{Any,Int}()
        alg_ids = Dict{Any,Int}()

        for (k, (_, (dif_fun, _))) in enumerate(dif_cons)
            dif_ids[dif_fun] = k
            push!(map_rows, (
                phase=phase,
                constraint_type="dif",
                constraint_id=k,
                label=_compact_constraint_label(dif_fun),
            ))
        end

        for (k, (_, (alg_fun, _))) in enumerate(alg_cons)
            alg_ids[alg_fun] = k
            push!(map_rows, (
                phase=phase,
                constraint_type="alg",
                constraint_id=k,
                label=_compact_constraint_label(alg_fun),
            ))
        end

        for (i, mesh_i) in enumerate(get_points_meshes(mesh))
            Δt = 0.5 * (mesh_i.t_b - mesh_i.t_a)
            Σt = 0.5 * (mesh_i.t_b + mesh_i.t_a)

            for (q_idx, (τ, ω)) in enumerate(zip(τ_nodes, τ_weights))
                t = Σt + Δt * τ
                w = Δt * ω

                for (_, (dif_fun, _)) in dif_cons
                    value = _evaluate_differential_residual(model, dif_fun, t)
                    magnitude = abs(value)
                    magnitude < min_magnitude && continue

                    push!(rows, (
                        phase=phase,
                        constraint_type="dif",
                        constraint_id=dif_ids[dif_fun],
                        i=i,
                        q=q_idx,
                        t=t,
                        magnitude=magnitude,
                        l1_term=w * magnitude,
                    ))
                end

                for (_, (alg_fun, _)) in alg_cons
                    value = _evaluate_dynamic_function(model, alg_fun, t)
                    magnitude = abs(value)
                    magnitude < min_magnitude && continue

                    push!(rows, (
                        phase=phase,
                        constraint_type="alg",
                        constraint_id=alg_ids[alg_fun],
                        i=i,
                        q=q_idx,
                        t=t,
                        magnitude=magnitude,
                        l1_term=w * magnitude,
                    ))
                end
            end
        end
    end

    sort!(rows; by=row -> row.l1_term, rev=true)

    if top_n !== nothing && length(rows) > top_n
        resize!(rows, top_n)
    end

    open(path, "w") do io
        println(io, "phase,constraint_type,constraint_id,i,q,t,magnitude,l1_term")
        for row in rows
            println(
                io,
                string(
                    row.phase, ",",
                    row.constraint_type, ",",
                    row.constraint_id, ",",
                    row.i, ",",
                    row.q, ",",
                    row.t, ",",
                    row.magnitude, ",",
                    row.l1_term,
                ),
            )
        end
    end

    open(map_path, "w") do io
        println(io, "phase,constraint_type,constraint_id,label")
        for row in map_rows
            println(
                io,
                string(
                    row.phase, ",",
                    row.constraint_type, ",",
                    row.constraint_id, ",",
                    _csv_escape(row.label),
                ),
            )
        end
    end

    return rows
end

function _compact_constraint_label(dif_fun::DOI.ExplicitDifferentialFunction)
    return string("d/dt(", _compact_string(dif_fun.dyn_var), ") - ", _compact_string(dif_fun.dyn_fun))
end

function _compact_constraint_label(fun)
    return _compact_string(fun)
end

function _compact_string(x)
    s = repr(x)
    s = replace(s, '\n' => ' ', '\r' => ' ', '\t' => ' ')
    s = replace(s, r"\s+" => " ")
    return strip(s)
end

function _csv_escape(s::AbstractString)
    if occursin(',', s) || occursin('"', s)
        return "\"" * replace(s, "\"" => "\"\"") * "\""
    end
    return s
end