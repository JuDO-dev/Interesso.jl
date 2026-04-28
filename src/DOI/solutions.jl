function update_mesh!(
    ::FixedIntervalsMesh,
    ::MOI.ModelLike,
    ::PHS_VARS,
    ::PHS,
    ::AbstractPoints,
)
    return nothing
end

function update_mesh!(
    mesh::FlexibleIntervalsMesh,
    solver::MOI.ModelLike,
    phase_vars::PHS_VARS,
    phase::PHS,
    points::AbstractPoints,
)
    n_h = get_intervals_length(mesh)

    flexed_vars = MOI.get(
        solver,
        MOI.VariablePrimal(),
        phase_vars[phase],
    )

    flexed_points = vcat(
        mesh.fixed.points_meshes[1].t_a,
        flexed_vars,
        mesh.fixed.points_meshes[end].t_b,
    )
        
    for i in 1:n_h
        mesh.fixed.points_meshes[i] = build_points_mesh(
            points,
            flexed_points[i],
            flexed_points[i+1],
        )
    end
    return nothing
end


function MOI.get(model::Optimizer, attr::MOI.VariablePrimal, var::VAR)
    return MOI.get(model.inner, attr, var)
end

function MOI.get(model::Optimizer, ::DOI.DynamicVariableSolution, dyn_var::DYN_VAR)

    phase = DOI.phase_index(dyn_var)

    if !haskey(model.sol_dyn_vars[phase], dyn_var)
        transcribe_sol_dyn_var!(model, model.time_vars[phase], phase, dyn_var)
    end
    return model.sol_dyn_vars[phase][dyn_var]
end

function save_solutions!(model::Optimizer)

    for phase in model.phases
        
        time_var = model.time_vars[phase]

        for dyn_var in model.dyn_vars[phase]
            transcribe_sol_dyn_var!(model, time_var, phase, dyn_var)
            if dyn_var in model.dif_dyn_vars
                transcribe_sol_derivative!(
                    model,
                    time_var,
                    phase,
                    DOI.Derivative(dyn_var),
                )
            end
        end
    end

    return nothing
end

function get_solutions(model::Optimizer)
    solutions = WSS{DOI.AbstractDynamicSolution}()
    for phase in model.phases
        sols = WS{DOI.AbstractDynamicSolution}()
        for dyn_var in model.dyn_vars[phase]
            name = get(model.dyn_var_names, dyn_var, nothing)
            if name !== nothing
                sol = MOI.get(model, DOI.DynamicVariableSolution(), dyn_var)
                sols[name] = sol
            end
        end
        solutions[phase] = sols
    end
    return solutions
end

function normalize_solutions!(model::Optimizer)
    for phase in model.phases
        for dyn_var in model.dyn_vars[phase]
            transcribe_sol_dyn_var!(model, 1.0, phase, dyn_var)
            if dyn_var in model.dif_dyn_vars
                transcribe_sol_derivative!(
                    model,
                    1.0,
                    phase,
                    DOI.Derivative(dyn_var),
                )
            end
        end
    end
    return nothing
end

function warmstart!(
    model::Optimizer,
    starts::WSS{T}
) where {T<:DOI.AbstractDynamicSolution}
    for (dyn_var, name) in model.dyn_var_names
        phase = DOI.phase_index(dyn_var)
        phase_dict = get(starts, phase, nothing)
        if phase_dict === nothing
            continue
        end
        sol = get(phase_dict, name, nothing)
        if sol !== nothing
            model.start_dyn_vars[phase][dyn_var] = sol
        end
    end
    return nothing
end

function get_primal(model::Interesso.Optimizer)
    vars = MOI.get(model.inner, MOI.ListOfVariableIndices())
    sort!(vars; by = v -> v.value)
    return MOI.get(model.inner, MOI.VariablePrimal(), vars)
end

const _NLPBLOCK_DUAL_KEY = (MOI.NLPBlock, Float64)

get_nlpblock_dual(model::Interesso.Optimizer) = Float64.(MOI.get(model.inner, MOI.NLPBlockDual()))

function get_dual(model::Interesso.Optimizer)

    dual = Dict{Tuple{DataType,DataType}, Vector{Float64}}()

    for (F,S) in MOI.get(model.inner, MOI.ListOfConstraintTypesPresent())
        cons = MOI.get(model.inner, MOI.ListOfConstraintIndices{F,S}())
        sort!(cons; by = c -> c.value)
        dual[(F,S)] = Float64.(MOI.get(model.inner, MOI.ConstraintDual(), cons))
    end

    dual[_NLPBLOCK_DUAL_KEY] = Float64.(MOI.get(model.inner, MOI.NLPBlockDual()))

    return dual
end

get_primal_dual(model::Interesso.Optimizer) = (get_primal(model), get_dual(model))

function set_primal_start!(model::Interesso.Optimizer, x0::AbstractVector{T}) where {T<:Real}
    vars = MOI.get(model.inner, MOI.ListOfVariableIndices())
    sort!(vars; by = v -> v.value)

    @assert length(vars) == length(x0) "Primal length mismatch."
    for (v, xv) in zip(vars, x0)
        MOI.set(model.inner, MOI.VariablePrimalStart(), v, Float64(xv))
    end
    return nothing
end

function set_dual_start!(model::Interesso.Optimizer, dual::Dict{Tuple{DataType,DataType},Vector{Float64}})
    if haskey(dual, _NLPBLOCK_DUAL_KEY)
        MOI.set(model.inner, MOI.NLPBlockDualStart(), dual[_NLPBLOCK_DUAL_KEY])
    end
    for (F,S) in MOI.get(model.inner, MOI.ListOfConstraintTypesPresent())
        vals = get(dual, (F,S), nothing)
        vals === nothing && continue

        cons = MOI.get(model.inner, MOI.ListOfConstraintIndices{F,S}())
        sort!(cons; by = c -> c.value)

        @assert length(cons) == length(vals) "Dual length mismatch for ($F,$S)."
        for (c, μ0) in zip(cons, vals)
            MOI.set(model.inner, MOI.ConstraintDualStart(), c, μ0)
        end
    end
    return nothing
end

function warmstart!(
    model::Optimizer;
    primal::Union{Nothing,Vector{Float64}}=nothing,
    dual::Union{Nothing,Dict{Tuple{DataType,DataType},Vector{Float64}}}=nothing
)
    if primal !== nothing
        set_primal_start!(model, primal)
    end
    if dual !== nothing
        set_dual_start!(model, dual)
    end
    return nothing
end