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
    warm_start = Dict{String, DOI.AbstractDynamicSolution}()

    for phase in model.phases
        for dyn_var in model.dyn_vars[phase]
            name = get(model.dyn_var_names, dyn_var, nothing)
            if name !== nothing
                if name == "t"
                    error("Avoid setting variables as name t.")
                end
                sol = MOI.get(model, DOI.DynamicVariableSolution(), dyn_var)
                warm_start[name] = sol
            end
        end

        # if model.time_vars[phase] isa VAR
        #     val = MOI.get(model.inner, MOI.VariablePrimal(), model.time_vars[phase])
        #     push!(sol_t, val)
        #     warm_start["t"] = sol_t
        # end

    end
    return warm_start
end

function warmstart!(
    model::Optimizer,
    starts::AbstractDict{String, T},
) where {T<:DOI.AbstractDynamicSolution}

    for (dyn_var, name) in model.dyn_var_names
        sol = get(starts, name, nothing)
        if !isnothing(sol)
            phase = DOI.phase_index(dyn_var)
            model.start_dyn_vars[phase][dyn_var] = starts[name]
        end
    end
    return nothing
end