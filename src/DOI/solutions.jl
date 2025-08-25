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