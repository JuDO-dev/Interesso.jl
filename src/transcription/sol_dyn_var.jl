function transcribe_sol_dyn_var!(
    sol_dyn_vars::SOLS{DYN_VAR},
    solver::MOI.ModelLike,
    dyn_var_vars::DYN_VAR_VARS,
    dyn_var::DYN_VAR,
    dif_dyn_vars::OrderedSet{DYN_VAR},
    mesh::AbstractIntervalsMesh,
)
    vars = dyn_var_vars[dyn_var]

    points_meshes = get_points_meshes(mesh)

    if dyn_var in dif_dyn_vars
        sol_dyn_vars[dyn_var] = PiecewiseInterpolant([
            LagrangeInterpolant(
                mesh_i.t_a,
                mesh_i.t_b,
                mesh_i.points_dif,
                mesh_i.bary_weights_dif,
                MOI.get(solver, MOI.VariablePrimal(), vars[i]),
            ) for (i, mesh_i) in enumerate(points_meshes)
        ])
    else
        sol_dyn_vars[dyn_var] = PiecewiseInterpolant([
            LagrangeInterpolant(
                mesh_i.t_a,
                mesh_i.t_b,
                mesh_i.points_alg,
                mesh_i.bary_weights_alg,
                MOI.get(solver, MOI.VariablePrimal(), vars[i]),
            ) for (i, mesh_i) in enumerate(points_meshes)
        ])
    end
    return nothing
end