function transcribe_sol_derivative!(
    sol_derivatives::SOLS{DOI.Derivative{DYN_VAR}},
    solver::MOI.ModelLike,
    dyn_var_vars::DYN_VAR_VARS,
    derivative::DOI.Derivative{DYN_VAR},
    mesh::AbstractIntervalsMesh,
)
    n_p_dif = get_points_dif_length(mesh)

    vars = dyn_var_vars[derivative.dyn_fun]

    points_meshes = get_points_meshes(mesh)

    sol_derivatives[derivative] = PiecewiseInterpolant([
        LagrangeInterpolant(
            mesh_i.t_a,
            mesh_i.t_b,
            mesh_i.points_dif,
            mesh_i.bary_weights_dif,
            [   
                sum(mesh_i.differentiation[j,k] * MOI.get(
                    solver,
                    MOI.VariablePrimal(),
                    vars[i][k],
                ) for k in 1:n_p_dif)
            for j in 1:n_p_dif],
        ) for (i, mesh_i) in enumerate(points_meshes)
    ])
    return nothing
end