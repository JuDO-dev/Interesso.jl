function transcribe_sol_derivative!(
    model::Optimizer,
    ::Float64,
    phase::PHS,
    derivative::DOI.Derivative{DYN_VAR},
)
    transcribe_sol_derivative!(
        model.sol_derivatives[phase],
        model.inner,
        0.0,
        1.0,
        model.dyn_var_vars,
        derivative,
        model.meshes[phase],
    )
    return nothing
end

function transcribe_sol_derivative!(
    model::Optimizer,
    ::VAR,
    phase::PHS,
    derivative::DOI.Derivative{DYN_VAR},
)
    phase_initials = OrderedDict{PHS,Float64}()
    Δt = OrderedDict{PHS,Float64}()
    t = model.phase_initials[first(model.phases)]

    for p in model.phases
        phase_initials[p] = t
        Δt[p] = MOI.get(model.inner, MOI.VariablePrimal(), model.time_vars[p])
        t += Δt[p]
    end

    transcribe_sol_derivative!(
        model.sol_derivatives[phase],
        model.inner,
        phase_initials[phase],
        Δt[phase],
        model.dyn_var_vars,
        derivative,
        model.meshes[phase],
    )
    return nothing
end

function transcribe_sol_derivative!(
    sol_derivatives::SOLS{DOI.Derivative{DYN_VAR}},
    solver::MOI.ModelLike,
    t_0::Float64,
    Δt::Float64,
    dyn_var_vars::DYN_VAR_VARS,
    derivative::DOI.Derivative{DYN_VAR},
    mesh::AbstractIntervalsMesh,
)
    n_p_dif = get_points_dif_length(mesh)

    vars = dyn_var_vars[derivative.dyn_fun]

    points_meshes = get_points_meshes(mesh)

    sol_derivatives[derivative] = PiecewiseInterpolant([
        LagrangeInterpolant(
            mesh_i.t_a * Δt + t_0,
            mesh_i.t_b * Δt + t_0,
            mesh_i.points_dif .* Δt .+ t_0,
            mesh_i.bary_weights_dif .* (Δt ^ (length(mesh_i.points_dif) - 1)),
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