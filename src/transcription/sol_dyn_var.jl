function transcribe_sol_dyn_var!(
    model::Optimizer,
    ::Float64,
    phase::PHS,
    dyn_var::DYN_VAR
)
    transcribe_sol_dyn_var!(
        model.sol_dyn_vars[phase],
        model.inner,
        0.0,
        1.0,
        model.dyn_var_vars,
        dyn_var,
        model.dif_dyn_vars,
        model.meshes[phase]
    )
    return nothing
end

function transcribe_sol_dyn_var!(
    model::Optimizer,
    ::VAR,
    phase::PHS,
    dyn_var::DYN_VAR
)
    phase_initials = OrderedDict{PHS,Float64}()
    Δt = OrderedDict{PHS,Float64}()
    t = model.phase_initials[first(model.phases)].value

    for p in model.phases
        phase_initials[p] = t
        Δt[p] = MOI.get(model.inner, MOI.VariablePrimal(), model.time_vars[p])
        t += Δt[p]
    end

    transcribe_sol_dyn_var!(
        model.sol_dyn_vars[phase],
        model.inner,
        phase_initials[phase],
        Δt[phase],
        model.dyn_var_vars,
        dyn_var,
        model.dif_dyn_vars,
        model.meshes[phase]
    )
    return nothing
end

function transcribe_sol_dyn_var!(
    sol_dyn_vars::SOLS{DYN_VAR},
    solver::MOI.ModelLike,
    t_0::Float64,
    Δt::Float64,
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
                mesh_i.t_a * Δt + t_0,
                mesh_i.t_b * Δt + t_0,
                mesh_i.points_dif .* Δt .+ t_0,
                mesh_i.bary_weights_dif .* (Δt ^ (length(mesh_i.points_dif) - 1)),
                MOI.get(solver, MOI.VariablePrimal(), vars[i]),
            ) for (i, mesh_i) in enumerate(points_meshes)
        ])
    else
        sol_dyn_vars[dyn_var] = PiecewiseInterpolant([
            LagrangeInterpolant(
                mesh_i.t_a * Δt + t_0,
                mesh_i.t_b * Δt + t_0,
                mesh_i.points_alg .* Δt .+ t_0,
                mesh_i.bary_weights_alg .* ((Δt) ^ (length(mesh_i.points_alg) - 1)),
                MOI.get(solver, MOI.VariablePrimal(), vars[i]),
            ) for (i, mesh_i) in enumerate(points_meshes)
        ])
    end
    return nothing
end