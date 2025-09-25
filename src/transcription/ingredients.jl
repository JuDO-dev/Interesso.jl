function transcribe_phase!(model::Optimizer, phase::PHS, ::FixedIntervalsMesh)
    
    if !haskey(model.phase_finals, phase)
        Δt = MOI.add_variable(model.inner)
        MOI.set(model.inner, MOI.VariablePrimalStart(), Δt, 1.0)
        MOI.add_constraint(model.inner, Δt, MOI.GreaterThan(1e-6))
        model.time_vars[phase] = Δt
    else
        model.time_vars[phase] = 1.0
    end
    return nothing
end

function transcribe_phase!(model::Optimizer, phase::PHS, mesh::FlexibleIntervalsMesh)

    if !haskey(model.phase_finals, phase)
        Δt = MOI.add_variable(model.inner)
        MOI.set(model.inner, MOI.VariablePrimalStart(), Δt, 1.0)
        MOI.add_constraint(model.inner, Δt, MOI.GreaterThan(1e-6))
        model.time_vars[phase] = Δt
    else
        model.time_vars[phase] = 1.0
    end

    n_h = get_intervals_length(mesh)
    t_0 = mesh.fixed.points_meshes[1].t_a
    t_f = mesh.fixed.points_meshes[end].t_b

    flex_vars = MOI.add_variables(model.inner, n_h - 1)

    for i in eachindex(flex_vars)
        MOI.set(
            model.inner,
            MOI.VariablePrimalStart(),
            flex_vars[i],
            mesh.fixed.points_meshes[i].t_b,
        )
    end

    MOI.add_constraint(
        model.inner, 
        1.0 * flex_vars[1] - t_0,
        MOI.Interval(mesh.Δt_min, mesh.Δt_max),
    )

    for i in 2:(n_h - 1)
        MOI.add_constraint(
            model.inner,
            1.0 * flex_vars[i] - 1.0 * flex_vars[i-1],
            MOI.Interval(mesh.Δt_min, mesh.Δt_max),
        )
    end

    MOI.add_constraint(
        model.inner,
        t_f - 1.0 * flex_vars[end],
        MOI.Interval(mesh.Δt_min, mesh.Δt_max)
    )

    model.phase_vars[phase] = flex_vars

    return nothing
end

function transcribe_dyn_vars!(model::Optimizer, phase::PHS, mesh::AbstractIntervalsMesh)
    
    n_h = get_intervals_length(mesh)
    n_p_dif = get_points_dif_length(mesh)
    n_p_alg = get_points_alg_length(mesh)

    for dyn_var in model.dyn_vars[phase]
        if dyn_var in model.dif_dyn_vars

            vars = [
                [MOI.add_variable(model.inner) for _ in 1:n_p_dif] for _ in 1:n_h
            ]

            for i in 2:n_h
                MOI.add_constraint(
                    model.inner,
                    1.0 * last(vars[i-1]) - 1.0 * first(vars[i]),
                    MOI.EqualTo(0.0),
                )
            end     

            model.dyn_var_vars[dyn_var] = vars
        else
            model.dyn_var_vars[dyn_var] = [
                [MOI.add_variable(model.inner) for _ in 1:n_p_alg] for _ in 1:n_h
            ]
        end
    end
    return nothing
end

function transcribe_dyn_var_starts!(
    model::Optimizer,
    phase::PHS,
    mesh::AbstractIntervalsMesh,
)
    n_h = get_intervals_length(mesh)
    n_p_dif = get_points_dif_length(mesh)
    n_p_alg = get_points_alg_length(mesh)

    points_meshes = get_points_meshes(mesh)

    for (dyn_var, start) in model.start_dyn_vars[phase]
        
        vars = model.dyn_var_vars[dyn_var]
        
        if dyn_var in model.dif_dyn_vars
            for i in 1:n_h
                for j in 1:n_p_dif
                    MOI.set(
                        model.inner,
                        MOI.VariablePrimalStart(),
                        vars[i][j],
                        start(points_meshes[i].points_dif[j]),
                    )
                end
            end
        else
            for i in 1:n_h
                for j in 1:n_p_alg
                    MOI.set(
                        model.inner,
                        MOI.VariablePrimalStart(),
                        vars[i][j],
                        start(points_meshes[i].points_alg[j]),
                    )
                end
            end
        end
    end
    return nothing
end

function transcribe_bounds!(
    model::Optimizer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM,BM<:ExactBoundsMesh}

    n_h = get_intervals_length(mesh)
    n_p_dif = get_points_dif_length(mesh)
    n_p_alg = get_points_alg_length(mesh)

    for (dyn_var, set) in model.dyn_var_bounds[phase]
        
        vars = model.dyn_var_vars[dyn_var]

        if dyn_var in model.dif_dyn_vars
            for i in 1:n_h
                for j in 1:n_p_dif
                    MOI.add_constraint(model.inner, vars[i][j], set)
                end
            end
        else
            for i in 1:n_h
                for j in 1:n_p_alg
                    MOI.add_constraint(model.inner, vars[i][j], set)
                end
            end
        end
    end
    return nothing
end

function transcribe_bounds!(
    model::Optimizer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM,BM<:SampledBoundsMesh}

    n_h = get_intervals_length(mesh)
    n_p_dif = get_points_dif_length(mesh)
    n_p_alg = get_points_alg_length(mesh)
    n_p_samp = get_points_samp_length(mesh)

    bounds_mesh = get_bounds_mesh(mesh)

    for (dyn_var, set) in model.dyn_var_bounds[phase]

        vars = model.dyn_var_vars[dyn_var]

        if dyn_var in model.dif_dyn_vars
            for i in 1:n_h
                for j in 1:n_p_samp
                    MOI.add_constraint(
                        model.inner,
                        sum(bounds_mesh.sampled_dif[j,k] * vars[i][k] for k in 1:n_p_dif),
                        set,
                    )
                end
            end
        else
            for i in 1:n_h
                for j in 1:n_p_samp
                    MOI.add_constraint(
                        model.inner,
                        sum(bounds_mesh.sampled_alg[j,k] * vars[i][k] for k in 1:n_p_alg),
                        set,
                    )
                end
            end
        end
    end
    return nothing
end

function transcribe_bounds!(
    model::Optimizer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM,BM<:BernsteinBoundsMesh}

    n_h = get_intervals_length(mesh)
    n_p_dif = get_points_dif_length(mesh)
    n_p_alg = get_points_alg_length(mesh)

    bounds_mesh = get_bounds_mesh(mesh)

    for (dyn_var, set) in model.dyn_var_bounds[phase]

        vars = model.dyn_var_vars[dyn_var]

        if dyn_var in model.dif_dyn_vars
            for i in 1:n_h
                for j in 1:n_p_dif
                    MOI.add_constraint(
                        model.inner,
                        sum(bounds_mesh.bernstein_dif[j,k] * vars[i][k] for k in 1:n_p_dif),
                        set,
                    )
                end
            end
        else
            for i in 1:n_h
                for j in 1:n_p_alg
                    MOI.add_constraint(
                        model.inner,
                        sum(bounds_mesh.bernstein_alg[j,k] * vars[i][k] for k in 1:n_p_alg),
                        set,
                    )
                end
            end
        end
    end
    return nothing
end

function transcribe_dif_cons!(
    model::Optimizer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:CollocationMesh,BM}

    dif_cons = model.dif_cons[phase]

    n_h = get_intervals_length(mesh)
    n_p_alg = get_points_alg_length(mesh)

    for (dif_fun, set) in values(dif_cons)
        for i in 1:n_h
            for j in 1:n_p_alg
                MOI.add_constraint(
                    model.inner,
                    transcribe_dyn_fun(dif_fun, i, j, model.phase_vars, model.time_vars[phase],
                        model.dyn_var_vars, model.dif_dyn_vars, mesh
                    ),
                    set,
                )
            end
        end
    end
    return nothing
end    

function transcribe_alg_cons!(
    model::Optimizer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:CollocationMesh,BM}

    alg_cons = model.alg_cons[phase]

    n_h = get_intervals_length(mesh)
    n_p_alg = get_points_alg_length(mesh)

    for (alg_fun, set) in values(alg_cons)
        for i in 1:n_h
            for q in 1:n_p_alg
                MOI.add_constraint(
                    model.inner,
                    transcribe_dyn_fun(alg_fun, i, q, model.phase_vars, model.time_vars[phase],
                        model.dyn_var_vars, model.dif_dyn_vars, mesh
                    ), 
                    set,
                )
            end
        end
    end
    return nothing
end

function transcribe_dif_cons!(
    model::Optimizer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:IntResidualMesh,BM}

    n_h = get_intervals_length(mesh)
    for i = 1:n_h
        grad_res_funcs = transcribe_grad_dyn_least_square(model, i, phase, mesh)
        for f in grad_res_funcs
            MOI.add_constraint(
                model.inner,
                f,
                MOI.Interval(-1e-4, 1e-4),
            )
            push!(model.dif_res_funcs, f)
        end
    end
    return nothing
end

function transcribe_alg_cons!(
    model::Optimizer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:IntResidualMesh,BM}

    n_h = get_intervals_length(mesh)
    for i = 1:n_h
        f = transcribe_dyn_least_square(model, i, phase, mesh)
        MOI.add_constraint(
            model.inner,
            f,
            MOI.LessThan(1e-2),
        )
        push!(model.res_funcs, f)
    end
    return nothing
end

function transcribe_initials!(model::Optimizer, phase::PHS, mesh::AbstractIntervalsMesh)

    for (dyn_var, set) in model.dyn_var_initials[phase]
        MOI.add_constraint(
            model.inner,
            transcribe_dyn_var_initial(dyn_var, model.dyn_var_vars, mesh),
            set,
        )
    end
    return nothing
end

function transcribe_finals!(model::Optimizer, phase::PHS, mesh::AbstractIntervalsMesh)

    for (dyn_var, set) in model.dyn_var_finals[phase]
        MOI.add_constraint(
            model.inner,
            transcribe_dyn_var_final(dyn_var, model.dyn_var_vars, mesh),
            set,
        )
    end
    return nothing
end

function transcribe_linkages!(
    model::Optimizer,
    meshes::MESHES,
)
    for (_, (linkage, set)) in model.linkages
        MOI.add_constraint(
            model.inner,
            transcribe_dyn_var_final(
                linkage.dyn_fun_final,
                model.dyn_var_vars,
                meshes[DOI.phase_index(linkage.dyn_fun_final)],
            ) - transcribe_dyn_var_initial(
                linkage.dyn_fun_initial,
                model.dyn_var_vars,
                meshes[DOI.phase_index(linkage.dyn_fun_initial)],
            ),
            set,
        )
    end
    return nothing
end