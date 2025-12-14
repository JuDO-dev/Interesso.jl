function transcribe_phase!(model::Optimizer, phase::PHS, ::FixedIntervalsMesh)

    function _time(phase::PHS)
        terms = MOI.ScalarAffineTerm{Float64}[]
        constant = 0.0

        for p in model.phases
            t = model.time_vars[p]
            if t isa MOI.VariableIndex
                push!(terms, MOI.ScalarAffineTerm(1.0, t))
            else
                constant += (model.phase_finals[phase].value - model.phase_initials[phase].value)
            end
            p == phase && break
        end

        return MOI.ScalarAffineFunction(terms, constant)
    end

    final = get(model.phase_finals, phase, nothing)

    if final isa MOI.EqualTo
        model.time_vars[phase] = 1.0
        return nothing
    end

    Δt = MOI.add_variable(model.inner)
    model.time_vars[phase] = Δt
    MOI.set(model.inner, MOI.VariablePrimalStart(), Δt, 1.0)

    if final === nothing
        MOI.add_constraint(model.inner, Δt, MOI.GreaterThan(0.0))
    elseif final isa MOI.LessThan
        MOI.add_constraint(model.inner, _time(phase), MOI.Interval(0.0, final.upper))
    elseif final isa MOI.GreaterThan
        MOI.add_constraint(model.inner, _time(phase), MOI.GreaterThan(max(0.0, final.lower)))
    elseif final isa MOI.Interval
        MOI.add_constraint(model.inner, _time(phase), final)
    end
    return nothing
end

function transcribe_phase!(model::Optimizer, phase::PHS, mesh::FlexibleIntervalsMesh)

    function _time(phase::PHS)
        terms = MOI.ScalarAffineTerm{Float64}[]
        constant = 0.0

        for p in model.phases
            t = model.time_vars[p]
            if t isa MOI.VariableIndex
                push!(terms, MOI.ScalarAffineTerm(1.0, t))
            else
                constant += (model.phase_finals[phase].value - model.phase_initials[phase].value)
            end
            p == phase && break
        end

        return MOI.ScalarAffineFunction(terms, constant)
    end

    final = get(model.phase_finals, phase, nothing)

    if final isa MOI.EqualTo
        model.time_vars[phase] = 1.0
        return nothing
    end

    Δt = MOI.add_variable(model.inner)
    model.time_vars[phase] = Δt
    MOI.set(model.inner, MOI.VariablePrimalStart(), Δt, 1.0)

    if final === nothing
        MOI.add_constraint(model.inner, Δt, MOI.GreaterThan(0.0))
    elseif final isa MOI.LessThan
        MOI.add_constraint(model.inner, _time(phase), MOI.Interval(0.0, final.upper))
    elseif final isa MOI.GreaterThan
        MOI.add_constraint(model.inner, _time(phase), MOI.GreaterThan(max(0.0, final.lower)))
    elseif final isa MOI.Interval
        MOI.add_constraint(model.inner, _time(phase), final)
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

function transcribe_dyn_vars!(
    model::Optimizer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM<:AbstractRadauMesh,MM,BM}
    n_h = get_intervals_length(mesh)
    n_p_dif = get_points_dif_length(mesh)
    n_p_alg = get_points_alg_length(mesh)

    for dyn_var in model.dyn_vars[phase]
        if dyn_var in model.dif_dyn_vars

            vars = Vector{Vector{VAR}}(undef, n_h)

            vars[1] = [MOI.add_variable(model.inner) for _ in 1:n_p_dif]
            for i in 2:n_h
                vars[i] = Vector{VAR}(undef, n_p_dif)
                vars[i][1] = vars[i-1][end]
                for j in 2:(n_p_dif)
                    vars[i][j] = MOI.add_variable(model.inner)
                end
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

function transcribe_dyn_vars!(
    model::Optimizer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM<:AbstractLobattoMesh,MM,BM}
    
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
) where {PM<:AbstractRadauMesh,MM,BM<:ExactBoundsMesh}

    n_h = get_intervals_length(mesh)
    n_p_dif = get_points_dif_length(mesh)
    n_p_alg = get_points_alg_length(mesh)

    for (dyn_var, set) in model.dyn_var_bounds[phase]
        
        vars = model.dyn_var_vars[dyn_var]

        if dyn_var in model.dif_dyn_vars
            MOI.add_constraint(model.inner, vars[1][1], set)
            for i in 1:n_h
                for j in 2:n_p_dif
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
) where {PM<:AbstractLobattoMesh,MM,BM<:ExactBoundsMesh}

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

# Collocation, dynamic equations
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
            for q in 1:n_p_alg

                dif_con = transcribe_dyn_fun(
                    dif_fun, i, q, model.phase_vars, model.time_vars[phase],
                    model.dyn_var_vars, model.dif_dyn_vars, mesh
                )

                MOI.add_constraint(model.inner, dif_con, set)
                push!(model.res_funcs, dif_con)
            end
        end
    end
    return nothing
end    

# Collocation, algebraic equations
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

                alg_con = transcribe_dyn_fun(
                    alg_fun, i, q, model.phase_vars, model.time_vars[phase],
                    model.dyn_var_vars, model.dif_dyn_vars, mesh
                )

                MOI.add_constraint(model.inner, alg_con, set)
                push!(model.res_funcs, alg_con)
            end
        end
    end
    return nothing
end

# Integrated Residual, transcription of differentiation of residuals
function transcribe_dif_cons!(
    ::Optimizer,
    ::PHS,
    ::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:Union{DAIRMesh,QPMMesh},BM}
    return nothing
end

function transcribe_dif_cons!(
    model::Optimizer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:Union{SAIRMesh,SAPMMesh},BM}

    grad_res_funcs = transcribe_grad_dyn_least_square(model, phase, mesh)

    for f in grad_res_funcs
        MOI.add_constraint(
            model.inner,
            f,
            # MOI.Interval(-1e-4, 1e-4),
            MOI.EqualTo(0.0)
        )
        push!(model.dif_res_funcs, f)
    end

    return nothing
end

# Integrated Residual, transcription of residuals
function transcribe_alg_cons!(
    model::Optimizer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:Union{DAIRMesh,SAIRMesh},BM}

    n_h = get_intervals_length(mesh)
    ϵ = 1e-4

    scale = length(model.dif_cons[phase]) + length(model.alg_cons[phase])
    ϵ *= scale
    for i = 1:n_h
        f = transcribe_dyn_least_square(model, i, phase, mesh)
        MOI.add_constraint(
            model.inner,
            f,
            MOI.LessThan(ϵ),
        )
        push!(model.res_funcs, f)
    end

    return nothing
end

function transcribe_alg_cons!(
    ::Optimizer,
    ::PHS,
    ::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:Union{QPMMesh,SAPMMesh},BM}
    return nothing
end

# Boundary constraints
function transcribe_initials!(model::Optimizer, meshes::MESHES)

    phase = model.phases[1]
    
    for (dyn_var, set) in model.dyn_var_initials[phase]
        MOI.add_constraint(
            model.inner,
            transcribe_dyn_var_initial(dyn_var, model.dyn_var_vars, meshes[phase]),
            set,
        )
    end
    return nothing
end

function transcribe_finals!(model::Optimizer, meshes::MESHES)

    phase = model.phases[end]

    for (dyn_var, set) in model.dyn_var_finals[phase]
        MOI.add_constraint(
            model.inner,
            transcribe_dyn_var_final(dyn_var, model.dyn_var_vars, meshes[phase]),
            set,
        )
    end
    return nothing
end

function transcribe_bou_cons!(model::Optimizer, meshes::MESHES)

    for (_, (fun, set)) in model.bou_cons
        MOI.add_constraint(
            model.inner,
            transcribe_bou_fun(fun, model, meshes),
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