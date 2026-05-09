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
    else
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
    else
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
        1.0 * flex_vars[1],
        MOI.Interval(mesh.Δt_min + t_0, mesh.Δt_max + t_0),
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
        - 1.0 * flex_vars[end],
        MOI.Interval(mesh.Δt_min - t_f, mesh.Δt_max - t_f)
    )

    model.phase_vars[phase] = flex_vars

    return nothing
end

function transcribe_dyn_vars!(
    model::Optimizer,
    phase::PHS,
    mesh::AbstractIntervalsMesh,
)
    n_h = get_intervals_length(mesh)
    n_p_dif = get_points_dif_length(mesh)
    n_p_alg = get_points_alg_length(mesh)

    for dyn_var in model.dyn_vars[phase]
        if dyn_var in model.dif_dyn_vars
            vars = [Vector{VAR}(undef, n_p_dif) for _ in 1:n_h]
            model.dyn_var_vars[dyn_var] = vars
        else
            vars = [Vector{VAR}(undef, n_p_alg) for _ in 1:n_h]
            model.dyn_var_vars[dyn_var] = vars
        end
    end
    return nothing
end

function transcribe_dyn_vars!(
    model::Optimizer,
    i::Integer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM<:AbstractRadauMesh,MM,BM}
    n_p_dif = get_points_dif_length(mesh)
    n_p_alg = get_points_alg_length(mesh)

    for j in 1:n_p_dif
        for dyn_var in model.dyn_vars[phase]
            vars = model.dyn_var_vars[dyn_var]
            if dyn_var in model.dif_dyn_vars
                if i > 1 && j == 1
                    vars[i][j] = vars[i-1][end]
                else
                    vars[i][j] = MOI.add_variable(model.inner)
                end
            elseif j <= n_p_alg
                vars[i][j] = MOI.add_variable(model.inner)
            end
        end
    end

    # add initial constraints
    if i == 1
        transcribe_initials!(model, phase, mesh)
    end
    # add final constraints
    if i == get_intervals_length(mesh)
        transcribe_finals!(model, phase, mesh)
    end

    return nothing
end

function transcribe_dyn_vars!(
    model::Optimizer,
    i::Integer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM<:AbstractLobattoMesh,MM,BM}
    n_p_dif = get_points_dif_length(mesh)
    n_p_alg = get_points_alg_length(mesh)

    for j in 1:n_p_dif
        for dyn_var in model.dyn_vars[phase]
            vars = model.dyn_var_vars[dyn_var]
            if (j <= n_p_alg) || (dyn_var in model.dif_dyn_vars)
                vars[i][j] = MOI.add_variable(model.inner)
            end
        end
    end

    # add initial constraints
    if i == 1
        transcribe_initials!(model, phase, mesh)
    end

    if i > 1
        for dyn_var in model.dyn_vars[phase]
            if dyn_var in model.dif_dyn_vars
                vars = model.dyn_var_vars[dyn_var]
                MOI.add_constraint(
                    model.inner,
                    1.0 * last(vars[i-1]) - 1.0 * first(vars[i]),
                    MOI.EqualTo(0.0),
                )
            end
        end
    end
    
    # add final constraints
    if i == get_intervals_length(mesh)
        transcribe_finals!(model, phase, mesh)
    end
    return nothing
end

function transcribe_dyn_var_starts!(
    model::Optimizer,
    i::Integer,
    phase::PHS,
    mesh::AbstractIntervalsMesh,
)
    n_p_dif = get_points_dif_length(mesh)
    n_p_alg = get_points_alg_length(mesh)

    points_meshes = get_points_meshes(mesh)

    for j in 1:n_p_dif
        for (dyn_var, start) in model.start_dyn_vars[phase]
            vars = model.dyn_var_vars[dyn_var]

            if dyn_var in model.dif_dyn_vars
                MOI.set(
                    model.inner,
                    MOI.VariablePrimalStart(),
                    vars[i][j],
                    start(points_meshes[i].points_dif[j]),
                )
            elseif j <= n_p_alg
                MOI.set(
                    model.inner,
                    MOI.VariablePrimalStart(),
                    vars[i][j],
                    start(points_meshes[i].points_alg[j]),
                )
            end
        end
    end
    return nothing
end

function transcribe_bounds!(
    model::Optimizer,
    i::Integer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM<:AbstractRadauMesh,MM,BM<:ExactBoundsMesh}

    n_p_dif = get_points_dif_length(mesh)
    n_p_alg = get_points_alg_length(mesh)

    for j in 1:n_p_dif
        for (dyn_var, set) in model.dyn_var_bounds[phase]
            vars = model.dyn_var_vars[dyn_var]
            if dyn_var in model.dif_dyn_vars
                if (j != 1) || (i == 1)
                    MOI.add_constraint(model.inner, vars[i][j], set)
                end
            else
                if j <= n_p_alg
                    MOI.add_constraint(model.inner, vars[i][j], set)
                end
            end
        end
    end
    return nothing
end

function transcribe_bounds!(
    model::Optimizer,
    i::Integer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM<:AbstractLobattoMesh,MM,BM<:ExactBoundsMesh}

    n_p_dif = get_points_dif_length(mesh)
    n_p_alg = get_points_alg_length(mesh)

    for j in 1:n_p_dif
        for (dyn_var, set) in model.dyn_var_bounds[phase]
            vars = model.dyn_var_vars[dyn_var]

            if (j <= n_p_alg) || (dyn_var in model.dif_dyn_vars)
                MOI.add_constraint(model.inner, vars[i][j], set)
            end
        end
    end
    return nothing
end

function transcribe_bounds!(
    model::Optimizer,
    i::Integer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM,BM<:SampledBoundsMesh}

    n_p_dif = get_points_dif_length(mesh)
    n_p_alg = get_points_alg_length(mesh)
    n_p_samp = get_points_samp_length(mesh)

    bounds_mesh = get_bounds_mesh(mesh)

    for j in 1:n_p_samp
        for (dyn_var, set) in model.dyn_var_bounds[phase]
            vars = model.dyn_var_vars[dyn_var]

            if dyn_var in model.dif_dyn_vars
                MOI.add_constraint(
                    model.inner,
                    sum(bounds_mesh.sampled_dif[j,k] * vars[i][k] for k in 1:n_p_dif),
                    set,
                )
            else
                MOI.add_constraint(
                    model.inner,
                    sum(bounds_mesh.sampled_alg[j,k] * vars[i][k] for k in 1:n_p_alg),
                    set,
                )
            end
        end
    end
    return nothing
end

function transcribe_bounds!(
    model::Optimizer,
    i::Integer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM,BM<:BernsteinBoundsMesh}

    n_p_dif = get_points_dif_length(mesh)
    n_p_alg = get_points_alg_length(mesh)

    bounds_mesh = get_bounds_mesh(mesh)

    for j in 1:n_p_dif
        for (dyn_var, set) in model.dyn_var_bounds[phase]
            vars = model.dyn_var_vars[dyn_var]

            if dyn_var in model.dif_dyn_vars
                if j <= n_p_dif
                    MOI.add_constraint(
                        model.inner,
                        sum(bounds_mesh.bernstein_dif[j,k] * vars[i][k] for k in 1:n_p_dif),
                        set,
                    )
                end
            elseif j <= n_p_alg
                MOI.add_constraint(
                    model.inner,
                    sum(bounds_mesh.bernstein_alg[j,k] * vars[i][k] for k in 1:n_p_alg),
                    set,
                )
            end
        end
    end
    return nothing
end

# Collocation, dynamic equations
function transcribe_dif_cons!(
    model::Optimizer,
    i::Integer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:CollocationMesh,BM}

    dif_cons = model.dif_cons[phase]
    n_p_alg = get_points_alg_length(mesh)

    for q in 1:n_p_alg
        for (dif_fun, _, scaling) in values(dif_cons)
            dif_con = transcribe_dyn_fun(
                dif_fun, i, q, model.phase_vars, model.time_vars[phase],
                model.dyn_var_vars, model.dif_dyn_vars, mesh
            )
            scaled_dif_con = apply_scaling(dif_con, scaling)

            MOI.add_constraint(
                model.inner,
                scaled_dif_con,
                MOI.EqualTo(0.0),
            )
            push!(model.res_funcs, scaled_dif_con)
        end
    end
    return nothing
end

# Collocation, algebraic equations
function transcribe_alg_cons!(
    model::Optimizer,
    i::Integer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:CollocationMesh,BM}

    alg_cons = model.alg_cons[phase]
    n_p_alg = get_points_alg_length(mesh)

    for q in 1:n_p_alg
        for (alg_fun, _, scaling) in values(alg_cons)
            alg_con = transcribe_dyn_fun(
                alg_fun, i, q, model.phase_vars, model.time_vars[phase],
                model.dyn_var_vars, model.dif_dyn_vars, mesh
            )
            scaled_alg_con = apply_scaling(alg_con, scaling)

            MOI.add_constraint(
                model.inner,
                scaled_alg_con,
                MOI.EqualTo(0.0),
            )
            push!(model.res_funcs, scaled_alg_con)
        end
    end
    return nothing
end

# Integrated Residual, transcription of differentiation of residuals
# Petrov-Galerkin, differential moment equations
function transcribe_dif_cons!(
    model::Optimizer,
    i::Integer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:GalerkinMesh,BM}

    method_mesh = get_method_mesh(mesh, i)

    for (dif_fun, _, scaling) in values(model.dif_cons[phase])
        for m in 1:size(method_mesh.test_values_dif, 1)
            moment = transcribe_dif_moment(
                model, dif_fun, scaling, i, m, phase, mesh
            )
            MOI.add_constraint(model.inner, moment, MOI.EqualTo(0.0))
            push!(model.res_funcs, moment)
        end
    end

    return nothing
end

function transcribe_dif_cons!(
    ::Optimizer,
    ::Integer,
    ::PHS,
    ::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:Union{AbstractDAIRMesh,QPMMesh},BM}
    return nothing
end

function transcribe_dif_cons!(
    model::Optimizer,
    i::Integer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:Union{SAIRMesh,SAPMMesh},BM}

    grad_res_funcs = transcribe_grad_dyn_least_square(model, i, phase, mesh)

    for f in grad_res_funcs
        MOI.add_constraint(
            model.inner,
            f,
            MOI.EqualTo(0.0)
        )
        push!(model.dif_res_funcs, f)
    end

    if i == get_intervals_length(mesh)
        if model.time_vars[phase] isa VAR
            dyn_res = transcribe_dyn_least_square(model, phase, mesh)
            var = model.time_vars[phase]
            func_t = MOI.Nonlinear.SymbolicAD.derivative(dyn_res, var)
            push!(grad_res_funcs, func_t)
            MOI.add_constraint(model.inner, func_t, MOI.EqualTo(0.0))
        end
    end

    return nothing
end

# Integrated Residual, transcription of residuals
# Petrov-Galerkin, algebraic moment equations
function transcribe_alg_cons!(
    model::Optimizer,
    i::Integer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:GalerkinMesh,BM}

    method_mesh = get_method_mesh(mesh, i)

    for (alg_fun, _, scaling) in values(model.alg_cons[phase])
        for m in 1:size(method_mesh.test_values_alg, 1)
            moment = transcribe_alg_moment(
                model, alg_fun, scaling, i, m, phase, mesh
            )
            MOI.add_constraint(model.inner, moment, MOI.EqualTo(0.0))
            push!(model.res_funcs, moment)
        end
    end

    return nothing
end

function transcribe_alg_cons!(
    ::Optimizer,
    ::Integer,
    ::PHS,
    ::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:Union{DAIRFeasMesh,QPMMesh,SAIRMesh,SAPMMesh},BM}
    return nothing
end

function transcribe_alg_cons!(
    model::Optimizer,
    i::Integer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:DAIROptiMesh,BM}

    method_mesh = get_method_mesh(mesh, i)

    for (dif_fun, _, scaling) in values(model.dif_cons[phase])
        f = transcribe_dif_least_square(model, dif_fun, scaling, i, phase, mesh)
        MOI.add_constraint(
            model.inner,
            f,
            MOI.LessThan(method_mesh.tolerance),
        )
        push!(model.res_funcs, f)
    end

    for (alg_fun, _, scaling) in values(model.alg_cons[phase])
        f = transcribe_alg_least_square(model, alg_fun, scaling, i, phase, mesh)
        MOI.add_constraint(
            model.inner,
            f,
            MOI.LessThan(method_mesh.tolerance),
        )
        push!(model.res_funcs, f)
    end

    return nothing
end

function transcribe_path_cons!(
    model::Optimizer,
    i::Integer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:CollocationMesh,BM}

    path_cons = model.path_cons[phase]
    n_p_alg = get_points_alg_length(mesh)

    for q in 1:n_p_alg
        for (path_fun, set) in values(path_cons)
            path_con = transcribe_dyn_fun(
                path_fun, i, q, model.phase_vars, model.time_vars[phase],
                model.dyn_var_vars, model.dif_dyn_vars, mesh
            )
            MOI.add_constraint(model.inner, path_con, set)
        end
    end
    return nothing
end

function transcribe_path_cons!(
    model::Optimizer,
    i::Integer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:AbstractIntResMesh,BM}

    path_cons = model.path_cons[phase]
    n_p_quad = get_points_quad_length(mesh)

    for q in 1:n_p_quad
        for (path_fun, set) in values(path_cons)
            path_con = transcribe_dyn_fun(
                path_fun, i, q, model.phase_vars, model.time_vars[phase],
                model.dyn_var_vars, model.dif_dyn_vars, mesh
            )
            MOI.add_constraint(model.inner, path_con, set)
        end
    end
    return nothing
end

# Boundary constraints
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
