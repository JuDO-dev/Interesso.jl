# Numbers
function transcribe_bou_fun(
    number::Real,
    ::Optimizer,
    ::MESHES,
)
    return number
end

# MOI Functions
function transcribe_bou_fun(
    fun::MOI.AbstractScalarFunction,
    ::Optimizer,
    ::MESHES,
)
    return fun
end

# Phase Boundaries
function transcribe_bou_fun(
    phase_bou::Union{DOI.Initial{PHS},DOI.Final{PHS}},
    model::Optimizer,
    meshes::MESHES,
)
    return transcribe_phase_bou(
        phase_bou,
        model.phase_vars,
        model.time_vars,
        meshes[phase_bou.dyn_fun],
    )
end

transcribe_phase_bou(
    ::DOI.Initial{PHS},
    ::PHS_VARS,
    ::TIME_VARS,
    mesh::FixedIntervalsMesh,
) = mesh.points_meshes[1].t_a

function transcribe_phase_bou(
    phase_final::DOI.Final{PHS},
    ::PHS_VARS,
    time_vars::TIME_VARS,
    mesh::FixedIntervalsMesh,
)
    phase = phase_final.dyn_fun
    if haskey(time_vars, phase)
        return 1.0 * time_vars[phase]
    else
        return mesh.points_meshes[end].t_b
    end
end

transcribe_phase_bou(
    ::DOI.Initial{PHS},
    ::PHS_VARS,
    ::TIME_VARS,
    mesh::FlexibleIntervalsMesh,
) = mesh.fixed.points_meshes[1].t_a

transcribe_phase_bou(
    ::DOI.Final{PHS},
    ::PHS_VARS,
    ::TIME_VARS,
    mesh::FlexibleIntervalsMesh,
) = mesh.fixed.points_meshes[end].t_b

# Dynamic Variable Boundaries
function transcribe_bou_fun(
    dyn_var_initial::DOI.Initial{DYN_VAR},
    model::Optimizer,
    meshes::MESHES,
)
    return transcribe_dyn_var_initial(
        dyn_var_initial.dyn_fun,
        model.dyn_var_vars,
        meshes[DOI.phase_index(dyn_var_initial.dyn_fun)]
    )
end

function transcribe_dyn_var_initial(
    dyn_var::DYN_VAR,
    dyn_var_vars::DYN_VAR_VARS,
    ::AbstractIntervalsMesh{PM,MM,BM},
) where {PM<:AbstractPointsMesh,MM,BM}

    return 1.0 * dyn_var_vars[dyn_var][1][1]
end

function transcribe_bou_fun(
    dyn_var_final::DOI.Final{DYN_VAR},
    model::Optimizer,
    meshes::MESHES,
)
    return transcribe_dyn_var_final(
        dyn_var_final.dyn_fun,
        model.dyn_var_vars,
        meshes[DOI.phase_index(dyn_var_final.dyn_fun)]
    )
end

function transcribe_dyn_var_final(
    dyn_var::DYN_VAR,
    dyn_var_vars::DYN_VAR_VARS,
    ::AbstractIntervalsMesh{PM,MM,BM},
) where {PM<:AbstractPointsMesh,MM,BM}

    return 1.0 * dyn_var_vars[dyn_var][end][end]
end

# NonlinearBoundaryFunction
function transcribe_bou_fun(
    nl_bou_fun::NBF,
    model::Optimizer,
    meshes::MESHES,
)
    return MOI.ScalarNonlinearFunction(
        nl_bou_fun.head,
        [transcribe_bou_fun(arg, model, meshes) for arg in nl_bou_fun.args],
    )
end

# Multi-Phase Integral
function transcribe_bou_fun(
    integrals::DOI.MultiPhaseIntegral{NDF},
    model::Optimizer,
    meshes::MESHES,
)
    return MOI.ScalarNonlinearFunction(
        :+,
        [
            transcribe_integral(dyn_fun, model, meshes[DOI.phase_index(dyn_fun)], model.time_vars[DOI.phase_index(dyn_fun)])
            for dyn_fun in integrals.dyn_funs
        ],
    )
end

function transcribe_integral(
    integrand::NDF,
    model::Optimizer,
    mesh::FixedIntervalsMesh{PM,MM,BM},
    time_var::Union{Float64, VAR}
) where {PM,MM<:CollocationMesh,BM}

    n_h = get_intervals_length(mesh)
    n_p_alg = get_points_alg_length(mesh)

    return MOI.ScalarNonlinearFunction(
        :+,
        [MOI.ScalarNonlinearFunction(
            :*,
            Any[time_var, MOI.ScalarNonlinearFunction(
                :+,
                [MOI.ScalarNonlinearFunction(
                    :*,
                    [
                        mesh.points_meshes[i].quad_weights[q],
                        transcribe_dyn_fun(
                            integrand, i, q, model.phase_vars, time_var,
                            model.dyn_var_vars, model.dif_dyn_vars, mesh,
                        ),
                    ],
                ) for q in 1:n_p_alg],
            )],
        ) for i in 1:n_h],
    )
end

function transcribe_integral(
    integrand::NDF,
    model::Optimizer,
    mesh::FlexibleIntervalsMesh{PM,MM,BM},
    time_var::Union{Float64, VAR}
) where {PM,MM<:CollocationMesh,BM}

    n_h = get_intervals_length(mesh)
    n_p_alg = get_points_alg_length(mesh)
    t_0 = mesh.fixed.points_meshes[1].t_a
    t_f = mesh.fixed.points_meshes[end].t_b

    flex_vars = model.phase_vars[DOI.phase_index(integrand)]

    Δt_1 = 1.0 * flex_vars[1] - t_0
    Δt_n_h = t_f - 1.0 * flex_vars[end]
    Δt_inner = [1.0 * flex_vars[i] - 1.0 * flex_vars[i-1] for i in 2:(n_h-1)]
    Δt = vcat(Δt_1, Δt_inner, Δt_n_h) .* time_var

    return MOI.ScalarNonlinearFunction(
        :+,
        [MOI.ScalarNonlinearFunction(
            :*,
            Any[0.5, Δt[i], MOI.ScalarNonlinearFunction(
                :+,
                [MOI.ScalarNonlinearFunction(
                    :*,
                    [
                        mesh.points_mesh.quad_weights[q],
                        transcribe_dyn_fun(
                            integrand, i, q, model.phase_vars, time_var,
                            model.dyn_var_vars, model.dif_dyn_vars, mesh,
                        ),
                    ],
                ) for q in 1:n_p_alg],  
            )],
        ) for i in 1:n_h],
    )
end

function transcribe_integral(
    integrand::NDF,
    model::Optimizer,
    mesh::FixedIntervalsMesh{PM,MM,BM},
    time_var::Union{Float64, VAR}
) where {PM,MM<:AbstractIntResMesh,BM}

    n_h = get_intervals_length(mesh)
    n_p_quad = get_points_quad_length(mesh)

    return MOI.ScalarNonlinearFunction(
        :+,
        [MOI.ScalarNonlinearFunction(
            :*,
            Any[time_var, MOI.ScalarNonlinearFunction(
                :+,
                [MOI.ScalarNonlinearFunction(
                    :*,
                    [
                        mesh.method_meshes[i].quad_points_mesh.quad_weights[q],
                        transcribe_dyn_fun(
                            integrand, i, q, model.phase_vars, time_var,
                            model.dyn_var_vars, model.dif_dyn_vars, mesh,
                        ),
                    ],
                ) for q in 1:n_p_quad],
            )],
        ) for i in 1:n_h],
    )
end

function transcribe_integral(
    integrand::NDF,
    model::Optimizer,
    mesh::FlexibleIntervalsMesh{PM,MM,BM},
    time_var::Union{Float64, VAR}
) where {PM,MM<:AbstractIntResMesh,BM}

    n_h = get_intervals_length(mesh)
    n_p_quad = get_points_quad_length(mesh)
    t_0 = mesh.fixed.points_meshes[1].t_a
    t_f = mesh.fixed.points_meshes[end].t_b

    flex_vars = model.phase_vars[DOI.phase_index(integrand)]

    Δt_1 = 1.0 * flex_vars[1] - t_0
    Δt_n_h = t_f - 1.0 * flex_vars[end]
    Δt_inner = [1.0 * flex_vars[i] - 1.0 * flex_vars[i-1] for i in 2:(n_h-1)]
    Δt = vcat(Δt_1, Δt_inner, Δt_n_h) .* time_var

    return MOI.ScalarNonlinearFunction(
        :+,
        [MOI.ScalarNonlinearFunction(
            :*,
            Any[0.5, Δt[i], MOI.ScalarNonlinearFunction(
                :+,
                [MOI.ScalarNonlinearFunction(
                    :*,
                    [
                        mesh.method_mesh.quad_points_mesh.quad_weights[q],
                        transcribe_dyn_fun(
                            integrand, i, q, model.phase_vars, time_var,
                            model.dyn_var_vars, model.dif_dyn_vars, mesh,
                        ),
                    ],
                ) for q in 1:n_p_quad],  
            )],
        ) for i in 1:n_h],
    )
end


# Petrov-Galerkin residual moments
function transcribe_dif_moment(
    model::Optimizer,
    dif_fun::DIF_FUN,
    scaling::Float64,
    i::Integer,
    m::Integer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:GalerkinMesh,BM}

    n_p_quad = get_points_quad_length(mesh)
    method_mesh = get_method_mesh(mesh, i)

    return MOI.ScalarNonlinearFunction(
        :+,
        [
            MOI.ScalarNonlinearFunction(:*, Any[
                method_mesh.quad_points_mesh.quad_weights[q] * method_mesh.test_values_dif[m, q],
                transcribe_dyn_fun(
                    dif_fun, i, q, model.phase_vars, model.time_vars[phase],
                    model.dyn_var_vars, model.dif_dyn_vars, mesh, scaling
                ),
            ]) for q in 1:n_p_quad
        ],
    )
end

function transcribe_alg_moment(
    model::Optimizer,
    alg_fun::NDF,
    scaling::Float64,
    i::Integer,
    m::Integer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:GalerkinMesh,BM}

    n_p_quad = get_points_quad_length(mesh)
    method_mesh = get_method_mesh(mesh, i)

    return MOI.ScalarNonlinearFunction(
        :+,
        [
            MOI.ScalarNonlinearFunction(:*, Any[
                method_mesh.quad_points_mesh.quad_weights[q] * method_mesh.test_values_alg[m, q],
                transcribe_dyn_fun(
                    alg_fun, i, q, model.phase_vars, model.time_vars[phase],
                    model.dyn_var_vars, model.dif_dyn_vars, mesh, scaling
                )
            ]) for q in 1:n_p_quad
        ],
    )
end

# Bolza
function transcribe_bou_fun(bolza::BOLZA, model::Optimizer, meshes::MESHES)
    return MOI.ScalarNonlinearFunction(
        :+,
        [
            transcribe_bou_fun(bolza.bou_fun, model, meshes),
            transcribe_bou_fun(bolza.integral, model, meshes),
        ],
    )
end

# least-square dynamics
function lift!(
    model::Optimizer,
    f::MOI.AbstractScalarFunction,
)
    s = MOI.add_variable(model.inner)
    MOI.set(model.inner, MOI.VariablePrimalStart(), s, 0.0)
    push!(model.lift_vars, s)
    MOI.add_constraint(
        model.inner,
        MOI.ScalarNonlinearFunction(:-, Any[s, f]),
        MOI.EqualTo(0.0),
    )
    return s
end

function transcribe_dif_least_square(
    model::Optimizer,
    dif_fun::DIF_FUN,
    scaling::Float64,
    i::Integer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:AbstractIntResMesh,BM}

    n_p_quad = get_points_quad_length(mesh)
    method_mesh = get_method_mesh(mesh, i)
    
    return MOI.ScalarNonlinearFunction(
        :+,
        [
            begin
                f = transcribe_dyn_fun(
                    dif_fun, i, q, model.phase_vars, model.time_vars[phase],
                    model.dyn_var_vars, model.dif_dyn_vars, mesh, scaling
                )
                s = lift!(model, f)
                MOI.ScalarNonlinearFunction(:*, Any[
                    method_mesh.quad_points_mesh.quad_weights[q],
                    MOI.ScalarNonlinearFunction(:^, [s, 2.0])
                ])
            end for q in 1:n_p_quad
        ]
    ) 
end

function transcribe_dif_least_square(
    model::Optimizer,
    i::Integer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:AbstractIntResMesh,BM}

    dif_cons = model.dif_cons[phase]

    return MOI.ScalarNonlinearFunction(
        :+,
        [
            transcribe_dif_least_square(
                model, dif_fun, scaling, i, phase, mesh
            ) for (dif_fun, _, scaling) in values(dif_cons)
        ]
    ) 
end

function transcribe_dif_least_square(
    model::Optimizer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:AbstractIntResMesh,BM}

    n_h = get_intervals_length(mesh)

    return MOI.ScalarNonlinearFunction(
        :+,
        [   
            transcribe_dif_least_square(model, i, phase, mesh) for i in 1:n_h
        ]
    ) 
end

function transcribe_alg_least_square(
    model::Optimizer,
    alg_fun::NDF,
    scaling::Float64,
    i::Integer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:AbstractIntResMesh,BM}

    n_p_quad = get_points_quad_length(mesh)
    method_mesh = get_method_mesh(mesh, i)
    
    return MOI.ScalarNonlinearFunction(
        :+,
        [
            begin
                f = transcribe_dyn_fun(
                    alg_fun, i, q, model.phase_vars, model.time_vars[phase],
                    model.dyn_var_vars, model.dif_dyn_vars, mesh, scaling
                )
                s = lift!(model, f)
                MOI.ScalarNonlinearFunction(:*, Any[
                    method_mesh.quad_points_mesh.quad_weights[q],
                    MOI.ScalarNonlinearFunction(:^, [s, 2.0])
                ])
            end for q in 1:n_p_quad
        ]
    ) 
end

function transcribe_alg_least_square(
    model::Optimizer,
    i::Integer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:AbstractIntResMesh,BM}

    alg_cons = model.alg_cons[phase]

    return MOI.ScalarNonlinearFunction(
        :+,
        [
            transcribe_alg_least_square(
                model, alg_fun, scaling, i, phase, mesh
            ) for (alg_fun, _, scaling) in values(alg_cons)
        ]
    )
end

function transcribe_alg_least_square(
    model::Optimizer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:AbstractIntResMesh,BM}

    n_h = get_intervals_length(mesh)

    return MOI.ScalarNonlinearFunction(
        :+,
        [   
            transcribe_alg_least_square(model, i, phase, mesh) for i in 1:n_h
        ]
    ) 
end

function transcribe_dyn_least_square(
    model::Optimizer,
    i::Integer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:AbstractIntResMesh,BM}

    return MOI.ScalarNonlinearFunction(
        :+,
        [
            transcribe_dif_least_square(model, i, phase, mesh),
            transcribe_alg_least_square(model, i, phase, mesh),
        ]
    )
end

function transcribe_dyn_least_square(
    model::Optimizer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:AbstractIntResMesh,BM}

    return MOI.ScalarNonlinearFunction(
        :+,
        [
            transcribe_dif_least_square(model, phase, mesh),
            transcribe_alg_least_square(model, phase, mesh),
        ]
    )
end

function transcribe_dyn_least_square(
    model::Optimizer,
    meshes::MESHES,
)
    return MOI.ScalarNonlinearFunction(
        :+,
        [
            transcribe_dyn_least_square(model, phase, mesh) for (phase, mesh) in meshes if mesh.method_mesh isa Union{QPMMesh,SAPMMesh}
        ]
    )
end

function transcribe_dyn_fun_derivative(
    f::MOI.AbstractScalarFunction,
    var::VAR,
)
    return MOI.Nonlinear.SymbolicAD.derivative(f, var)
end

function transcribe_grad_dif_dyn(
    model::Optimizer,
    i::Integer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:AbstractIntResMesh,BM}

    grad_res_funcs = Vector{MOI.AbstractFunction}()
    vars = _get_interval_dyn_vars(model, i, phase)
    grad_terms = [MOI.AbstractFunction[] for _ in vars]
    path_terms = _path_multiplier_terms!(model, i, phase, mesh, vars)

    n_p_quad = get_points_quad_length(mesh)
    method_mesh = get_method_mesh(mesh, i)

    for cons in (values(model.dif_cons[phase]), values(model.alg_cons[phase]))
        for (dyn_fun, _, scaling) in cons
            for q in 1:n_p_quad
                f = transcribe_dyn_fun(
                    dyn_fun, i, q, model.phase_vars, model.time_vars[phase],
                    model.dyn_var_vars, model.dif_dyn_vars, mesh, scaling
                )
                s = lift!(model, f)

                for (j, dyn_var) in enumerate(vars)
                    df = transcribe_dyn_fun_derivative(f, dyn_var)
                    push!(
                        grad_terms[j],
                        MOI.ScalarNonlinearFunction(:*, Any[
                            method_mesh.quad_points_mesh.quad_weights[q],
                            s, df,
                        ])
                    )
                end
            end
        end
    end

    for (j, dyn_var) in enumerate(vars)
        terms = grad_terms[j]
        if haskey(path_terms, dyn_var)
            push!(terms, path_terms[dyn_var])
        end
        push!(grad_res_funcs, MOI.ScalarNonlinearFunction(:+, terms))
    end

    return grad_res_funcs
end

function transcribe_grad_dyn_least_square(
    model::Optimizer,
    i::Integer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:AbstractIntResMesh,BM}

    grad_res_funcs = Vector{MOI.AbstractFunction}()
    funcs_dif = transcribe_grad_dif_dyn(model, i, phase, mesh)
    append!(grad_res_funcs, funcs_dif)

    return grad_res_funcs
end

function transcribe_grad_dyn_least_square(
    model::Optimizer,
    var::VAR,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:AbstractIntResMesh,BM}

    n_h = get_intervals_length(mesh)
    grad_terms = MOI.AbstractFunction[]

    for i in 1:n_h
        n_p_quad = get_points_quad_length(mesh)
        method_mesh = get_method_mesh(mesh, i)

        for cons in (values(model.dif_cons[phase]), values(model.alg_cons[phase]))
            for (dyn_fun, _, scaling) in cons
                for q in 1:n_p_quad
                    f = transcribe_dyn_fun(
                        dyn_fun, i, q, model.phase_vars, model.time_vars[phase],
                        model.dyn_var_vars, model.dif_dyn_vars, mesh, scaling
                    )
                    s = lift!(model, f)
                    df = transcribe_dyn_fun_derivative(f, var)

                    push!(
                        grad_terms,
                        MOI.ScalarNonlinearFunction(:*, Any[
                            method_mesh.quad_points_mesh.quad_weights[q],
                            s, df,
                        ])
                    )
                end
            end
        end
    end

    return MOI.ScalarNonlinearFunction(:+, grad_terms)
end

function transcribe_grad_dyn_least_square(
    model::Optimizer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:AbstractIntResMesh,BM}

    n_h = get_intervals_length(mesh)
    grad_res_funcs = Vector{MOI.AbstractFunction}()

    for i = 1:n_h
        append!(
            grad_res_funcs,
            transcribe_grad_dyn_least_square(model, i, phase, mesh)
        )
    end

    if model.time_vars[phase] isa VAR
        push!(
            grad_res_funcs,
            transcribe_grad_dyn_least_square(
                model, model.time_vars[phase], phase, mesh
            )
        )
    end

    return grad_res_funcs
end

function _get_interval_dyn_vars(
    model::Optimizer,
    i::Integer,
    phase::PHS
)

    interval_vars = VAR[]
    for dyn_var in model.dyn_vars[phase]
        if dyn_var in model.dif_dyn_vars
            append!(interval_vars, model.dyn_var_vars[dyn_var][i])
            filter!(var -> (var != model.dyn_var_vars[dyn_var][i][1]), interval_vars)  # remove the first one for continuity
        elseif !(dyn_var in model.ctrl_vars)
            append!(interval_vars, model.dyn_var_vars[dyn_var][i])  # algebraic variable: differentiated, no shared continuity node
        end
    end

    unique!(interval_vars)

    return unique!(interval_vars)
end

function _add_nonnegative_variable!(model::Optimizer)
    μ = MOI.add_variable(model.inner)
    MOI.set(model.inner, MOI.VariablePrimalStart(), μ, 0.0)
    MOI.add_constraint(model.inner, μ, MOI.GreaterThan(0.0))
    return MOI.ScalarNonlinearFunction(:+, [μ])
end

# function _add_nonnegative_variable!(model::Optimizer)
#     μ = MOI.add_variable(model.inner)
#     return MOI.ScalarNonlinearFunction(:^, [μ, 2.0])
# end

function _add_complementarity!(
    model::Optimizer,
    μ::MOI.ScalarNonlinearFunction,
    residual::MOI.ScalarNonlinearFunction,
)
    comp = MOI.ScalarNonlinearFunction(:*, Any[μ, residual])
    MOI.add_constraint(model.inner, comp, MOI.EqualTo(0.0))
    return comp
end

function _path_upper_residual(path_con::MOI.ScalarNonlinearFunction, upper::Float64)
    return MOI.ScalarNonlinearFunction(:-, Any[path_con, upper])
end

function _path_lower_residual(path_con::MOI.ScalarNonlinearFunction, lower::Float64)
    return MOI.ScalarNonlinearFunction(:-, Any[lower, path_con])
end

function _accumulate_path_derivative_term!(
    terms::Dict{VAR,MOI.AbstractFunction},
    var::VAR,
    op::Symbol,
    term::MOI.AbstractFunction,
)
    if haskey(terms, var)
        terms[var] = MOI.ScalarNonlinearFunction(op, Any[terms[var], term])
    elseif op === :+
        terms[var] = term
    elseif op === :-
        terms[var] = MOI.ScalarNonlinearFunction(:-, Any[term])
    else
        error("Unsupported path derivative operator: $(op)")
    end

    return nothing
end

function _add_path_derivative_term!(
    terms::Dict{VAR,MOI.AbstractFunction},
    var::VAR,
    μ::MOI.ScalarNonlinearFunction,
    op::Symbol,
    path_con::MOI.ScalarNonlinearFunction,
)
    dpath = MOI.Nonlinear.SymbolicAD.derivative(path_con, var)
    term = MOI.ScalarNonlinearFunction(:*, Any[μ, dpath])

    return _accumulate_path_derivative_term!(terms, var, op, term)
end

function _path_multiplier_terms!(
    model::Optimizer,
    i::Integer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
    vars::Vector{VAR},
) where {PM,MM<:AbstractIntResMesh,BM}

    terms = Dict{VAR,MOI.AbstractFunction}()
    n_p_quad = get_points_quad_length(mesh)

    for q in 1:n_p_quad
        for (path_fun, set) in values(model.path_cons[phase])
            path_con = transcribe_dyn_fun(
                path_fun, i, q, model.phase_vars, model.time_vars[phase],
                model.dyn_var_vars, model.dif_dyn_vars, mesh,
            )

            if set isa LE64
                μ = _add_nonnegative_variable!(model)
                residual = _path_upper_residual(path_con, set.upper)
                _add_complementarity!(model, μ, residual)
                for var in vars
                    _add_path_derivative_term!(terms, var, μ, :+, path_con)
                end

            elseif set isa GE64
                μ = _add_nonnegative_variable!(model)
                residual = _path_lower_residual(path_con, set.lower)
                _add_complementarity!(model, μ, residual)
                for var in vars
                    _add_path_derivative_term!(terms, var, μ, :-, path_con)
                end

            elseif set isa IV64
                μ = _add_nonnegative_variable!(model)
                residual = _path_lower_residual(path_con, set.lower)
                _add_complementarity!(model, μ, residual)
                for var in vars
                    _add_path_derivative_term!(terms, var, μ, :-, path_con)
                end

                μ = _add_nonnegative_variable!(model)
                residual = _path_upper_residual(path_con, set.upper)
                _add_complementarity!(model, μ, residual)
                for var in vars
                    _add_path_derivative_term!(terms, var, μ, :+, path_con)
                end
            end
        end
    end

    return terms
end
