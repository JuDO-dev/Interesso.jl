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
    ::Optimizer,
    meshes::MESHES,
)
    return transcribe_phase_bou(phase_bou, meshes[phase_bou.dyn_fun])
end

transcribe_phase_bou(::DOI.Initial{PHS}, mesh::FixedIntervalsMesh) = mesh.points_meshes[1].t_a
transcribe_phase_bou(::DOI.Final{PHS}, mesh::FixedIntervalsMesh) = mesh.points_meshes[end].t_b
transcribe_phase_bou(::DOI.Initial{PHS}, mesh::FlexibleIntervalsMesh) = mesh.fixed.points_meshes[1].t_a
transcribe_phase_bou(::DOI.Final{PHS}, mesh::FlexibleIntervalsMesh) = mesh.fixed.points_meshes[end].t_b

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
) where {PM<:LGRPointsMesh,MM,BM}

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
) where {PM<:LGRPointsMesh,MM,BM}

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
            transcribe_integral(dyn_fun, model, meshes[DOI.phase_index(dyn_fun)])
            for dyn_fun in integrals.dyn_funs
        ],
    )
end

function transcribe_integral(integrand::NDF, model::Optimizer, mesh::FixedIntervalsMesh)

    n_h = get_intervals_length(mesh)
    n_p_alg = get_points_alg_length(mesh)

    return MOI.ScalarNonlinearFunction(
        :+,
        [MOI.ScalarNonlinearFunction(
            :+,
            [MOI.ScalarNonlinearFunction(
                :*,
                [
                    mesh.points_meshes[i].quad_weights[q],
                    transcribe_dyn_fun(
                        integrand,
                        i,
                        q,
                        model.phase_vars,
                        model.dyn_var_vars,
                        model.dif_dyn_vars,
                        mesh,
                    ),
                ],
            ) for q in 1:n_p_alg],
        ) for i in 1:n_h],
    )
end

function transcribe_integral(integrand::NDF, model::Optimizer, mesh::FlexibleIntervalsMesh)

    n_h = get_intervals_length(mesh)
    n_p_alg = get_points_alg_length(mesh)
    t_0 = mesh.fixed.points_meshes[1].t_a
    t_f = mesh.fixed.points_meshes[end].t_b

    flex_vars = model.phase_vars[DOI.phase_index(integrand)]

    Δt_1 = 1.0 * flex_vars[1] - t_0
    Δt_n_h = t_f - 1.0 * flex_vars[end]
    Δt_inner = [1.0 * flex_vars[i] - 1.0 * flex_vars[i-1] for i in 2:(n_h-1)]
    Δt = vcat(Δt_1, Δt_inner, Δt_n_h)

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
                            integrand, i, q, model.phase_vars, model.dyn_var_vars,
                            model.dif_dyn_vars, mesh,
                        ),
                    ],
                ) for q in 1:n_p_alg],  
            )],
        ) for i in 1:n_h],
    )
end

# Bolza
function transcribe_bou_fun(bolza::OBJ, model::Optimizer, meshes::MESHES)
    return MOI.ScalarNonlinearFunction(
        :+,
        [
            transcribe_bou_fun(bolza.bou_fun, model, meshes),
            transcribe_bou_fun(bolza.integral, model, meshes),
        ],
    )
end