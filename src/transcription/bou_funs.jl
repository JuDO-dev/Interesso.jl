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

function transcribe_integral(
    integrand::NDF,
    model::Optimizer,
    mesh::FixedIntervalsMesh{PM,MM,BM}
) where {PM,MM<:CollocationMesh,BM}

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
                        integrand, i, q, model.phase_vars, model.dyn_var_vars,
                        model.dif_dyn_vars, mesh,
                    ),
                ],
            ) for q in 1:n_p_alg],
        ) for i in 1:n_h],
    )
end

function transcribe_integral(
    integrand::NDF,
    model::Optimizer,
    mesh::FlexibleIntervalsMesh{PM,MM,BM}
) where {PM,MM<:CollocationMesh,BM}

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

function transcribe_integral(
    integrand::NDF,
    model::Optimizer,
    mesh::FixedIntervalsMesh{PM,MM,BM}
) where {PM,MM<:PenaltyIRMesh,BM}

    n_h = get_intervals_length(mesh)
    n_p_quad = get_points_quad_length(mesh)

    return MOI.ScalarNonlinearFunction(
        :+,
        [MOI.ScalarNonlinearFunction(
            :+,
            [MOI.ScalarNonlinearFunction(
                :*,
                [
                    mesh.points_meshes[i].quad_weights[q],
                    transcribe_dyn_fun(
                        integrand, i, q, model.phase_vars, model.dyn_var_vars,
                        model.dif_dyn_vars, mesh,
                    ),
                ],
            ) for q in 1:n_p_quad],
        ) for i in 1:n_h],
    )
end

function transcribe_integral(
    integrand::NDF,
    model::Optimizer,
    mesh::FlexibleIntervalsMesh{PM,MM,BM}
) where {PM,MM<:PenaltyIRMesh,BM}

    n_h = get_intervals_length(mesh)
    n_p_quad = get_points_quad_length(mesh)
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
                        mesh.method_mesh.quad_points_mesh.quad_weights[q],
                        transcribe_dyn_fun(
                            integrand, i, q, model.phase_vars, model.dyn_var_vars,
                            model.dif_dyn_vars, mesh,
                        ),
                    ],
                ) for q in 1:n_p_quad],  
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

function transcribe_dif_least_square(
    model::Optimizer,
    i::Integer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:PenaltyIRMesh,BM}

    dif_cons = model.dif_cons[phase]
    n_p_quad = get_points_quad_length(mesh)

    return MOI.ScalarNonlinearFunction(
        :+,
        [   
            MOI.ScalarNonlinearFunction(
                :^,
                [
                    transcribe_dyn_fun(dif_fun, i, q, model.phase_vars,
                        model.dyn_var_vars, model.dif_dyn_vars, mesh
                    ),
                    2.0
                ] 
            ) for (dif_fun, _) in values(dif_cons) for q in 1:n_p_quad
        ]
    ) 
end

function transcribe_dif_least_square(
    model::Optimizer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:PenaltyIRMesh,BM}

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
    i::Integer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:PenaltyIRMesh,BM}

    alg_cons = model.alg_cons[phase]
    n_p_quad = get_points_quad_length(mesh)

    return MOI.ScalarNonlinearFunction(
        :+,
        [
            MOI.ScalarNonlinearFunction(
                :^,
                [
                    transcribe_dyn_fun(alg_fun, i, q, model.phase_vars,
                        model.dyn_var_vars, model.dif_dyn_vars, mesh
                    ), 
                    2.0
                ]
            ) for (alg_fun, _) in values(alg_cons) for q in 1:n_p_quad
        ]
    )
end

function transcribe_alg_least_square(
    model::Optimizer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:PenaltyIRMesh,BM}

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
) where {PM,MM<:PenaltyIRMesh,BM}

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
) where {PM,MM<:PenaltyIRMesh,BM}

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
            transcribe_dyn_least_square(model, phase, mesh) for (phase, mesh) in meshes if mesh.method_mesh isa PenaltyIRMesh
        ]
    )
end

function transcribe_grad_dyn_least_square(
    model::Optimizer,
    i::Integer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:PenaltyIRMesh,BM}

"""
    for an interval i, it will include unique(dyn_var_vars[all dyn_vars][i]) optimizers
    
    the residual sum of an interval is, providing i, (phase, mesh) from meshes

    MOI.ScalarNonlinearFunction(
        :+,
        [
            transcribe_dyn_least_square(model, i, phase, mesh)
        ]
    )
"""

    grad_res_funcs = Vector{MOI.AbstractFunction}()
    interval_vars = _get_interval_dyn_vars(model, i, phase)

    for dyn_var in interval_vars
        
        func = MOI.Nonlinear.SymbolicAD.derivative(
                    transcribe_dyn_least_square(model, i, phase, mesh),
                    dyn_var
                )

        push!(grad_res_funcs, func)
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
        append!(interval_vars, model.dyn_var_vars[dyn_var][i])
    end

    return unique!(interval_vars)
end