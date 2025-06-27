# Fallback
function transcribe_dyn_fun(
    ::Any,
    ::Integer,
    ::Integer,
    ::PHS_VARS,
    ::DYN_VAR_VARS,
    ::AbstractSet{DYN_VAR},
    ::AbstractIntervalsMesh,
)
    return throw(ArgumentError("Unsupported"))
end

# Number
function transcribe_dyn_fun(
    number::Real,
    ::Integer,
    ::Integer,
    ::PHS_VARS,
    ::DYN_VAR_VARS,
    ::AbstractSet{DYN_VAR},
    ::AbstractIntervalsMesh,
)
    return number
end

# MOI Functions
function transcribe_dyn_fun(
    fun::MOI.AbstractScalarFunction,
    ::Integer,
    ::Integer,
    ::PHS_VARS,
    ::DYN_VAR_VARS,
    ::AbstractSet{DYN_VAR},
    ::AbstractIntervalsMesh,
)
    return fun
end

# Phase
function transcribe_dyn_fun(
    ::PHS,
    i::Integer,
    q::Integer,
    ::PHS_VARS,
    ::DYN_VAR_VARS,
    ::AbstractSet{DYN_VAR},
    mesh::FixedIntervalsMesh,
)
    return mesh.points_meshes[i].points_alg[q]
end

function transcribe_dyn_fun(
    phase::PHS,
    i::Integer,
    q::Integer,
    phase_vars::PHS_VARS,
    ::DYN_VAR_VARS,
    ::AbstractSet{DYN_VAR},
    mesh::FlexibleIntervalsMesh,
)
    n_h = get_intervals_length(mesh)
    flex_vars = phase_vars[phase]

    if i == 1
        t_a = mesh.fixed.points_meshes[1].t_a
        t_b = flex_vars[1]

    elseif i == n_h
        t_a = flex_vars[end]
        t_b = mesh.fixed.points_meshes[end].t_b
    else
        t_a = flex_vars[i]
        t_b = flex_vars[i + 1]
    end

    return 0.5 * (t_a + t_b) + 0.5 * (t_b - t_a) * mesh.points_mesh.points_alg[q]
end

# Dynamic Variable
function transcribe_dyn_fun(
    dyn_var::DYN_VAR,
    i::Integer,
    q::Integer,
    ::PHS_VARS,
    dyn_var_vars::DYN_VAR_VARS,
    ::AbstractSet{DYN_VAR},
    ::AbstractIntervalsMesh,
)
    return dyn_var_vars[dyn_var][i][q]
end

# Nonlinear Dynamic Function
function transcribe_dyn_fun(
    nl_dyn_fun::DOI.NonlinearDynamicFunction,
    i::Integer,
    q::Integer,
    phase_vars::PHS_VARS,
    dyn_var_vars::DYN_VAR_VARS,
    dif_dyn_vars::AbstractSet{DYN_VAR},
    mesh::AbstractIntervalsMesh,
)
    return MOI.ScalarNonlinearFunction(
        nl_dyn_fun.head,
        [transcribe_dyn_fun(
            arg, i, q, phase_vars, dyn_var_vars, dif_dyn_vars, mesh,
        ) for arg in nl_dyn_fun.args],
    )
end

# Explicit Differential Function
function transcribe_dyn_fun(
    dif_fun::DIF_FUN,
    i::Integer,
    q::Integer,
    phase_vars::PHS_VARS,
    dyn_var_vars::DYN_VAR_VARS,
    dif_dyn_vars::AbstractSet{DYN_VAR},
    mesh::FixedIntervalsMesh,
)
    vars = dyn_var_vars[dif_fun.dyn_var]
    n_p_dif = get_points_dif_length(mesh)

    return MOI.ScalarNonlinearFunction(:-, Any[
        sum(mesh.points_meshes[i].differentiation[q,k] * vars[i][k] for k in 1:n_p_dif),
        transcribe_dyn_fun(dif_fun.dyn_fun, i, q, phase_vars, dyn_var_vars, dif_dyn_vars, 
            mesh
        ),
    ])
end

function transcribe_dyn_fun(
    dif_fun::DIF_FUN,
    i::Integer,
    q::Integer,
    phase_vars::PHS_VARS,
    dyn_var_vars::DYN_VAR_VARS,
    dif_dyn_vars::AbstractSet{DYN_VAR},
    mesh::FlexibleIntervalsMesh,
)
    vars = dyn_var_vars[dif_fun.dyn_var]
    n_p_dif = get_points_dif_length(mesh)

    flex_vars = phase_vars[DOI.phase_index(dif_fun)]
    n_h = get_intervals_length(mesh)
    t_0 = mesh.fixed.points_meshes[1].t_a
    t_f = mesh.fixed.points_meshes[end].t_b

    if i == 1
        Δt = 1.0 * first(flex_vars) - t_0
    elseif i == n_h
        Δt = t_f - 1.0 * last(flex_vars)
    else
        Δt = 1.0 * flex_vars[i] - 1.0 * flex_vars[i-1]
    end

    return MOI.ScalarNonlinearFunction(:-, Any[
        sum(2.0 * mesh.points_mesh.differentiation[q,k] * vars[i][k] for k in 1:n_p_dif),
        MOI.ScalarNonlinearFunction(:*, Any[
            Δt,
            transcribe_dyn_fun(
                dif_fun.dyn_fun, i, q, phase_vars, dyn_var_vars, dif_dyn_vars, mesh,
            ),
        ]),
    ])
end