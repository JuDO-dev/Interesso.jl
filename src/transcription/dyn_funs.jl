# Fallback
function transcribe_dyn_fun(
    ::Any,
    ::Integer,
    ::Integer,
    ::PHS_VARS,
    ::TIME_VAR,
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
    ::TIME_VAR,
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
    ::TIME_VAR,
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
    ::TIME_VAR,
    ::DYN_VAR_VARS,
    ::AbstractSet{DYN_VAR},
    mesh::FixedIntervalsMesh,
)
    return mesh.method_meshes[i].quad_points_mesh.points_alg[q]
end

function transcribe_dyn_fun(
    phase::PHS,
    i::Integer,
    q::Integer,
    phase_vars::PHS_VARS,
    ::TIME_VAR,
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
        t_a = flex_vars[i - 1]
        t_b = flex_vars[i]
    end

    return (0.5 * t_a + 0.5 * t_b) + (0.5 * t_b - 0.5 * t_a) * mesh.method_mesh.quad_points_mesh.points_alg[q]
end

# AbstractInterpolant (fixed mesh only)
function transcribe_dyn_fun(
    p::T,
    i::Integer,
    q::Integer,
    ::PHS_VARS,
    ::TIME_VAR,
    ::DYN_VAR_VARS,
    ::AbstractSet{DYN_VAR},
    mesh::FixedIntervalsMesh,
) where {T<:AbstractInterpolant}
    s_val = mesh.method_meshes[i].quad_points_mesh.points_alg[q]
    return p(s_val)
end

# Dynamic Variable
function transcribe_dyn_fun(
    dyn_var::DYN_VAR,
    i::Integer,
    q::Integer,
    ::PHS_VARS,
    ::TIME_VAR,
    dyn_var_vars::DYN_VAR_VARS,
    dif_dyn_vars::AbstractSet{DYN_VAR},
    mesh::FixedIntervalsMesh,
)

    if dyn_var in dif_dyn_vars
        points_quad = mesh.method_meshes[i].interpolant.interpolant_dif * dyn_var_vars[dyn_var][i]
    else
        points_quad = mesh.method_meshes[i].interpolant.interpolant_alg * dyn_var_vars[dyn_var][i]
    end
    return points_quad[q]
end

function transcribe_dyn_fun(
    dyn_var::DYN_VAR,
    i::Integer,
    q::Integer,
    ::PHS_VARS,
    ::TIME_VAR,
    dyn_var_vars::DYN_VAR_VARS,
    dif_dyn_vars::AbstractSet{DYN_VAR},
    mesh::FlexibleIntervalsMesh,
)
    """
    need to do interpolations to dyn_var_vars, if least-square
    interpolant dyn_var_vars[dyn_var][i], where should be a vector of length(points.alg), to become a vector of quad
    need to save the interpolation matrix somewhere
    """
    if dyn_var in dif_dyn_vars
        points_quad = mesh.method_mesh.interpolant.interpolant_dif * dyn_var_vars[dyn_var][i]
    else
        points_quad = mesh.method_mesh.interpolant.interpolant_alg * dyn_var_vars[dyn_var][i]
    end
    return points_quad[q]
end


# Derivative of Dynamic Variable
function transcribe_dyn_fun(
    derivative::DOI.Derivative{DYN_VAR},
    i::Integer,
    q::Integer,
    ::PHS_VARS,
    time_var::TIME_VAR,
    dyn_var_vars::DYN_VAR_VARS,
    ::AbstractSet{DYN_VAR},
    mesh::FixedIntervalsMesh,
)
    vars = dyn_var_vars[derivative.dyn_fun]
    n_p_dif = get_points_dif_length(mesh)

    differentiation =
        mesh.method_meshes[i].interpolant.interpolant_dif *
        mesh.points_meshes[i].differentiation

    numer = sum(differentiation[q, k] * vars[i][k] for k in 1:n_p_dif)
    return MOI.ScalarNonlinearFunction(:/, Any[numer, time_var])
end

function transcribe_dyn_fun(
    derivative::DOI.Derivative{DYN_VAR},
    i::Integer,
    q::Integer,
    phase_vars::PHS_VARS,
    time_var::TIME_VAR,
    dyn_var_vars::DYN_VAR_VARS,
    ::AbstractSet{DYN_VAR},
    mesh::FlexibleIntervalsMesh,
)
    vars = dyn_var_vars[derivative.dyn_fun]
    n_p_dif = get_points_dif_length(mesh)

    differentiation =
        mesh.method_mesh.interpolant.interpolant_dif *
        mesh.points_mesh.differentiation

    numer = sum(2.0 * differentiation[q, k] * vars[i][k] for k in 1:n_p_dif)

    flex_vars = phase_vars[DOI.phase_index(derivative.dyn_fun)]
    n_h = get_intervals_length(mesh)
    t_0 = mesh.fixed.points_meshes[1].t_a
    t_f = mesh.fixed.points_meshes[end].t_b

    if i == 1
        Δt = 1.0 * flex_vars[1] - t_0
    elseif i == n_h
        Δt = t_f - 1.0 * flex_vars[end]
    else
        Δt = 1.0 * flex_vars[i] - 1.0 * flex_vars[i - 1]
    end

    denom = MOI.ScalarNonlinearFunction(:*, Any[time_var, Δt])
    return MOI.ScalarNonlinearFunction(:/, Any[numer, denom])
end

# Nonlinear Dynamic Function
function transcribe_dyn_fun(
    nl_dyn_fun::NDF,
    i::Integer,
    q::Integer,
    phase_vars::PHS_VARS,
    time_var::TIME_VAR,
    dyn_var_vars::DYN_VAR_VARS,
    dif_dyn_vars::AbstractSet{DYN_VAR},
    mesh::AbstractIntervalsMesh,
)
    return MOI.ScalarNonlinearFunction(
        nl_dyn_fun.head,
        [transcribe_dyn_fun(
            arg, i, q, phase_vars, time_var, dyn_var_vars, dif_dyn_vars, mesh,
        ) for arg in nl_dyn_fun.args],
    )
end

# Explicit Differential Function
function transcribe_dyn_fun(
    dif_fun::DIF_FUN,
    i::Integer,
    q::Integer,
    phase_vars::PHS_VARS,
    time_var::TIME_VAR,
    dyn_var_vars::DYN_VAR_VARS,
    dif_dyn_vars::AbstractSet{DYN_VAR},
    mesh::FixedIntervalsMesh,
)
    vars = dyn_var_vars[dif_fun.dyn_var]
    n_p_dif = get_points_dif_length(mesh)

    differentiation = mesh.method_meshes[i].interpolant.interpolant_dif * mesh.points_meshes[i].differentiation

    return MOI.ScalarNonlinearFunction(:-, Any[
        sum(differentiation[q,k] * vars[i][k] for k in 1:n_p_dif),
        MOI.ScalarNonlinearFunction(:*, Any[
            time_var,
            transcribe_dyn_fun(
                dif_fun.dyn_fun, i, q, phase_vars, time_var, dyn_var_vars, dif_dyn_vars, mesh,
            ),
        ]),
    ])
end

function transcribe_dyn_fun(
    dif_fun::DIF_FUN,
    i::Integer,
    q::Integer,
    phase_vars::PHS_VARS,
    time_var::TIME_VAR,
    dyn_var_vars::DYN_VAR_VARS,
    dif_dyn_vars::AbstractSet{DYN_VAR},
    mesh::FlexibleIntervalsMesh,
)
    vars = dyn_var_vars[dif_fun.dyn_var]
    n_p_dif = get_points_dif_length(mesh)

    Δt = get_time_length(phase_vars, i, DOI.phase_index(dif_fun), mesh)

    """
    differentiation matrix here should be equivalent to mesh.Dx * mesh.QX in Tapir

    where Dx is a square matrix and QX is a matrix expanding X to the number of Q, so (num_quad, num_dif)
    """
    differentiation = mesh.method_mesh.interpolant.interpolant_dif * mesh.points_mesh.differentiation

    return MOI.ScalarNonlinearFunction(:-, Any[
        sum(2.0 * differentiation[q,k] * vars[i][k] for k in 1:n_p_dif),
        MOI.ScalarNonlinearFunction(:*, Any[
            time_var,
            Δt,
            transcribe_dyn_fun(
                dif_fun.dyn_fun, i, q, phase_vars, time_var, dyn_var_vars, dif_dyn_vars, mesh,
            ),
        ]),
    ])
end

function transcribe_dyn_fun(
    nl_fun::T,
    i::Integer,
    q::Integer,
    phase_vars::PHS_VARS,
    time_var::TIME_VAR,
    dyn_var_vars::DYN_VAR_VARS,
    dif_dyn_vars::AbstractSet{DYN_VAR},
    mesh::AbstractIntervalsMesh,
    scaling::Float64
) where {T<:Union{NDF,DIF_FUN}}

    f = transcribe_dyn_fun(nl_fun, i, q, phase_vars, time_var, dyn_var_vars, dif_dyn_vars, mesh)
    if scaling == 1.0
        return f
    else
        return MOI.ScalarNonlinearFunction(:*, [scaling, f])
    end
end

function get_time_length(
    ::PHS_VARS,
    i::Integer,
    ::PHS,
    mesh::FixedIntervalsMesh,
)
    t_0 = mesh.points_meshes[i].t_a
    t_f = mesh.points_meshes[i].t_b

    return t_f - t_0
end

function get_time_length(
    phase_vars::PHS_VARS,
    i::Integer,
    phase::PHS,
    mesh::FlexibleIntervalsMesh,
)
    flex_vars = phase_vars[phase]
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

    return Δt
end