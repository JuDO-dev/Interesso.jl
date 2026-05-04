using Plots

include(joinpath(@__DIR__, "../tracks/get_track.jl"))
include(joinpath(@__DIR__, "race_car_plot.jl"))

# ────────────────────────────────────────────────────────────────────
# NDF helpers
# ────────────────────────────────────────────────────────────────────

_sub(a, b, ph) = NDF(:-, [a, b], ph)
_div(a, b, ph) = NDF(:/, [a, b], ph)
_sq(a, ph)     = NDF(:^, [a, 2.0], ph)

function _add(args...)
    @assert length(args) >= 3
    ph = args[end]
    return NDF(:+, collect(args[1:end-1]), ph)
end

function _mul(args...)
    @assert length(args) >= 3
    ph = args[end]
    return NDF(:*, collect(args[1:end-1]), ph)
end


# ────────────────────────────────────────────────────────────────────
# Race-car OCP — spatial formulation
# ────────────────────────────────────────────────────────────────────

function race_car(
    model::Interesso.Optimizer,
    trackfile::String;
    starts::Interesso.WSS = Interesso.WSS{DOI.AbstractDynamicSolution}(),
)
    MOI.empty!(model)

    # Vehicle parameters
    m   = 0.043
    C1  = 0.5
    C2  = 15.5
    Cm1 = 0.28
    Cm2 = 0.05
    Cr0 = 0.011
    Cr2 = 0.006

    # Track splines
    sref, _, _, _, κref, nlref, nrref, _ = getTrack(trackfile)
    κ = CubicInterpolant(sref, κref)
    nl = CubicInterpolant(sref, nlref)
    nr = CubicInterpolant(sref, nrref)

    # Phase: s ∈ [0, L]
    s = DOI.add_phase(model)
    s_0 = sref[1]
    s_f = sref[end]
    MOI.add_constraint(model, DOI.Initial(s), MOI.EqualTo(s_0))
    MOI.add_constraint(model, DOI.Final(s),   MOI.EqualTo(s_f))

    # Controls
    @variable(model, derD, s)
    MOI.add_constraint(model, derD, MOI.Interval(-10.0, 10.0))

    @variable(model, derδ, s)
    MOI.add_constraint(model, derδ, MOI.Interval(-2.0, 2.0))

    # States
    @variable(model, t, s)

    @variable(model, n, s)
    MOI.add_constraint(model, _sub(n, nl, s), MOI.LessThan(0.0))
    MOI.add_constraint(model, _add(n, nr, s), MOI.GreaterThan(0.0))

    @variable(model, α, s)

    @variable(model, v, s)

    @variable(model, D, s)
    MOI.add_constraint(model, D, MOI.Interval(-1.0, 1.0))

    @variable(model, δ, s)
    MOI.add_constraint(model, δ, MOI.Interval(-0.40, 0.40))

    # Initial conditions
    MOI.add_constraint(model, DOI.Initial(t), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Initial(v), MOI.EqualTo(0.5))

    # Nonlinear sub-expressions
    v2      = _sq(v, s)
    C1δ     = _mul(C1, δ, s)
    β       = _add(α, C1δ, s)
    cos_β   = NDF(:cos, [β], s)
    sin_β   = NDF(:sin, [β], s)
    cos_C1δ = NDF(:cos, [C1δ], s)
    sin_C1δ = NDF(:sin, [C1δ], s)

    drive_gain = _sub(Cm1, _mul(Cm2, v, s), s)
    drive_term = _mul(drive_gain, D, s)
    drag_quad  = _mul(-Cr2, v2, s)
    drag_roll  = _mul(-Cr0, NDF(:tanh, [_mul(5.0, v, s)], s), s)
    Fxd = _add(drive_term, drag_quad, drag_roll, s)

    # Spatial dynamics
    one_minus_κn = _sub(1.0, _mul(κ, n, s), s)
    v_cos_β      = _mul(v, cos_β, s)
    dtds         = _div(one_minus_κn, v_cos_β, s)
    v_dtds       = _div(one_minus_κn, cos_β, s)

    MOI.add_constraint(model, one_minus_κn, MOI.GreaterThan(0.0))
    MOI.add_constraint(model, v_cos_β,      MOI.GreaterThan(0.0))

    dnds = _mul(sin_β, v_dtds, s)
    dαds = _sub(_mul(C2, δ, v_dtds, s), κ, s)
    dvds = _mul(Fxd, (1.0 / m), cos_C1δ, dtds, s)
    dDds = _mul(derD, dtds, s)
    dδds = _mul(derδ, dtds, s)

    # Acceleration constraints
    a_long = _mul(Fxd, (1.0 / m), s)
    a_lat  = _add(
        _mul(C2, v2, δ, s),
        _mul(Fxd, (1.0 / m), sin_C1δ, s),
        s,
    )
    MOI.add_constraint(model, a_long, MOI.Interval(-4.0, 4.0))
    MOI.add_constraint(model, a_lat,  MOI.Interval(-4.0, 4.0))

    # Dynamics
    MOI.add_constraint(model, DOI.ExplicitDifferentialFunction(t, dtds), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.ExplicitDifferentialFunction(n, dnds), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.ExplicitDifferentialFunction(α, dαds), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.ExplicitDifferentialFunction(v, dvds), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.ExplicitDifferentialFunction(D, dDds), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.ExplicitDifferentialFunction(δ, dδds), MOI.EqualTo(0.0))

    # Objective
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    obj_fun = DOI.NonlinearBoundaryFunction(:+, [DOI.Final(t)])
    MOI.set(model, MOI.ObjectiveFunction{typeof(obj_fun)}(), obj_fun)

    # Warm start
    if starts == Interesso.WSS{DOI.AbstractDynamicSolution}()
        MOI.set(model, DOI.DynamicVariableStart(), t, LinearInterpolant(0.1, 10.0, s_0, s_f))
        MOI.set(model, DOI.DynamicVariableStart(), n, LinearInterpolant(0.0, 0.0, s_0, s_f))
        MOI.set(model, DOI.DynamicVariableStart(), α, LinearInterpolant(1.0, 1.0, s_0, s_f))
        MOI.set(model, DOI.DynamicVariableStart(), v, LinearInterpolant(1.0, 1.0, s_0, s_f))
        MOI.set(model, DOI.DynamicVariableStart(), D, LinearInterpolant(1.0, 1.0, s_0, s_f))
        MOI.set(model, DOI.DynamicVariableStart(), δ, LinearInterpolant(1.0, 1.0, s_0, s_f))
    else
        Interesso.warmstart!(model, starts)
    end

    return nothing
end
