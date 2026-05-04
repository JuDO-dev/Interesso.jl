"""
    bicycle_control_LT(model::Interesso.Optimizer; param=VehicleParams(), starts=Interesso.WSS{DOI.AbstractDynamicSolution}())

Interesso/DOI formulation of a constant-curvature vehicle minimum-time problem.

State (dynamic variables):
  t      time [s]
  e_y    lateral deviation from reference [m]
  v_x    body-frame longitudinal speed [m/s]
  v_y    body-frame lateral speed [m/s]
  ξ      heading error to reference [rad]
  dψ      yaw rate [rad/s]
  ω_f    front wheel angular rate [rad/s]
  ω_r    rear  wheel angular rate [rad/s]

Controls / algebraic variables (dynamic variables):
  δ           steering angle [rad]
  u_T, u_B    throttle and brake [0–1]
  κ_fx, κ_rx  longitudinal slip ratios front/rear
  κ_fy, κ_ry  lateral slip ratios front/rear

Slip variables are enforced through algebraic equalities (no divisions).
"""

include(joinpath(@__DIR__, "vehicle_param.jl"))

function linear_bicycle(
    model::Interesso.Optimizer,
    trackfile::String;
    starts::Interesso.WSS = Interesso.WSS{DOI.AbstractDynamicSolution}(),
)

    MOI.empty!(model)

    param = VehicleParams()

    # Track splines
    sref, _, _, _, κref, nlref, nrref, _ = getTrack(trackfile)
    κ_c = CubicInterpolant(sref, κref)
    n_l = CubicInterpolant(sref, nlref)
    n_r = CubicInterpolant(sref, nrref)

    # -----------------------------
    # Phase (time)
    # -----------------------------
    s = DOI.add_phase(model)
    s_0 = sref[1]
    s_f = sref[end]
    MOI.add_constraint(model, DOI.Initial(s), MOI.EqualTo(s_0))
    MOI.add_constraint(model, DOI.Final(s), MOI.EqualTo(s_f))

    # -----------------------------
    # Controls / algebraic variables
    # -----------------------------
    @variable(model, δ, s)
    MOI.add_constraint(model, δ, MOI.Interval(-π / 6, π / 6))

    @variable(model, u_T, s) # throttle
    MOI.add_constraint(model, u_T, MOI.Interval(0.0, 1.0))

    @variable(model, u_B, s) # brake
    MOI.add_constraint(model, u_B, MOI.Interval(0.0, 1.0))

    # No simultaneous throttle + brake
    # MOI.add_constraint(model, NDF(:*, [u_T, u_B], s), MOI.EqualTo(0.0))

    @variable(model, κ_fx, s)
    @variable(model, κ_rx, s)
    @variable(model, κ_fy, s)
    @variable(model, κ_ry, s)

    κ_lim = 1.0
    MOI.add_constraint(model, κ_fx, MOI.Interval(-κ_lim, κ_lim))
    MOI.add_constraint(model, κ_rx, MOI.Interval(-κ_lim, κ_lim))
    MOI.add_constraint(model, κ_fy, MOI.Interval(-κ_lim, κ_lim))
    MOI.add_constraint(model, κ_ry, MOI.Interval(-κ_lim, κ_lim))

    # -----------------------------
    # States
    # -----------------------------
    @variable(model, t, s)

    @variable(model, e_y, s)

    @variable(model, v_x, s)

    @variable(model, v_y, s)

    @variable(model, ξ, s)

    @variable(model, dψ, s)

    @variable(model, ω_f, s)
    MOI.add_constraint(model, ω_f, MOI.GreaterThan(0.0))

    @variable(model, ω_r, s)
    MOI.add_constraint(model, ω_r, MOI.GreaterThan(0.0))

    # -----------------------------
    # Boundary conditions
    # -----------------------------
    vxi = 5.0
    scale = 0.5 * (param.J_wf + param.J_wr) / param.J_zz
    MOI.add_constraint(model, DOI.Initial(t), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Initial(v_x), MOI.EqualTo(vxi))
    MOI.add_constraint(model, DOI.Initial(v_y), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Initial(dψ), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Initial(ω_f), MOI.EqualTo(vxi * scale / param.R_e))
    MOI.add_constraint(model, DOI.Initial(ω_r), MOI.EqualTo(vxi * scale / param.R_e))

    MOI.add_constraint(model, NDF(:-, [n_l, e_y], s), MOI.GreaterThan(0.0))
    MOI.add_constraint(model, NDF(:+, [n_r, e_y], s), MOI.GreaterThan(0.0))

    # -----------------------------
    # Common nonlinear expressions
    # -----------------------------
    cosδ = NDF(:cos, [δ], s)
    sinδ = NDF(:sin, [δ], s)
    cosξ = NDF(:cos, [ξ], s)
    sinξ = NDF(:sin, [ξ], s)

    vx2 = NDF(:^, [v_x, 2.0], s)

    # Aero and normal loads
    F_d  = NDF(:*, [param.kFd,  vx2], s)
    F_lf = NDF(:*, [param.kFlf, vx2], s)
    F_lr = NDF(:*, [param.kFlr, vx2], s)

    F_zf = NDF(:+, [param.kWf, F_lf], s)
    F_zr = NDF(:+, [param.kWr, F_lr], s)

    # Contact patch velocities
    v_yf = NDF(:+, [v_y, NDF(:*, [param.l_f, dψ], s)], s)
    v_yr = NDF(:-, [v_y, NDF(:*, [param.l_r, dψ], s)], s)
    v_fx = NDF(:+, [NDF(:*, [v_x, cosδ], s), NDF(:*, [v_yf, sinδ], s)], s)
    v_fy = NDF(:+, [NDF(:*, [-1.0, NDF(:*, [v_x, sinδ], s)], s), NDF(:*, [v_yf, cosδ], s)], s)

    # -----------------------------
    # Algebraic (path equality) constraints: slip definitions
    #   (κ+1)*v - R_e*ω = 0 ;  κ*v + v_lat = 0
    # -----------------------------
    MOI.add_constraint(
        model,
        NDF(:-, [NDF(:*, [NDF(:+, [κ_fx, 1.0], s), v_fx], s), NDF(:*, [(param.R_e / scale), ω_f], s)], s),
        MOI.EqualTo(0.0),
        1e-2
    )

    MOI.add_constraint(
        model,
        NDF(:-, [NDF(:*, [NDF(:+, [κ_rx, 1.0], s), v_x], s), NDF(:*, [(param.R_e / scale), ω_r], s)], s),
        MOI.EqualTo(0.0),
        1e-2
    )

    MOI.add_constraint(
        model,
        NDF(:+, [NDF(:*, [κ_fy, v_fx], s), v_fy], s),
        MOI.EqualTo(0.0),
        1e-2
    )

    MOI.add_constraint(
        model,
        NDF(:+, [NDF(:*, [κ_ry, v_x], s), v_yr], s),
        MOI.EqualTo(0.0),
        1e-2
    )

    nF_xf = NDF(:*, [20.0, κ_fx], s)
    nF_xr = NDF(:*, [20.0, κ_rx], s)
    nF_yf = NDF(:*, [15.0, κ_fy], s)
    nF_yr = NDF(:*, [15.0, κ_ry], s)

    F_xf = NDF(:*, [F_zf, nF_xf], s)
    F_xr = NDF(:*, [F_zr, nF_xr], s)
    F_yf = NDF(:*, [F_zf, nF_yf], s)
    F_yr = NDF(:*, [F_zr, nF_yr], s)

    nF_xf2 = NDF(:^, [nF_xf, 2.0], s)
    nF_xr2 = NDF(:^, [nF_xr, 2.0], s)
    nF_yf2 = NDF(:^, [nF_yf, 2.0], s)
    nF_yr2 = NDF(:^, [nF_yr, 2.0], s)

    MOI.add_constraint(model, NDF(:+, [nF_xf2, nF_yf2], s), MOI.LessThan(1.5))
    MOI.add_constraint(model, NDF(:+, [nF_xr2, nF_yr2], s), MOI.LessThan(1.5))

    # -----------------------------
    # Differential equations
    # -----------------------------
    t_nom = NDF(:-, [1.0, NDF(:*, [κ_c, e_y], s)], s)
    t_denom = NDF(:-, [NDF(:*, [v_x, cosξ], s), NDF(:*, [v_y, sinξ], s)], s)
    dt = NDF(:/, [t_nom, t_denom], s)

    MOI.add_constraint(model, t_nom,   MOI.GreaterThan(0.0))
    MOI.add_constraint(model, t_denom, MOI.GreaterThan(0.0))

    de = NDF(:*, [NDF(:+, [NDF(:*, [v_x, sinξ], s), NDF(:*, [v_y, cosξ], s)], s), dt], s)
    dξ = NDF(:-, [NDF(:*, [dψ, dt], s), κ_c], s)

    Fxf_cosδ = NDF(:*, [F_xf, cosδ], s)
    Fyf_sinδ = NDF(:*, [F_yf, sinδ], s)
    Fxf_sinδ = NDF(:*, [F_xf, sinδ], s)
    Fyf_cosδ = NDF(:*, [F_yf, cosδ], s)

    # (v_y*dψ + (F_xf*cosδ + F_xr - F_yf*sinδ - F_d)/m) * dt
    dv_x = NDF(:*, [
        NDF(:+, [
            NDF(:*, [
                NDF(:-, [
                    NDF(:+, [Fxf_cosδ, F_xr], s),
                    NDF(:+, [Fyf_sinδ, F_d], s)
                ], s),
                (1.0 / param.m)
            ], s),
            NDF(:*, [v_y, dψ], s)
        ], s),
        dt
    ], s)

    # (-v_x*dψ + (F_xf*sinδ + F_yf*cosδ + F_yr)/m) * dt
    dv_y = NDF(:*, [
        NDF(:-, [
            NDF(:*, [
                NDF(:+, [Fxf_sinδ, Fyf_cosδ, F_yr], s),
                (1.0 / param.m)
            ], s),
            NDF(:*, [v_x, dψ], s)
        ], s),
        dt
    ], s)

    # ((F_xf*sinδ + F_yf*cosδ)*l_f - F_yr*l_r)/J_zz * dt
    ddψ = NDF(:*, [
        NDF(:-, [
            NDF(:*, [NDF(:+, [Fxf_sinδ, Fyf_cosδ], s), param.l_f], s),
            NDF(:*, [F_yr, param.l_r], s),
        ], s),
        (1.0 / param.J_zz),
        dt
    ], s)

    # (-F_xf*R_e - u_B*B_b*B_kf)/J_wf * dt
    dω_f = NDF(:*, [
        NDF(:-, [
            NDF(:-, [NDF(:*, [F_xf, param.R_e], s)], s),
            NDF(:*, [u_B, NDF(:*, [param.B_b, param.B_kf], s)], s),
        ], s),
        (scale / param.J_wf),
        dt
    ], s)

    # (-F_xr*R_e + u_T*T_e*τ_g - u_B*(1-B_b)*B_kr)/J_wr * dt
    dω_r = NDF(:*, [
        NDF(:-, [
            NDF(:-, [
                NDF(:*, [u_T, (param.T_e * param.τ_g)], s),
                NDF(:*, [F_xr, param.R_e], s),
            ], s),
            NDF(:*, [u_B, ((1.0 - param.B_b) * param.B_kr)], s),
        ], s),
        (scale / param.J_wr),
        dt
    ], s)

    # Dynamics
    MOI.add_constraint(model, DOI.ExplicitDifferentialFunction(t,   dt),   MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.ExplicitDifferentialFunction(e_y, de),   MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.ExplicitDifferentialFunction(v_x, dv_x), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.ExplicitDifferentialFunction(v_y, dv_y), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.ExplicitDifferentialFunction(ξ,   dξ),   MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.ExplicitDifferentialFunction(dψ,  ddψ),  MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.ExplicitDifferentialFunction(ω_f, dω_f), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.ExplicitDifferentialFunction(ω_r, dω_r), MOI.EqualTo(0.0))

    # -----------------------------
    # Objective: minimum final time
    # -----------------------------
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    obj_fun = DOI.NonlinearBoundaryFunction(:+, [DOI.Final(t)])
    # obj_fun = DOI.MultiPhaseIntegral([NDF(:+, [t], s)])
    MOI.set(model, MOI.ObjectiveFunction{typeof(obj_fun)}(), obj_fun)

    # -----------------------------
    # Warm start (simple defaults)
    # -----------------------------
    if starts == Interesso.WSS{DOI.AbstractDynamicSolution}()
        MOI.set(model, DOI.DynamicVariableStart(), t,   LinearInterpolant(0.1, 10.0, s_0, s_f))
        MOI.set(model, DOI.DynamicVariableStart(), e_y, LinearInterpolant(0.0, 0.0,  s_0, s_f))
        MOI.set(model, DOI.DynamicVariableStart(), v_x, LinearInterpolant(0.1, 50.0, s_0, s_f))
        MOI.set(model, DOI.DynamicVariableStart(), v_y, LinearInterpolant(0.0, 0.0,  s_0, s_f))
        MOI.set(model, DOI.DynamicVariableStart(), ξ,   LinearInterpolant(0.0, 0.0,  s_0, s_f))
        MOI.set(model, DOI.DynamicVariableStart(), dψ,  LinearInterpolant(0.0, 0.0,  s_0, s_f))
        MOI.set(model, DOI.DynamicVariableStart(), ω_f, LinearInterpolant(0.1, 50.0, s_0, s_f))
        MOI.set(model, DOI.DynamicVariableStart(), ω_r, LinearInterpolant(0.1, 50.0, s_0, s_f))

        # MOI.set(model, DOI.DynamicVariableStart(), δ,    LinearInterpolant(0.0, 0.0, s_0, s_f))
        # MOI.set(model, DOI.DynamicVariableStart(), u_T,  LinearInterpolant(1.0, 1.0, s_0, s_f))
        # MOI.set(model, DOI.DynamicVariableStart(), u_B,  LinearInterpolant(0.0, 0.0, s_0, s_f))
        MOI.set(model, DOI.DynamicVariableStart(), κ_fx, LinearInterpolant(0.0, 0.0, s_0, s_f))
        MOI.set(model, DOI.DynamicVariableStart(), κ_rx, LinearInterpolant(0.0, 0.0, s_0, s_f))
        MOI.set(model, DOI.DynamicVariableStart(), κ_fy, LinearInterpolant(0.0, 0.0, s_0, s_f))
        MOI.set(model, DOI.DynamicVariableStart(), κ_ry, LinearInterpolant(0.0, 0.0, s_0, s_f))
    else
        Interesso.warmstart!(model, starts)
    end

    return nothing
end
