"""
    bicycle_control_LT(model::Interesso.Optimizer; param=VehicleParams(), starts=Interesso.WSS{DOI.AbstractDynamicSolution}())

Interesso/DOI formulation of a constant-curvature vehicle minimum-time problem.

State (dynamic variables):
  s      distance along reference [m]
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
    model::Interesso.Optimizer;
    starts::Interesso.WSS = Interesso.WSS{DOI.AbstractDynamicSolution}(),
)

    MOI.empty!(model)

    param = VehicleParams()
    # -----------------------------
    # Phase (time)
    # -----------------------------
    t = DOI.add_phase(model)
    t_0 = 0.0
    t_f = 5.0
    MOI.add_constraint(model, DOI.Initial(t), MOI.EqualTo(t_0))
    MOI.add_constraint(model, DOI.Final(t), MOI.EqualTo(t_f))

    # -----------------------------
    # Controls / algebraic variables
    # -----------------------------
    @variable(model, δ, t)
    MOI.add_constraint(model, δ, MOI.Interval(-π / 6, π / 6))

    @variable(model, u_T, t) # throttle
    MOI.add_constraint(model, u_T, MOI.Interval(0.0, 1.0))

    @variable(model, u_B, t) # brake
    MOI.add_constraint(model, u_B, MOI.Interval(0.0, 1.0))

    # No simultaneous throttle + brake
    MOI.add_constraint(model, NDF(:*, [u_T, u_B], t), MOI.EqualTo(0.0))

    @variable(model, κ_fx, t)
    @variable(model, κ_rx, t)
    @variable(model, κ_fy, t)
    @variable(model, κ_ry, t)

    κ_lim = 1.0
    MOI.add_constraint(model, κ_fx, MOI.Interval(-κ_lim, κ_lim))
    MOI.add_constraint(model, κ_rx, MOI.Interval(-κ_lim, κ_lim))
    MOI.add_constraint(model, κ_fy, MOI.Interval(-κ_lim, κ_lim))
    MOI.add_constraint(model, κ_ry, MOI.Interval(-κ_lim, κ_lim))

    # -----------------------------
    # States
    # -----------------------------
    @variable(model, s, t)

    @variable(model, e_y, t)

    @variable(model, v_x, t)

    @variable(model, v_y, t)

    @variable(model, ξ, t)

    @variable(model, dψ, t)

    @variable(model, ω_f, t)
    MOI.add_constraint(model, ω_f, MOI.GreaterThan(0.0))

    @variable(model, ω_r, t)
    MOI.add_constraint(model, ω_r, MOI.GreaterThan(0.0))

    # -----------------------------
    # Boundary conditions
    # -----------------------------
    vxi = 5.0
    scale = 0.5 * (param.J_wf + param.J_wr) / param.J_zz
    MOI.add_constraint(model, DOI.Initial(s), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Initial(v_x), MOI.EqualTo(vxi))
    MOI.add_constraint(model, DOI.Initial(v_y), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Initial(dψ), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Initial(ω_f), MOI.EqualTo(vxi * scale / param.R_e))
    MOI.add_constraint(model, DOI.Initial(ω_r), MOI.EqualTo(vxi * scale / param.R_e))

    # -----------------------------
    # Common nonlinear expressions
    # -----------------------------
    cosδ = NDF(:cos, [δ], t)
    sinδ = NDF(:sin, [δ], t)
    cosξ = NDF(:cos, [ξ], t)
    sinξ = NDF(:sin, [ξ], t)

    vx2 = NDF(:^, [v_x, 2.0], t)

    # Aero and normal loads
    F_d  = NDF(:*, [param.kFd,  vx2], t)
    F_lf = NDF(:*, [param.kFlf, vx2], t)
    F_lr = NDF(:*, [param.kFlr, vx2], t)

    F_zf = NDF(:+, [param.kWf, F_lf], t)
    F_zr = NDF(:+, [param.kWr, F_lr], t)

    # Contact patch velocities
    v_yf = NDF(:+, [v_y, NDF(:*, [param.l_f, dψ], t)], t)
    v_yr = NDF(:-, [v_y, NDF(:*, [param.l_r, dψ], t)], t)
    v_fx = NDF(:+, [NDF(:*, [v_x, cosδ], t), NDF(:*, [v_yf, sinδ], t)], t)
    v_fy = NDF(:+, [NDF(:*, [-1.0, NDF(:*, [v_x, sinδ], t)], t), NDF(:*, [v_yf, cosδ], t)], t)

    # -----------------------------
    # Algebraic (path equality) constraints: slip definitions
    #   (κ+1)*v - R_e*ω = 0 ;  κ*v + v_lat = 0
    # -----------------------------
    MOI.add_constraint(
        model,
        NDF(:-, [NDF(:*, [NDF(:+, [κ_fx, 1.0], t), v_fx], t), NDF(:*, [(param.R_e / scale), ω_f], t)], t),
        MOI.EqualTo(0.0),
    )

    MOI.add_constraint(
        model,
        NDF(:-, [NDF(:*, [NDF(:+, [κ_rx, 1.0], t), v_x], t), NDF(:*, [(param.R_e / scale), ω_r], t)], t),
        MOI.EqualTo(0.0),
    )

    MOI.add_constraint(
        model,
        NDF(:+, [NDF(:*, [κ_fy, v_fx], t), v_fy], t),
        MOI.EqualTo(0.0),
    )

    MOI.add_constraint(
        model,
        NDF(:+, [NDF(:*, [κ_ry, v_x], t), v_yr], t),
        MOI.EqualTo(0.0),
    )

    F_xf = NDF(:*, [F_zf, 20.0, κ_fx], t)
    F_xr = NDF(:*, [F_zr, 20.0, κ_rx], t)
    F_yf = NDF(:*, [F_zf, 15.0, κ_fy], t)
    F_yr = NDF(:*, [F_zr, 15.0, κ_ry], t)

    # -----------------------------
    # Differential equations
    # -----------------------------
    s_denom = NDF(:-, [1.0, NDF(:*, [param.κ_c, e_y], t)], t)
    s_nom = NDF(:-, [NDF(:*, [v_x, cosξ], t), NDF(:*, [v_y, sinξ], t)], t)
    ds = NDF(:/, [s_nom, s_denom], t)

    de = NDF(:+, [NDF(:*, [v_x, sinξ], t), NDF(:*, [v_y, cosξ], t)], t)
    dξ = NDF(:-, [dψ, NDF(:*, [param.κ_c, ds], t)], t)

    Fxf_cosδ = NDF(:*, [F_xf, cosδ], t)
    Fyf_sinδ = NDF(:*, [F_yf, sinδ], t)
    Fxf_sinδ = NDF(:*, [F_xf, sinδ], t)
    Fyf_cosδ = NDF(:*, [F_yf, cosδ], t)

    # v_y*dψ + (F_xf*cosδ + F_xr - F_yf*sinδ - F_d)/m
    dv_x = NDF(:+, [
        NDF(:*, [
            NDF(:-, [
                NDF(:+, [Fxf_cosδ, F_xr], t),
                NDF(:+, [Fyf_sinδ, F_d], t)
            ], t),
            (1.0 / param.m),
        ], t),
        NDF(:*, [v_y, dψ], t)
    ], t)

    # -v_x*dψ + (F_xf*sinδ + F_yf*cosδ + F_yr)/m
    dv_y = NDF(:-, [
        NDF(:*, [
            NDF(:+, [Fxf_sinδ, Fyf_cosδ, F_yr], t),
            (1.0 / param.m),
        ], t),
        NDF(:*, [v_x, dψ], t)
    ], t)

    # ((F_xf*sinδ + F_yf*cosδ)*l_f - F_yr*l_r)/J_zz
    ddψ = NDF(:*, [
        NDF(:-, [
            NDF(:*, [NDF(:+, [Fxf_sinδ, Fyf_cosδ], t), param.l_f], t),
            NDF(:*, [F_yr, param.l_r], t),
        ], t),
        (1.0 / param.J_zz)
    ], t)

    # (-F_xf*R_e - u_B*B_b*B_kf)/J_wf
    dω_f = NDF(:*, [
        NDF(:-, [
            NDF(:-, [NDF(:*, [F_xf, param.R_e], t)], t),
            NDF(:*, [u_B, NDF(:*, [param.B_b, param.B_kf], t)], t),
        ], t),
        (scale / param.J_wf)
    ], t)

    # (-F_xr*R_e + u_T*T_e*τ_g - u_B*(1-B_b)*B_kr)/J_wr
    dω_r = NDF(:*, [
        NDF(:-, [
            NDF(:-, [
                NDF(:*, [u_T, (param.T_e * param.τ_g)], t),
                NDF(:*, [F_xr, param.R_e], t),
            ], t),
            NDF(:*, [u_B, ((1.0 - param.B_b) * param.B_kr)], t),
        ], t),
        (scale / param.J_wr)
    ], t)

    # Dynamics
    MOI.add_constraint(model, DOI.ExplicitDifferentialFunction(s,   ds),   MOI.EqualTo(0.0))
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
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    obj_fun = DOI.NonlinearBoundaryFunction(:+, [DOI.Final(s)])
    # obj_fun = DOI.MultiPhaseIntegral([NDF(:+, [s], t)])
    MOI.set(model, MOI.ObjectiveFunction{typeof(obj_fun)}(), obj_fun)

    # -----------------------------
    # Warm start (simple defaults)
    # -----------------------------
    if starts == Interesso.WSS{DOI.AbstractDynamicSolution}()
        MOI.set(model, DOI.DynamicVariableStart(), s,   LinearInterpolant(0.0, 300.0, t_0, t_f))
        MOI.set(model, DOI.DynamicVariableStart(), e_y, LinearInterpolant(0.0, 0.0,  t_0, t_f))
        MOI.set(model, DOI.DynamicVariableStart(), v_x, LinearInterpolant(vxi, 100.0, t_0, t_f))
        MOI.set(model, DOI.DynamicVariableStart(), v_y, LinearInterpolant(0.0, 0.0,  t_0, t_f))
        MOI.set(model, DOI.DynamicVariableStart(), ξ,   LinearInterpolant(0.0, 0.0,  t_0, t_f))
        MOI.set(model, DOI.DynamicVariableStart(), dψ,  LinearInterpolant(0.0, 0.0,  t_0, t_f))
        MOI.set(model, DOI.DynamicVariableStart(), ω_f, LinearInterpolant(vxi * scale / param.R_e, 100.0 * scale / param.R_e, t_0, t_f))
        MOI.set(model, DOI.DynamicVariableStart(), ω_r, LinearInterpolant(vxi * scale / param.R_e, 100.0 * scale / param.R_e, t_0, t_f))

        # MOI.set(model, DOI.DynamicVariableStart(), δ,    LinearInterpolant(0.0, 0.0, t_0, t_f))
        # MOI.set(model, DOI.DynamicVariableStart(), u_T,  LinearInterpolant(1.0, 1.0, t_0, t_f))
        # MOI.set(model, DOI.DynamicVariableStart(), u_B,  LinearInterpolant(0.0, 0.0, t_0, t_f))
        MOI.set(model, DOI.DynamicVariableStart(), κ_fx, LinearInterpolant(0.0, 0.0, t_0, t_f))
        MOI.set(model, DOI.DynamicVariableStart(), κ_rx, LinearInterpolant(0.0, 0.0, t_0, t_f))
        MOI.set(model, DOI.DynamicVariableStart(), κ_fy, LinearInterpolant(0.0, 0.0, t_0, t_f))
        MOI.set(model, DOI.DynamicVariableStart(), κ_ry, LinearInterpolant(0.0, 0.0, t_0, t_f))
    else
        Interesso.warmstart!(model, starts)
    end

    return nothing
end
