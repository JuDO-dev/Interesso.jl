""" 
    bicycle_control_LT(model::Interesso.Optimizer; param=VehicleParams(), starts=Interesso.WSS{DOI.AbstractDynamicSolution}())

Interesso/DOI formulation of a constant-curvature vehicle minimum-time problem.
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
    @control(model, δ, t)
    MOI.add_constraint(model, δ, MOI.Interval(-π / 6, π / 6))

    @control(model, u_T, t) # throttle
    MOI.add_constraint(model, u_T, MOI.Interval(0.0, 1.0))

    @control(model, u_B, t) # brake
    MOI.add_constraint(model, u_B, MOI.Interval(0.0, 1.0))

    # No simultaneous throttle + brake
    # MOI.add_constraint(model, NDF(:*, [u_T, u_B], t), MOI.EqualTo(0.0))

    # -----------------------------
    # States
    # -----------------------------
    @variable(model, s, t)

    @variable(model, e_y, t)

    @variable(model, v_x, t)

    @variable(model, v_y, t)

    @variable(model, ξ, t)

    @variable(model, dψ, t)

    # -----------------------------
    # Boundary conditions
    # -----------------------------
    vxi = 5.0

    MOI.add_constraint(model, DOI.Initial(s), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Initial(v_x), MOI.EqualTo(vxi))
    MOI.add_constraint(model, DOI.Initial(v_y), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Initial(dψ), MOI.EqualTo(0.0))

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

    # Wheel-center lateral velocities (single-track)
    v_yf = NDF(:+, [v_y, NDF(:*, [param.l_f, dψ], t)], t)
    v_yr = NDF(:-, [v_y, NDF(:*, [param.l_r, dψ], t)], t)
    κ_fy = NDF(:-, [δ, NDF(:/, [v_yf, v_x], t)], t)
    κ_ry = NDF(:-, [NDF(:/, [v_yr, v_x], t)], t)

    # -----------------------------
    # Direct-force tyre model (NO ω, NO κ):
    #   - Longitudinal: throttle/brake mapped to wheel torques / radius
    #   - Lateral: small-angle slip angles (α_f, α_r)
    #   - Combined limit: friction circle
    # -----------------------------
    # Longitudinal forces (RWD traction, fixed brake bias param.B_b)
    F_xf = NDF(:*, [u_B, (- param.B_b * param.B_kf / param.R_e)], t)
    F_xr = NDF(:-, [
        NDF(:*, [u_T, (param.T_e * param.τ_g / param.R_e)], t),
        NDF(:*, [u_B, ((1.0 - param.B_b) * param.B_kr / param.R_e)], t)
    ], t)    

    # Linear lateral forces scaled by normal load (reuse the same gains as the original κ-based model)
    F_yf = NDF(:*, [F_zf, 15.0, κ_fy], t)
    F_yr = NDF(:*, [F_zr, 15.0, κ_ry], t)

    # Friction circle (set μ here; if your VehicleParams has μ_f/μ_r fields, replace accordingly)
    μf2 = 1.2^2
    μr2 = 1.2^2

    # Front: F_xf^2 + F_yf^2 <= (μf*F_zf)^2
    Fx_f2 = NDF(:^, [F_xf, 2.0], t)
    Fy_f2 = NDF(:^, [F_yf, 2.0], t)
    Fz_f2 = NDF(:^, [F_zf, 2.0], t)
    MOI.add_constraint(model, NDF(:-, [NDF(:+, [Fx_f2, Fy_f2], t), NDF(:*, [μf2, Fz_f2], t)], t), MOI.LessThan(0.0))

    # Rear: F_xr^2 + F_yr^2 <= (μr*F_zr)^2
    Fx_r2 = NDF(:^, [F_xr, 2.0], t)
    Fy_r2 = NDF(:^, [F_yr, 2.0], t)
    Fz_r2 = NDF(:^, [F_zr, 2.0], t)
    MOI.add_constraint(model, NDF(:-, [NDF(:+, [Fx_r2, Fy_r2], t), NDF(:*, [μr2, Fz_r2], t)], t), MOI.LessThan(0.0))

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

    # Dynamics
    MOI.add_constraint(model, DOI.ExplicitDifferentialFunction(s,   ds),   MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.ExplicitDifferentialFunction(e_y, de),   MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.ExplicitDifferentialFunction(v_x, dv_x), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.ExplicitDifferentialFunction(v_y, dv_y), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.ExplicitDifferentialFunction(ξ,   dξ),   MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.ExplicitDifferentialFunction(dψ,  ddψ),  MOI.EqualTo(0.0))

    # -----------------------------
    # Objective: minimum final time
    # -----------------------------
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    obj_fun = DOI.NonlinearBoundaryFunction(:+, [DOI.Final(s)])
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
    else
        Interesso.warmstart!(model, starts)
    end

    return nothing
end
