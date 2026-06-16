function two_link_robot_arm(
    model::Interesso.Optimizer;
    starts::Interesso.WSS = Interesso.WSS{DOI.AbstractDynamicSolution}()
)

    MOI.empty!(model)

    # Keep the same shorthand used in orbit_raising/bang_bang
    NDF = DOI.NonlinearDynamicFunction

    # ---------------------------
    # Time as a phase (free tf)
    # ---------------------------
    t = DOI.add_phase(model)
    MOI.add_constraint(model, DOI.Initial(t), MOI.EqualTo(0.0))

    # ---------------------------
    # Controls
    # ---------------------------
    @control(model, u1, t)
    @control(model, u2, t)
    MOI.add_constraint(model, u1, MOI.Interval(-1.0, 1.0))
    MOI.add_constraint(model, u2, MOI.Interval(-1.0, 1.0))

    # ---------------------------
    # States
    # ---------------------------
    @variable(model, x1, t)
    @variable(model, x2, t)
    @variable(model, x3, t)
    @variable(model, x4, t)

    # Boundary conditions (from Two_Link_Robot_Arm.jl)
    MOI.add_constraint(model, DOI.Initial(x1), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Initial(x2), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Initial(x3), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Initial(x4), MOI.EqualTo(0.0))

    MOI.add_constraint(model, DOI.Final(x1), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Final(x2), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Final(x3), MOI.EqualTo(0.5))
    MOI.add_constraint(model, DOI.Final(x4), MOI.EqualTo(0.522))

    # ---------------------------
    # Dynamics (exactly as Tapir F)
    # ---------------------------
    sx3 = NDF(:sin, [x3], t)
    cx3 = NDF(:cos, [x3], t)

    x1_sq = NDF(:^, [x1, 2.0], t)
    x2_sq = NDF(:^, [x2, 2.0], t)
    sx3_sq = NDF(:^, [sx3, 2.0], t)

    den = NDF(:+, [31.0 / 36.0, NDF(:*, [9.0 / 4.0, sx3_sq], t)], t)

    u1_minus_u2 = NDF(:-, [u1, u2], t)

    # x1dot = ( sin(x3)*(9/4*cos(x3)*x1^2) + 2*x2^2 + 4/3*(u1-u2) - 3/2*cos(x3)*u2 ) / den
    cosx3_x1sq = NDF(:*, [cx3, x1_sq], t)
    term1 = NDF(:*, [sx3, NDF(:*, [9.0 / 4.0, cosx3_x1sq], t)], t)
    term2 = NDF(:*, [2.0, x2_sq], t)
    term3 = NDF(:*, [4.0 / 3.0, u1_minus_u2], t)
    term4 = NDF(:*, [-3.0 / 2.0, NDF(:*, [cx3, u2], t)], t)
    rhs1 = NDF(:/, [NDF(:+, [term1, term2, term3, term4], t), den], t)

    MOI.add_constraint(
        model,
        DOI.ExplicitDifferentialFunction(x1, rhs1),
        MOI.EqualTo(0.0),
    )

    # x2dot = - ( sin(x3)*(9/4*cos(x3)*x2^2) + 7/2*x1^2 - 7/3*u2 + 3/2*cos(x3)*(u1-u2) ) / den
    cosx3_x2sq = NDF(:*, [cx3, x2_sq], t)
    t2_1 = NDF(:*, [sx3, NDF(:*, [9.0 / 4.0, cosx3_x2sq], t)], t)
    t2_2 = NDF(:*, [7.0 / 2.0, x1_sq], t)
    t2_3 = NDF(:*, [-7.0 / 3.0, u2], t)
    t2_4 = NDF(:*, [3.0 / 2.0, NDF(:*, [cx3, u1_minus_u2], t)], t)

    num2 = NDF(:+, [t2_1, t2_2, t2_3, t2_4], t)
    rhs2 = NDF(:/, [NDF(:*, [-1.0, num2], t), den], t)

    MOI.add_constraint(
        model,
        DOI.ExplicitDifferentialFunction(x2, rhs2),
        MOI.EqualTo(0.0),
    )

    # x3dot = x2 - x1
    MOI.add_constraint(
        model,
        DOI.ExplicitDifferentialFunction(x3, NDF(:-, [x2, x1], t)),
        MOI.EqualTo(0.0),
    )

    # x4dot = x1
    MOI.add_constraint(
        model,
        DOI.ExplicitDifferentialFunction(x4, NDF(:+, Any[x1], t)),
        MOI.EqualTo(0.0),
    )

    # ---------------------------
    # Objective (single integral)
    # ---------------------------
    # tf + ∫ 0.01*(u1^2+u2^2) dt  ==  ∫ (1 + 0.01*(u1^2+u2^2)) dt
    u1_sq = NDF(:^, [u1, 2.0], t)
    u2_sq = NDF(:^, [u2, 2.0], t)
    quad = NDF(:*, [0.01, NDF(:+, [u1_sq, u2_sq], t)], t)
    integrand = NDF(:+, [1.0, quad], t)

    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    obj_fun = DOI.MultiPhaseIntegral([integrand])
    MOI.set(model, MOI.ObjectiveFunction{typeof(obj_fun)}(), obj_fun)

    Interesso.warmstart!(model, starts)

    return nothing
end