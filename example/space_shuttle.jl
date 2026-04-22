function space_shuttle_reentry(
    model::Interesso.Optimizer;
    starts::Interesso.WSS=Interesso.WSS{DOI.AbstractDynamicSolution}()
)

    MOI.empty!(model)

    # Problem constants
    m_val = 203000 / 32.174
    ρ_0   = 0.002378
    h_r   = 23800.0
    R_e   = 20902900.0
    μ_val = 0.14076539e17
    a_0   = -0.20704
    a_1   = 0.029244
    b_0   = 0.07854
    b_1   = -0.61592e-2
    b_2   = 0.621408e-3
    S_val = 2690.0

    t_0     = 0.0
    t_f_max = 2500.0

    ## Phase (free final time, tf ≤ 2500)
    t = DOI.add_phase(model)
    MOI.add_constraint(model, DOI.Initial(t), MOI.EqualTo(t_0))
    MOI.add_constraint(model, DOI.Final(t),   MOI.LessThan(t_f_max))

    ## States
    @variable(model, scaled_h, t)
    MOI.add_constraint(model, scaled_h, MOI.GreaterThan(0.0))

    @variable(model, θ, t)
    MOI.add_constraint(model, θ, MOI.Interval(deg2rad(-89.0), deg2rad(89.0)))

    @variable(model, Φ, t)

    @variable(model, scaled_v, t)
    MOI.add_constraint(model, scaled_v, MOI.GreaterThan(1e-4)) #scaled by 1e4

    @variable(model, γ, t)
    MOI.add_constraint(model, γ, MOI.Interval(deg2rad(-89.0), deg2rad(89.0)))

    @variable(model, ψ, t)

    ## Controls
    @variable(model, α, t)
    MOI.add_constraint(model, α, MOI.Interval(deg2rad(-90.0), deg2rad(90.0)))

    @variable(model, β, t)
    MOI.add_constraint(model, β, MOI.Interval(deg2rad(-90.0), deg2rad(1.0)))

    ## Boundary Conditions
    MOI.add_constraint(model, DOI.Initial(scaled_h), MOI.EqualTo(2.6)) #scaled by 1e5
    MOI.add_constraint(model, DOI.Final(scaled_h),   MOI.EqualTo(0.8)) #scaled by 1e5
    MOI.add_constraint(model, DOI.Initial(θ), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Initial(Φ), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Initial(scaled_v), MOI.EqualTo(2.56)) #scaled by 1e4
    MOI.add_constraint(model, DOI.Final(scaled_v),   MOI.EqualTo(0.25)) #scaled by 1e4
    MOI.add_constraint(model, DOI.Initial(γ), MOI.EqualTo(deg2rad(-1.0)))
    MOI.add_constraint(model, DOI.Final(γ),   MOI.EqualTo(deg2rad(-5.0)))
    MOI.add_constraint(model, DOI.Initial(ψ), MOI.EqualTo(deg2rad(90.0)))

    ## Intermediate expressions (reused across dynamics)
    h = NDF(:*, [scaled_h, 1e5], t)
    v = NDF(:*, [scaled_v, 1e4], t)

    r_expr   = NDF(:+, [R_e, h], t)
    g_expr   = NDF(:/, [μ_val, NDF(:^, [r_expr, 2.0], t)], t)
    ρ_expr   = NDF(:*, [ρ_0, NDF(:exp, [NDF(:*, [-1.0 / h_r, h], t)], t)], t)
    α_deg    = NDF(:*, [180.0 / pi, α], t)
    cl       = NDF(:+, [a_0, NDF(:*, [a_1, α_deg], t)], t)
    cd       = NDF(:+, [b_0, NDF(:*, [b_1, α_deg], t), NDF(:*, [b_2, NDF(:^, [α_deg, 2.0], t)], t)], t)
    dyn_pres = NDF(:*, [0.5 * S_val, ρ_expr, NDF(:^, [v, 2.0], t)], t)
    L_expr   = NDF(:*, [dyn_pres, cl], t)
    D_expr   = NDF(:*, [dyn_pres, cd], t)

    ## Differential Equations

    # d(scaled_h)/dt = v * sin(γ) / 1e5
    MOI.add_constraint(
        model,
        DOI.ExplicitDifferentialFunction(scaled_h,
            NDF(:*, [1.0 / 1e5, v, NDF(:sin, [γ], t)], t)
        ),
        MOI.EqualTo(0.0),
    )

    # θ̇ = v * cos(γ) * cos(ψ) / r
    MOI.add_constraint(
        model,
        DOI.ExplicitDifferentialFunction(θ,
            NDF(:/, [NDF(:*, [v, NDF(:cos, [γ], t), NDF(:cos, [ψ], t)], t), r_expr], t)
        ),
        MOI.EqualTo(0.0),
    )

    # Φ̇ = v * cos(γ) * sin(ψ) / (r * cos(θ))
    MOI.add_constraint(
        model,
        DOI.ExplicitDifferentialFunction(Φ,
            NDF(:/, [
                NDF(:*, [v, NDF(:cos, [γ], t), NDF(:sin, [ψ], t)], t),
                NDF(:*, [r_expr, NDF(:cos, [θ], t)], t)
            ], t)
        ),
        MOI.EqualTo(0.0),
    )

    # d(scaled_v)/dt = (-D/m - g*sin(γ)) / 1e4
    MOI.add_constraint(
        model,
        DOI.ExplicitDifferentialFunction(scaled_v,
            NDF(:*, [1.0 / 1e4,
                NDF(:+, [
                    NDF(:*, [-1.0 / m_val, D_expr], t),
                    NDF(:*, [-1.0, g_expr, NDF(:sin, [γ], t)], t)
                ], t)
            ], t)
        ),
        MOI.EqualTo(0.0),
    )

    # γ̇ = L*cos(β)/(m*v) + cos(γ)*(v/r - g/v)
    MOI.add_constraint(
        model,
        DOI.ExplicitDifferentialFunction(γ,
            NDF(:+, [
                NDF(:/, [
                    NDF(:*, [L_expr, NDF(:cos, [β], t)], t),
                    NDF(:*, [m_val, v], t)
                ], t),
                NDF(:*, [
                    NDF(:cos, [γ], t),
                    NDF(:-, [NDF(:/, [v, r_expr], t), NDF(:/, [g_expr, v], t)], t)
                ], t)
            ], t)
        ),
        MOI.EqualTo(0.0),
    )

    # ψ̇ = L*sin(β)/(m*v*cos(γ)) + v*cos(γ)*sin(ψ)*sin(θ)/(r*cos(θ))
    MOI.add_constraint(
        model,
        DOI.ExplicitDifferentialFunction(ψ,
            NDF(:+, [
                NDF(:/, [
                    NDF(:*, [L_expr, NDF(:sin, [β], t)], t),
                    NDF(:*, [m_val, v, NDF(:cos, [γ], t)], t)
                ], t),
                NDF(:/, [
                    NDF(:*, [v, NDF(:cos, [γ], t), NDF(:sin, [ψ], t), NDF(:sin, [θ], t)], t),
                    NDF(:*, [r_expr, NDF(:cos, [θ], t)], t)
                ], t)
            ], t)
        ),
        MOI.EqualTo(0.0),
    )

    ## Warm-starts
    if starts == Interesso.WSS{DOI.AbstractDynamicSolution}()
        MOI.set(model, DOI.DynamicVariableStart(), scaled_h, LinearInterpolant(2.6,           0.8,          t_0, t_f_max))
        MOI.set(model, DOI.DynamicVariableStart(), θ,        LinearInterpolant(0.0,            deg2rad(45.0),  t_0, t_f_max))
        MOI.set(model, DOI.DynamicVariableStart(), Φ,        LinearInterpolant(0.0,            deg2rad(50.0),  t_0, t_f_max))
        MOI.set(model, DOI.DynamicVariableStart(), scaled_v, LinearInterpolant(2.56,           0.25,         t_0, t_f_max))
        MOI.set(model, DOI.DynamicVariableStart(), γ, LinearInterpolant(deg2rad(-1.0),  deg2rad(-5.0),  t_0, t_f_max))
        MOI.set(model, DOI.DynamicVariableStart(), ψ, LinearInterpolant(deg2rad(90.0),  deg2rad(-20.0), t_0, t_f_max))
        MOI.set(model, DOI.DynamicVariableStart(), α, LinearInterpolant(0.0,            0.0,            t_0, t_f_max))
        MOI.set(model, DOI.DynamicVariableStart(), β, LinearInterpolant(0.0,            0.0,            t_0, t_f_max))
    else
        Interesso.warmstart!(model, starts)
    end

    ## Objective: maximize final(θ)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    obj_fun = DOI.NonlinearBoundaryFunction(:+, [DOI.Final(θ)])
    MOI.set(model, MOI.ObjectiveFunction{typeof(obj_fun)}(), obj_fun)

    return nothing
end
