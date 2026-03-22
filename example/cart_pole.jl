function cart_pole(
    model::Interesso.Optimizer;
    starts::Interesso.WSS=Interesso.WSS{DOI.AbstractDynamicSolution}()
) 

    MOI.empty!(model)

    ## Time as a phase
    t = DOI.add_phase(model)
    t_0 = 0.0
    t_f = 2.0
    MOI.add_constraint(model, DOI.Initial(t), MOI.EqualTo(t_0))
    MOI.add_constraint(model, DOI.Final(t), MOI.EqualTo(t_f))

    # Problem Constants
    g = 9.81
    l = 0.5
    m_1 = 1.0
    m_2 = 0.3
    u_max = 20.0
    r_max = 2.0

    ## Input Dynamic Variable
    @variable(model, u, t)

    ## State Dynamic Variables
    @variable(model, r, t)
    @variable(model, θ, t)
    @variable(model, v, t)
    @variable(model, ω, t)

    ## Inequality constraint
    MOI.add_constraint(model, u, MOI.Interval(-u_max, u_max))
    MOI.add_constraint(model, r, MOI.Interval(0.0, r_max))

    ## Boundary Conditions
    MOI.add_constraint(model, DOI.Initial(r), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Initial(θ), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Initial(v), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Initial(ω), MOI.EqualTo(0.0))

    MOI.add_constraint(model, DOI.Final(r), MOI.EqualTo(1.0))
    MOI.add_constraint(model, DOI.Final(θ), MOI.EqualTo(1.0 * pi))
    MOI.add_constraint(model, DOI.Final(v), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Final(ω), MOI.EqualTo(0.0))

    # Starts
    if starts == Interesso.WSS{DOI.AbstractDynamicSolution}()
        MOI.set(model, DOI.DynamicVariableStart(), r, LinearInterpolant(0.0, 1.0, t_0, t_f))
        MOI.set(model, DOI.DynamicVariableStart(), θ, LinearInterpolant(0.0, 1.0 * pi, t_0, t_f))
        MOI.set(model, DOI.DynamicVariableStart(), u, LinearInterpolant(0.0, 0.0, t_0, t_f))
        MOI.set(model, DOI.DynamicVariableStart(), v, LinearInterpolant(0.0, 0.0, t_0, t_f))
        MOI.set(model, DOI.DynamicVariableStart(), ω, LinearInterpolant(0.0, 0.0, t_0, t_f))
    else
        Interesso.warmstart!(model, starts)
    end

    ## Differential Equations
    sinθ = NDF(:sin, [θ], t)
    cosθ = NDF(:cos, [θ], t)

    num_v = NDF(:+, [
        NDF(:*, [l * m_2, sinθ, NDF(:^, [ω, 2], t)], t),
        u,
        NDF(:*, [m_2 * g, cosθ, sinθ], t)
    ], t)
    den_v = NDF(:+, [
        m_1,
        NDF(:*, [m_2, NDF(:^, [sinθ, 2], t)], t),
    ], t)

    num_ω = NDF(:+, [
        NDF(:*, [-1.0 * l * m_2, cosθ, sinθ, NDF(:^, [ω, 2], t)], t),
        NDF(:*, [-1.0, u, cosθ], t),
        NDF(:*, [-1.0 * (m_1 + m_2) * g, sinθ], t),
    ], t)
    den_ω = NDF(:*, [
        l,
        NDF(:+, [
            m_1,
            NDF(:*, [m_2, NDF(:^, [sinθ, 2], t)], t),
        ], t)
    ], t)

    MOI.add_constraint(
        model,
        DOI.ExplicitDifferentialFunction(
            r,
            NDF(:*, [1.0, v], t),            
        ),
        MOI.EqualTo(0.0),
    )
    MOI.add_constraint(
        model,
        DOI.ExplicitDifferentialFunction(
            v,
            NDF(:/, [num_v, den_v], t),            
        ),
        MOI.EqualTo(0.0),
    )
    MOI.add_constraint(
        model,
        DOI.ExplicitDifferentialFunction(
            θ,
            NDF(:*, [1.0, ω], t),            
        ),
        MOI.EqualTo(0.0),
    )
    MOI.add_constraint(
        model,
        DOI.ExplicitDifferentialFunction(
            ω,
            NDF(:/, [num_ω, den_ω], t),            
        ),
        MOI.EqualTo(0.0),
    )

    ## Objective Function
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    obj_fun = DOI.MultiPhaseIntegral([NDF(:^, [u, 2], t)])
    MOI.set(model, MOI.ObjectiveFunction{typeof(obj_fun)}(), obj_fun)

    return nothing
end