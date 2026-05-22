function goddard_rocket(
    model::Interesso.Optimizer;
    starts::Interesso.WSS=Interesso.WSS{DOI.AbstractDynamicSolution}()
)

    MOI.empty!(model)

    ## Problem constants
    D0 = 5.4915e-05
    H = 23800.0
    T_min = 0.0
    T_max = 193.0
    g = 32.174
    c = 1580.9425

    h_0 = 0.0
    v_0 = 0.0
    m_0 = 3.0

    h_min = 0.0
    h_max = 23800.0
    v_min = -10.0
    v_max = 900.0
    m_min = 1.0
    m_max = 3.0

    t_0 = 0.0
    t_f = 45.0

    ## Time as a phase
    t = DOI.add_phase(model)
    MOI.add_constraint(model, DOI.Initial(t), MOI.EqualTo(t_0))
    MOI.add_constraint(model, DOI.Final(t), MOI.EqualTo(t_f))

    ## Input Dynamic Variable
    @variable(model, T, t)
    MOI.add_constraint(model, T, MOI.Interval(T_min, T_max))

    ## State Dynamic Variables
    @variable(model, h, t)
    MOI.add_constraint(model, h, MOI.Interval(h_min, h_max))

    @variable(model, v, t)
    MOI.add_constraint(model, v, MOI.Interval(v_min, v_max))

    @variable(model, m, t)
    MOI.add_constraint(model, m, MOI.Interval(m_min, m_max))

    ## Boundary Conditions
    MOI.add_constraint(model, DOI.Initial(h), MOI.EqualTo(h_0))
    MOI.add_constraint(model, DOI.Initial(v), MOI.EqualTo(v_0))
    MOI.add_constraint(model, DOI.Initial(m), MOI.EqualTo(m_0))

    # MOI.add_constraint(model, DOI.Final(m), MOI.EqualTo(m_min))

    ## Differential Equations
    drag = NDF(
        :*,
        [
            D0,
            NDF(:^, [v, 2.0], t),
            NDF(:exp, [NDF(:/, [NDF(:*, [-1.0, h], t), H], t)], t),
        ],
        t,
    )

    MOI.add_constraint(
        model,
        DOI.ExplicitDifferentialFunction(
            h,
            NDF(:+, [v], t),
        ),
        MOI.EqualTo(0.0),
    )

    MOI.add_constraint(
        model,
        DOI.ExplicitDifferentialFunction(
            v,
            NDF(
                :+,
                [
                    NDF(:/, [NDF(:-, [T, drag], t), m], t),
                    -g,
                ],
                t,
            ),
        ),
        MOI.EqualTo(0.0),
    )

    MOI.add_constraint(
        model,
        DOI.ExplicitDifferentialFunction(
            m,
            NDF(:/, [NDF(:*, [-1.0, T], t), c], t),
        ),
        MOI.EqualTo(0.0),
    )

    ## Objective Function
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    obj_fun = DOI.NonlinearBoundaryFunction(:*, [DOI.Final(h), 1e-4])
    MOI.set(model, MOI.ObjectiveFunction{typeof(obj_fun)}(), obj_fun)

    ## Warm-starts
    if starts == Interesso.WSS{DOI.AbstractDynamicSolution}()
        MOI.set(model, DOI.DynamicVariableStart(), T, LinearInterpolant(T_max, T_min, t_0, t_f))
        MOI.set(model, DOI.DynamicVariableStart(), h, LinearInterpolant(h_0, 18000.0, t_0, t_f))
        MOI.set(model, DOI.DynamicVariableStart(), v, LinearInterpolant(v_0, 0.0, t_0, t_f))
        MOI.set(model, DOI.DynamicVariableStart(), m, LinearInterpolant(m_0, m_min, t_0, t_f))
    else
        Interesso.warmstart!(model, starts)
    end

    return nothing
end
