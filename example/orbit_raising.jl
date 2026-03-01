function orbit_raising(
    model::Interesso.Optimizer;
    starts::Interesso.WSS=Interesso.WSS{DOI.AbstractDynamicSolution}()
) 

    MOI.empty!(model)

    ## Time as a phase
    t = DOI.add_phase(model)
    t_0 = 0.0
    t_f = 3.32
    MOI.add_constraint(model, DOI.Initial(t), MOI.EqualTo(t_0))
    MOI.add_constraint(model, DOI.Final(t), MOI.EqualTo(t_f))

    ## Input Dynamic Variable
    @variable(model, u1, t)
    @variable(model, u2, t)
    MOI.add_constraint(
        model, 
        NDF(
            :+, 
            [
                NDF(:^, [u1, 2.0], t),
                NDF(:^, [u2, 2.0], t),
                -1.0
            ],
            t
        ),
        MOI.EqualTo(0.0)
    )

    ## State Dynamic Variables
    @variable(model, r, t)
    MOI.add_constraint(model, r, MOI.Interval(0.0, 2.0))
    @variable(model, θ, t)
    MOI.add_constraint(model, θ, MOI.Interval(0.0, π))
    @variable(model, v_r, t)
    MOI.add_constraint(model, v_r, MOI.Interval(0.0, 2.0))
    @variable(model, v_θ, t)
    MOI.add_constraint(model, v_θ, MOI.Interval(0.0, 2.0))

    ## Boundary Conditions
    MOI.add_constraint(model, DOI.Initial(r), MOI.EqualTo(1.0))
    MOI.add_constraint(model, DOI.Initial(θ), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Initial(v_r), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Initial(v_θ), MOI.EqualTo(1.0))

    MOI.add_constraint(model, DOI.Final(v_r), MOI.EqualTo(0.0))
    MOI.add_constraint(
        model, 
        DOI.NonlinearBoundaryFunction(
            :-,
            [
                DOI.NonlinearBoundaryFunction(
                    :*,
                    [
                        DOI.Final(v_θ),
                        DOI.NonlinearBoundaryFunction(:sqrt, [DOI.Final(r)]),
                    ],
                ),
                1.0,
            ],
        ),
        MOI.EqualTo(0.0))

    ## Differential Equations

    MOI.add_constraint(
        model,
        DOI.ExplicitDifferentialFunction(
            r,
            NDF(:+, [v_r], t),            
        ),
        MOI.EqualTo(0.0),
    )

    MOI.add_constraint(
        model,
        DOI.ExplicitDifferentialFunction(
            θ,
            NDF(:/, [v_θ, r], t),            
        ),
        MOI.EqualTo(0.0),
    )
    
    a_t = NDF(:/, [0.1405, NDF(:-, [1.0, NDF(:*, [0.0749, t], t)], t)], t)

    MOI.add_constraint(
        model,
        DOI.ExplicitDifferentialFunction(
            v_r,
            NDF(:+,
                [
                    NDF(:*, [a_t, u1], t),
                    NDF(:/, [NDF(:^, [v_θ, 2.0], t), r], t),
                    NDF(:/, [-1.0, NDF(:^, [r, 2.0], t)], t),
                ],
                t,
            )            
        ),
        MOI.EqualTo(0.0),
    )
    
    MOI.add_constraint(
        model,
        DOI.ExplicitDifferentialFunction(
            v_θ,
            NDF(:-,
                [
                    NDF(:*, [a_t, u2], t),
                    NDF(:/, [NDF(:*, [v_θ, v_r], t), r], t)
                ], t),            
        ),
        MOI.EqualTo(0.0),
    )

    ## Objective Function
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    obj_fun = DOI.NonlinearBoundaryFunction(:+, [DOI.Final(r)])
    MOI.set(model, MOI.ObjectiveFunction{typeof(obj_fun)}(), obj_fun)

    if starts == Interesso.WSS{DOI.AbstractDynamicSolution}()
        MOI.set(model, DOI.DynamicVariableStart(), u1, LinearInterpolant(1.0, 1.0, t_0, t_f))
        MOI.set(model, DOI.DynamicVariableStart(), u2, LinearInterpolant(1.0, 1.0, t_0, t_f))
        MOI.set(model, DOI.DynamicVariableStart(), r, LinearInterpolant(1.0, 1.0, t_0, t_f))
        MOI.set(model, DOI.DynamicVariableStart(), θ, LinearInterpolant(1.0, 1.0, t_0, t_f))
        MOI.set(model, DOI.DynamicVariableStart(), v_r, LinearInterpolant(1.0, 1.0, t_0, t_f))
        MOI.set(model, DOI.DynamicVariableStart(), v_θ, LinearInterpolant(1.0, 1.0, t_0, t_f))
    else
        Interesso.warmstart!(model, starts)
    end

    return nothing
end