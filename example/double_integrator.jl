function double_integrator(
    model::Interesso.Optimizer;
    starts::Interesso.WSS=Interesso.WSS{DOI.AbstractDynamicSolution}()
) 

    MOI.empty!(model)

    ## Time as a phase
    t = DOI.add_phase(model)
    MOI.add_constraint(model, DOI.Initial(t), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Final(t), MOI.EqualTo(10.0))

    ## Input Dynamic Variable
    @variable(model, u, t)
    
    ## State Dynamic Variables
    @variable(model, x, t)
    @variable(model, v, t)

    ## Inequality constraint
    MOI.add_constraint(model, u, MOI.Interval(-10.0, 10.0))
    MOI.add_constraint(model, x, MOI.Interval(-6.0, 6.0))
    MOI.add_constraint(model, v, MOI.Interval(-10.0, 10.0))

    ## Boundary Conditions
    MOI.add_constraint(model, DOI.Initial(x), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Initial(v), MOI.EqualTo(5.0))

    ## Differential Equations
    sin_t = NDF(:sin, [t], t)
    cos_t = NDF(:cos, [t], t)

    MOI.add_constraint(
        model,
        DOI.ExplicitDifferentialFunction(
            x,
            NDF(:+, Any[v], t),            
        ),
        MOI.EqualTo(0.0),
    )

    MOI.add_constraint(
        model,
        DOI.ExplicitDifferentialFunction(
            v,
            NDF(:+, Any[u], t),            
        ),
        MOI.EqualTo(0.0),
    )
 

    ## Objective Function
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    obj_fun = DOI.MultiPhaseIntegral(
        [
            NDF(:+, [
                NDF(:^, [NDF(:-, [x, NDF(:*, [5.0, sin_t], t)], t), 2.0], t),
                NDF(:^, [NDF(:-, [v, NDF(:*, [5.0, cos_t], t)], t), 2.0], t),
                NDF(:*, [0.0001, NDF(:^, [u, 2.0], t)], t)
            ], t)
        ]
    )
    MOI.set(model, MOI.ObjectiveFunction{typeof(obj_fun)}(), obj_fun)

    Interesso.warmstart!(model, starts)

    MOI.optimize!(model)

    return nothing
end