function bang_bang(
    model::Interesso.Optimizer;
    starts::Interesso.WSS=Interesso.WSS{DOI.AbstractDynamicSolution}()
) 

    MOI.empty!(model)

    ## Time as a phase
    t = DOI.add_phase(model)
    MOI.add_constraint(model, DOI.Initial(t), MOI.EqualTo(0.0))

    ## Input Dynamic Variable
    @control(model, u, t)
    MOI.add_constraint(model, u, MOI.Interval(-2.0, 1.0))
    
    ## State Dynamic Variables
    @variable(model, x, t)
    @variable(model, v, t)

    ## Boundary Conditions
    MOI.add_constraint(model, DOI.Initial(x), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Initial(v), MOI.EqualTo(0.0))

    MOI.add_constraint(model, DOI.Final(x), MOI.EqualTo(300.0))
    MOI.add_constraint(model, DOI.Final(v), MOI.EqualTo(0.0))

    ## Differential Equations

    MOI.add_constraint(
        model,
        DOI.ExplicitDifferentialFunction(
            v,
            NDF(:+, Any[u], t),            
        ),
        MOI.EqualTo(0.0),
    )

    MOI.add_constraint(
        model,
        DOI.ExplicitDifferentialFunction(
            x,
            NDF(:+, Any[v], t),            
        ),
        MOI.EqualTo(0.0),
    )

    ## Objective Function
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    obj_fun = DOI.MultiPhaseIntegral([NDF(:+, [1.0], t)])
    MOI.set(model, MOI.ObjectiveFunction{typeof(obj_fun)}(), obj_fun)

    Interesso.warmstart!(model, starts)

    return nothing
end