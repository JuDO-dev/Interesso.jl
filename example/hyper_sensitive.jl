function hyper_sensitive(
    model::Interesso.Optimizer;
    starts::Interesso.WSS=Interesso.WSS{DOI.AbstractDynamicSolution}()
) 

    MOI.empty!(model)

    ## Time as a phase
    t = DOI.add_phase(model)
    MOI.add_constraint(model, DOI.Initial(t), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Final(t), MOI.EqualTo(10000.0))

    ## Input Dynamic Variable
    @variable(model, u, t)
    
    ## State Dynamic Variables
    @variable(model, x, t)

    ## Boundary Conditions
    MOI.add_constraint(model, DOI.Initial(x), MOI.EqualTo(1.0))
    MOI.add_constraint(model, DOI.Final(x), MOI.EqualTo(1.5))

    # # Starts

    ## Differential Equations

    MOI.add_constraint(
        model,
        DOI.ExplicitDifferentialFunction(
            x,
            NDF(:-, [u, NDF(:^, [x, 3], t)], t),            
        ),
        MOI.EqualTo(0.0),
    )

    ## Objective Function
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    obj_fun = DOI.MultiPhaseIntegral([
        NDF(:+, [
            NDF(:^, [x, 2.0], t),
            NDF(:^, [u, 2.0], t),
        ], t)
    ])
    MOI.set(model, MOI.ObjectiveFunction{typeof(obj_fun)}(), obj_fun)

    Interesso.warmstart!(model, starts)

    return nothing
end