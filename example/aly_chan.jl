function aly_chan(
    model::Interesso.Optimizer;
    starts::Interesso.WSS=Interesso.WSS{DOI.AbstractDynamicSolution}()
)

    MOI.empty!(model)

    ## Time as a phase
    t = DOI.add_phase(model)
    MOI.add_constraint(model, DOI.Initial(t), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Final(t), MOI.EqualTo(π/2))

    ## Input Dynamic Variable
    @control(model, u, t)
    MOI.add_constraint(model, u, MOI.Interval(-1.0, 1.0))

    ## State Dynamic Variables
    @variable(model, x, t)
    @variable(model, v, t)
    @variable(model, cost, t)

    ## Boundary Conditions
    MOI.add_constraint(model, DOI.Initial(x), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Initial(v), MOI.EqualTo(1.0))
    MOI.add_constraint(model, DOI.Initial(cost), MOI.EqualTo(0.0))

    # # Starts

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

    MOI.add_constraint(
        model,
        DOI.ExplicitDifferentialFunction(
            cost,
            NDF(:-, [
                    NDF(:*, [0.5, NDF(:^, [v, 2.0], t)], t),
                    NDF(:*, [0.5, NDF(:^, [x, 2.0], t)], t),
                ], t)           
        ),
        MOI.EqualTo(0.0),
    )

    ## Objective Function
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    obj_fun = DOI.NonlinearBoundaryFunction(:+, [DOI.Final(cost)])
    MOI.set(model, MOI.ObjectiveFunction{typeof(obj_fun)}(), obj_fun)

    Interesso.warmstart!(model, starts)

    return nothing
end