function aly_chan(model::Interesso.Optimizer)

    @assert MOI.is_empty(model)

    ## Time as a phase
    t = DOI.add_phase(model)
    MOI.add_constraint(model, DOI.Initial(t), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Final(t), MOI.EqualTo(π/2))

    ## Input Dynamic Variable
    u = DOI.add_dynamic_variable(model, t)
    MOI.add_constraint(model, u, MOI.Interval(-0.1, 0.1))

    ## State Dynamic Variables
    x = DOI.add_dynamic_variable(model, t)
    v = DOI.add_dynamic_variable(model, t)
    cost = DOI.add_dynamic_variable(model, t)

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

    MOI.optimize!(model)

    # Retrieve solutions
    u_sol = MOI.get(model, DOI.DynamicVariableSolution(), u)
    x_sol = MOI.get(model, DOI.DynamicVariableSolution(), x)
    v_sol = MOI.get(model, DOI.DynamicVariableSolution(), v)

    return u_sol, x_sol, v_sol
end