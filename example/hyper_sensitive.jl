function hyper_sensitive(model::Interesso.Optimizer)

    @assert MOI.is_empty(model)

    ## Time as a phase
    t = DOI.add_phase(model)
    MOI.add_constraint(model, DOI.Initial(t), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Final(t), MOI.EqualTo(10000.0))

    ## Input Dynamic Variable
    u = DOI.add_dynamic_variable(model, t)

    ## State Dynamic Variables
    x = DOI.add_dynamic_variable(model, t)

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

    MOI.optimize!(model)

    # Retrieve solutions
    u_sol = MOI.get(model, DOI.DynamicVariableSolution(), u)
    x_sol = MOI.get(model, DOI.DynamicVariableSolution(), x)

    return u_sol, x_sol
end