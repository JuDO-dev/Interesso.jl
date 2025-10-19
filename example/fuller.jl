function fuller(model::Interesso.Optimizer)

    @assert MOI.is_empty(model)

    ## Time as a phase
    t = DOI.add_phase(model)
    MOI.add_constraint(model, DOI.Initial(t), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Final(t), MOI.EqualTo(30.0))

    ## Input Dynamic Variable
    u = DOI.add_dynamic_variable(model, t)
    MOI.add_constraint(model, u, MOI.Interval(-0.1, 0.1))

    ## State Dynamic Variables
    x = DOI.add_dynamic_variable(model, t)
    v = DOI.add_dynamic_variable(model, t)

    ## Boundary Conditions
    MOI.add_constraint(model, DOI.Initial(x), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Initial(v), MOI.EqualTo(1.0))

    MOI.add_constraint(model, DOI.Final(x), MOI.EqualTo(0.0))
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
    obj_fun = DOI.MultiPhaseIntegral([NDF(:^, [x, 2], t)])
    MOI.set(model, MOI.ObjectiveFunction{typeof(obj_fun)}(), obj_fun)

    MOI.optimize!(model)

    # Retrieve solutions
    u_sol = MOI.get(model, DOI.DynamicVariableSolution(), u)
    x_sol = MOI.get(model, DOI.DynamicVariableSolution(), x)
    v_sol = MOI.get(model, DOI.DynamicVariableSolution(), v)

    return u_sol, x_sol, v_sol
end