function orbit_raising(model::Interesso.Optimizer)

    @assert MOI.is_empty(model)

    ## Time as a phase
    t = DOI.add_phase(model)
    MOI.add_constraint(model, DOI.Initial(t), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Final(t), MOI.EqualTo(3.32))

    ## Input Dynamic Variable
    u1 = DOI.add_dynamic_variable(model, t)
    u2 = DOI.add_dynamic_variable(model, t)
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
    r = DOI.add_dynamic_variable(model, t)
    MOI.add_constraint(model, r, MOI.Interval(0.0, 2.0))
    θ = DOI.add_dynamic_variable(model, t)
    MOI.add_constraint(model, θ, MOI.Interval(0.0, π))
    v_r = DOI.add_dynamic_variable(model, t)
    MOI.add_constraint(model, v_r, MOI.Interval(0.0, 2.0))
    v_θ = DOI.add_dynamic_variable(model, t)
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

    MOI.optimize!(model)

    # Retrieve solutions
    u1_sol = MOI.get(model, DOI.DynamicVariableSolution(), u1)
    u2_sol = MOI.get(model, DOI.DynamicVariableSolution(), u2)
    r_sol = MOI.get(model, DOI.DynamicVariableSolution(), r)
    θ_sol = MOI.get(model, DOI.DynamicVariableSolution(), θ)
    vr_sol = MOI.get(model, DOI.DynamicVariableSolution(), v_r)
    vθ_sol = MOI.get(model, DOI.DynamicVariableSolution(), v_θ)

    return u1_sol, u2_sol, r_sol, θ_sol, vr_sol, vθ_sol
end