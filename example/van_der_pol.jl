function van_der_pol(
    model::Interesso.Optimizer;
    starts::AbstractDict{String,<:DOI.AbstractDynamicSolution}=Dict{String,DOI.AbstractDynamicSolution}()
) 

    @assert MOI.is_empty(model)

    ## Time as a phase
    t = DOI.add_phase(model)
    MOI.add_constraint(model, DOI.Initial(t), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Final(t), MOI.EqualTo(4.0))

    ## Input Dynamic Variable
    @variable(model, u, t)
    
    ## State Dynamic Variables
    @variable(model, x, t)
    @variable(model, v, t)

    ## Inequality constraint
    MOI.add_constraint(model, u, MOI.Interval(-1.0, 1.0))

    ## Boundary Conditions
    MOI.add_constraint(model, DOI.Initial(x), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Initial(v), MOI.EqualTo(1.0))

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
            NDF(:-, [NDF(:+, [NDF(:*, [NDF(:-, [1.0, NDF(:^, [x, 2.0], t)], t), v], t), u], t), x], t),           
        ),
        MOI.EqualTo(0.0),
    )

    ## Objective Function
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    obj_fun = DOI.MultiPhaseIntegral(
        [
            NDF(:+, [
                NDF(:*, [0.5, NDF(:^, [x, 2.0], t)], t),
                NDF(:*, [0.5, NDF(:^, [v, 2.0], t)], t),
            ], t)
        ]
    )
    MOI.set(model, MOI.ObjectiveFunction{typeof(obj_fun)}(), obj_fun)

    Interesso.warmstart!(model, starts)

    MOI.optimize!(model)

    return nothing
end