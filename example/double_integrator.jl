import MathOptInterface as MOI
import DynOptInterface as DOI
using Interesso

# Problem Constants
const NDF = DOI.NonlinearDynamicFunction

# Warm-starts
struct LinearInterpolant <: DOI.AbstractDynamicSolution
    y_a::Float64
    y_b::Float64
end
(li::LinearInterpolant)(t::Real) = li.y_a + (t - t_0) * (li.y_b - li.y_a) / (t_f - t_0)

# Problem Solver

function double_integrator(model::Interesso.Optimizer)

    @assert MOI.is_empty(model)

    ## Time as a phase
    t = DOI.add_phase(model)
    MOI.add_constraint(model, DOI.Initial(t), MOI.EqualTo(10.0))
    MOI.add_constraint(model, DOI.Final(t), MOI.EqualTo(20.0))

    ## Input Dynamic Variable
    u = DOI.add_dynamic_variable(model, t)

    ## State Dynamic Variables
    x = DOI.add_dynamic_variable(model, t)
    v = DOI.add_dynamic_variable(model, t)
    τ = DOI.add_dynamic_variable(model, t)

    ## Inequality constraint
    MOI.add_constraint(model, u, MOI.Interval(-10.0, 10.0))
    MOI.add_constraint(model, x, MOI.Interval(-6.0, 6.0))
    MOI.add_constraint(model, v, MOI.Interval(-10.0, 10.0))

    ## Boundary Conditions
    MOI.add_constraint(model, DOI.Initial(x), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Initial(v), MOI.EqualTo(5.0))
    MOI.add_constraint(model, DOI.Initial(τ), MOI.EqualTo(0.0))

    ## Differential Equations
    sinτ = NDF(:sin, [τ], t)
    cosτ = NDF(:cos, [τ], t)

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

    MOI.add_constraint(
        model,
        DOI.ExplicitDifferentialFunction(
            τ,
            NDF(:+, [1.0], t),            
        ),
        MOI.EqualTo(0.0),
    )
 

    ## Objective Function
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    obj_fun = DOI.Bolza(
        DOI.NonlinearBoundaryFunction(:+, [0.0]),
        DOI.MultiPhaseIntegral(
            [
                NDF(:+, [
                    NDF(:^, [NDF(:-, [x, NDF(:*, [5.0, sinτ], t)], t), 2.0], t),
                    NDF(:^, [NDF(:-, [v, NDF(:*, [5.0, cosτ], t)], t), 2.0], t),
                    NDF(:*, [0.0001, NDF(:^, [u, 2.0], t)], t)
                ], t)
            ]
        ))
    MOI.set(model, MOI.ObjectiveFunction{typeof(obj_fun)}(), obj_fun)

    MOI.optimize!(model)

    ## Retrieve solutions
    u_sol = MOI.get(model, DOI.DynamicVariableSolution(), u)
    x_sol = MOI.get(model, DOI.DynamicVariableSolution(), x)
    v_sol = MOI.get(model, DOI.DynamicVariableSolution(), v)

    return u_sol, x_sol, v_sol
end