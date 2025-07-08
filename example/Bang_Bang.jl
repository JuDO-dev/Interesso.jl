import MathOptInterface as MOI
import DynOptInterface as DOI
using Interesso
using Plots
using SLOW

# Problem Constants
const NDF = DOI.NonlinearDynamicFunction

# Warm-starts
struct LinearInterpolant <: DOI.AbstractDynamicSolution
    y_a::Float64
    y_b::Float64
end
(li::LinearInterpolant)(t::Real) = li.y_a + (t - t_0) * (li.y_b - li.y_a) / (t_f - t_0)

# Problem Solver

function bang_bang(model::Interesso.Optimizer)

    @assert MOI.is_empty(model)

    ## Time as a phase
    t = DOI.add_phase(model)
    MOI.add_constraint(model, DOI.Initial(t), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Final(t), MOI.EqualTo(tf))
    # MOI.add_constraint(model, DOI.Final(t), MOI.EqualTo(2.0))

    ## Input Dynamic Variable
    u = DOI.add_dynamic_variable(model, t)
    MOI.add_constraint(model, u, MOI.Interval(-2.0, 1.0))

    ## State Dynamic Variables
    x = DOI.add_dynamic_variable(model, t)
    v = DOI.add_dynamic_variable(model, t)

    ## Boundary Conditions
    MOI.add_constraint(model, DOI.Initial(x), MOI.EqualTo(0.0))
    MOI.add_constraint(model, DOI.Initial(v), MOI.EqualTo(0.0))

    MOI.add_constraint(model, DOI.Final(x), MOI.EqualTo(300.0))
    MOI.add_constraint(model, DOI.Final(v), MOI.EqualTo(0.0))

    # # Starts
    # MOI.set(model, DOI.DynamicVariableStart(), r, LinearInterpolant(0.0, 1.0))
    # MOI.set(model, DOI.DynamicVariableStart(), θ, LinearInterpolant(0.0, 1.0 * pi))

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
    obj_fun = DOI.Bolza(
        DOI.NonlinearBoundaryFunction(:+, [0.0]),
        DOI.MultiPhaseIntegral([NDF(:+, [1.0], t)]))#
    MOI.set(model, MOI.ObjectiveFunction{typeof(obj_fun)}(), obj_fun)

    MOI.optimize!(model)

    ## Retrieve solutions
    u_sol = MOI.get(model, DOI.DynamicVariableSolution(), u)
    x_sol = MOI.get(model, DOI.DynamicVariableSolution(), x)
    v_sol = MOI.get(model, DOI.DynamicVariableSolution(), v)

    return u_sol, x_sol, v_sol
end

model = Interesso.Optimizer(
    # inner = SLOW.Optimizer(),
    default_intervals=FlexibleIntervals(4, 0.5),
    default_points=LGRPoints(4),
    default_method=PenaltyIR(10),
)
u_sol, x_sol, v_sol= bang_bang(model)

# plot(tau -> u_sol(tau), xlims=(t_sol[1], t_sol[end]))
