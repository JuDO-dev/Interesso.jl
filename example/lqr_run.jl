import MathOptInterface as MOI
import DynOptInterface as DOI
using Interesso
using Plots
using SLOW
using Uno


const NDF = DOI.NonlinearDynamicFunction

# Warm-starts
struct LinearInterpolant <: DOI.AbstractDynamicSolution
    y_a::Float64
    y_b::Float64
end
(li::LinearInterpolant)(t::Real) = li.y_a + (t) * (li.y_b - li.y_a)

include(joinpath(@__DIR__, "..", "example", "lqr.jl"))


# Problem Solver
optimizer = SLOW.Optimizer()
MOI.set(optimizer, MOI.RawOptimizerAttribute("dual"), false)
MOI.set(optimizer, MOI.RawOptimizerAttribute("ρ0"), 1000)
MOI.set(optimizer, MOI.RawOptimizerAttribute("h_norm"), 1)
MOI.set(optimizer, MOI.RawOptimizerAttribute("logging"), 2)
MOI.set(optimizer, MOI.RawOptimizerAttribute("max_iter"),3000)

model = Interesso.Optimizer(
    inner=optimizer,
    # inner = Uno.Optimizer(preset="filtersqp"),
    # default_intervals=FlexibleIntervals(30, 0.1),
    default_intervals=FixedIntervals(10),
    default_points=LGRPoints(3),
    # default_method=Collocation(),
    default_method=SAIR(5)
    # default_method=QPM(5; pen_param=100),
    # default_bounds=SampledBounds(5),
)
u_sol, x_sol, v_sol, model = lqr(model)

assess_solution(model; q=20)

plot(tau -> u_sol(tau), x_sol.initial, x_sol.final)