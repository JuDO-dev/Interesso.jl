import MathOptInterface as MOI
import DynOptInterface as DOI
using Interesso
using Plots
using SLOW

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
MOI.set(optimizer, MOI.RawOptimizerAttribute("ρ0"), 10)
MOI.set(optimizer, MOI.RawOptimizerAttribute("solver"), "Clarabel")
# MOI.set(optimizer, MOI.RawOptimizerAttribute("logging"), 0)
# MOI.set(optimizer, MOI.RawOptimizerAttribute("verbose"), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("max_iter"),500)

model = Interesso.Optimizer(
    # inner=optimizer,
    # default_intervals=FlexibleIntervals(10, 0.1),
    default_intervals=FixedIntervals(50),
    default_points=LGLPoints(3),
    default_method=Collocation(),
    # default_method=DAIR(5),
    # default_method=QPM(5; pen_param=100),
    # default_bounds=SampledBounds(9),
)

lqr(model)

assess_solution(model; q=20)

sols = get_solutions(model)
sol = sols["u"]
plot(tau -> sol(tau), sol.initial, sol.final)