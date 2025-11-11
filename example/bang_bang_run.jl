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

include(joinpath(@__DIR__, "..", "example", "bang_bang.jl"))


# Problem Solver
optimizer = SLOW.Optimizer()
MOI.set(optimizer, MOI.RawOptimizerAttribute("dual"), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("ρ0"), 10.0)
MOI.set(optimizer, MOI.RawOptimizerAttribute("h_norm"), 1)
MOI.set(optimizer, MOI.RawOptimizerAttribute("logging"), 2)
MOI.set(optimizer, MOI.RawOptimizerAttribute("max_iter"), 10000)

model = Interesso.Optimizer(
    inner=optimizer,
    # inner = Uno.Optimizer(preset="filtersqp"),
    # default_intervals=FlexibleIntervals(30, 0.1),
    default_intervals=FixedIntervals(30),
    default_points=LGRPoints(3),
    # default_method=Collocation(),
    default_method=QPM(5),
    # default_bounds=SampledBounds(5),
)
u_sol, x_sol, v_sol= bang_bang(model)

# __eval = eval_funcs(model.inner, model.dif_res_funcs)
# abs__eval = abs.(__eval)
# println("maximum 1-norm residual gradient")
# println(maximum(abs__eval))
# println("average 1-norm residual gradient")
# println(sum(abs__eval) / length(__eval))

# _eval = eval_funcs(model.inner, model.res_funcs)
# abs_eval = abs.(_eval)
# println("maximum 1-norm residual")
# println(maximum(abs_eval))
# println("average 1-norm residual")
# println(sum(abs_eval) / length(_eval))

plot(tau -> u_sol(tau), x_sol.initial, x_sol.final)