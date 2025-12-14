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
(li::LinearInterpolant)(t::Real) = li.y_a + (t - t_0) * (li.y_b - li.y_a) / (t_f - t_0)

include(joinpath(@__DIR__, "..", "example", "cart_pole.jl"))
include(joinpath(@__DIR__, "..", "example", "cart_pole_implicit.jl"))

optimizer = SLOW.Optimizer()
MOI.set(optimizer, MOI.RawOptimizerAttribute("dual"), false)
MOI.set(optimizer, MOI.RawOptimizerAttribute("ρ0"), 10.0)
MOI.set(optimizer, MOI.RawOptimizerAttribute("h_norm"), 1)
MOI.set(optimizer, MOI.RawOptimizerAttribute("max_iter"), 300)
# MOI.set(optimizer, MOI.RawOptimizerAttribute("max_time"), 60.0)
MOI.set(optimizer, MOI.RawOptimizerAttribute("logging"), 2)


model = Interesso.Optimizer(
    inner=optimizer,
    # inner = Uno.Optimizer(preset="filtersqp"),
    # default_intervals=FlexibleIntervals(30, 0.1),
    default_intervals=FixedIntervals(50),
    default_points=LGRPoints(3),
    # default_method=Collocation(),
    default_method=SAIR(7),
    # default_method=QPM(5; pen_param=100),
    default_bounds=SampledBounds(9),
)

u_sol, r_sol, v_sol, model = cart_pole(model)

assess_solution(model;q=20)

# ws = Interesso.get_solutions(model)
# MOI.empty!(model)

# cart_pole(model; starts = ws)

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

# ws = perturb_solutions(Interesso.get_solutions(model), 0.0)

# MOI.empty!(model)

# MOI.set(model.inner, MOI.RawOptimizerAttribute("max_iter"), 0)
# u_sol, r_sol, v_sol, model = cart_pole(model; starts = ws)

# plot(tau -> r_sol(tau), r_sol.initial, r_sol.final)


# open("ws.txt", "w") do io
#     println(io, ws)
# end

# open("ws1.txt", "w") do io
#     println(io, ws1)
# end