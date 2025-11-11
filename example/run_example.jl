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
(li::LinearInterpolant)(t::Real) = li.y_a + (t - t_0) * (li.y_b - li.y_a) / (t_f - t_0)

# include(joinpath(@__DIR__, "..", "example", "aly_chan.jl"))
include(joinpath(@__DIR__, "..", "example", "bang_bang.jl"))
include(joinpath(@__DIR__, "..", "example", "cart_pole.jl"))
# include(joinpath(@__DIR__, "..", "example", "cart_pole_implicit.jl"))
# include(joinpath(@__DIR__, "..", "example", "double_integrator.jl"))
# include(joinpath(@__DIR__, "..", "example", "fuller.jl"))
# include(joinpath(@__DIR__, "..", "example", "hyper_sensitive.jl"))
include(joinpath(@__DIR__, "..", "example", "orbit_raising.jl"))
# include(joinpath(@__DIR__, "..", "example", "van_der_pol.jl"))
# include(joinpath(@__DIR__, "..", "example", "vehicle.jl"))
# include(joinpath(@__DIR__, "..", "example", "vehicle_2.jl"))

optimizer = SLOW.Optimizer()
MOI.set(optimizer, MOI.RawOptimizerAttribute("dual"), true)
MOI.set(optimizer, MOI.RawOptimizerAttribute("ρ0"), 10.0)
MOI.set(optimizer, MOI.RawOptimizerAttribute("h_norm"), 1)
MOI.set(optimizer, MOI.RawOptimizerAttribute("max_iter"), 1000)
MOI.set(optimizer, MOI.RawOptimizerAttribute("max_time"), Inf)
MOI.set(optimizer, MOI.RawOptimizerAttribute("scaling"), "none")
MOI.set(optimizer, MOI.RawOptimizerAttribute("logging"), 2)

model = Interesso.Optimizer(
    inner=optimizer,
    # default_intervals=FlexibleIntervals(50, 0.01),
    default_intervals=FixedIntervals(50),
    # default_points=LGRPoints(3, 1),
    default_points=LGRPoints(3),
    default_method=DAIR(5),
    # default_method=Collocation(),
    # default_bounds=SampledBounds(5)
)
# u_sol, x_sol, v_sol = cart_pole_im(model)
u1_sol, u2_sol, r_sol, θ_sol, vr_sol, vθ_sol, model = orbit_raising(model)

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

# x_sol = sol.s

# plot(tau -> u1_sol(tau), u1_sol.initial, u1_sol.final)

# ws = perturb_solutions(Interesso.get_solutions(model), 0.01)

# println(ws)

# model2 = Interesso.Optimizer(
#     inner=optimizer,
#     # default_intervals=FlexibleIntervals(50, 0.01),
#     default_intervals=FixedIntervals(50),
#     default_points=LGRPoints(3),
#     default_method=DAIR(5),
# )

# orbit_raising(model2; starts = ws)