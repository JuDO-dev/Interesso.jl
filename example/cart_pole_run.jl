using Plots
using SLOW

include(joinpath(@__DIR__, "..", "example", "cart_pole.jl"))


optimizer = SLOW.Optimizer()
MOI.set(optimizer, MOI.RawOptimizerAttribute("λ0"), nothing)
MOI.set(optimizer, MOI.RawOptimizerAttribute("ρ0"), 10.0)
MOI.set(optimizer, MOI.RawOptimizerAttribute("h_norm"), 1)
MOI.set(optimizer, MOI.RawOptimizerAttribute("max_iter"), 1000)
MOI.set(optimizer, MOI.RawOptimizerAttribute("max_time"), 60.0)

model = Interesso.Optimizer(
    # inner=optimizer,
    default_intervals=FlexibleIntervals(20, 0.0),
    # default_intervals=FixedIntervals(10),
    default_points=LGRPoints(3; order_control=2),
    default_method=IntResidual(5),
    default_bounds=SampledBounds(10)
)
u_sol, r_sol, v_sol, model = cart_pole(model)

__eval = eval_funcs(model.inner, model.dif_res_funcs)
println("maximum residual gradient")
println(maximum(__eval))
println("average residual gradient")
println(sum(__eval) / length(__eval))

_eval = eval_funcs(model.inner, model.res_funcs)
println("maximum residual")
println(maximum(_eval))
println("average residual")
println(sum(_eval) / length(_eval))

plot(tau -> r_sol(tau), r_sol.initial, r_sol.final)