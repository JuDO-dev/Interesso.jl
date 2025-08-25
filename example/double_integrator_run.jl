using Plots
using SLOW

include(joinpath(@__DIR__, "..", "example", "double_integrator.jl"))


model = Interesso.Optimizer(
    inner = SLOW.Optimizer(),
    default_intervals=FlexibleIntervals(4, 0.5),
    default_points=LGRPoints(8),
    default_method=IntResidual(10),
)
u_sol, x_sol, v_sol = double_integrator(model)

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

plot(tau -> x_sol(tau), x_sol.initial, x_sol.final)