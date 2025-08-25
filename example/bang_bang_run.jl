using Plots
using SLOW

include(joinpath(@__DIR__, "..", "example", "bang_bang.jl"))


model = Interesso.Optimizer(
    inner = SLOW.Optimizer(),
    default_intervals=FlexibleIntervals(20, 0.1),
    default_points=LGRPoints(3; order_control=2),
    default_method=IntResidual(5),
    default_bounds=SampledBounds(10)
)
u_sol, x_sol, v_sol= bang_bang(model)

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

plot(tau -> u_sol(tau), u_sol.initial, u_sol.final)