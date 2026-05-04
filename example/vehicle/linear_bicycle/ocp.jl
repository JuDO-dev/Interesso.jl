import MathOptInterface as MOI
import DynOptInterface as DOI
using Interesso
using Plots
using SLOW
using Ipopt
using JLD2


include(joinpath(@__DIR__, "linear_bicycle.jl"))
include(joinpath(@__DIR__, "linear_bicycle_plot.jl"))

const NDF = DOI.NonlinearDynamicFunction

# Warm-starts
struct LinearInterpolant <: DOI.AbstractDynamicSolution
    y_a::Float64
    y_b::Float64
    t_0::Float64
    t_f::Float64
end
(li::LinearInterpolant)(t::Real) = li.y_a + (t - li.t_0) * (li.y_b - li.y_a) / (li.t_f - li.t_0)


# optimizer = SLOW.Optimizer()
# MOI.set(optimizer,
#     # "globalization" => "filter",
#     "dual"     => false,
#     "ρ0"       => 100,
#     "h_norm"   => 1,
#     "γ"        => 1.0,
#     "solver"   => "FBstab",
#     "max_iter" => 10,
#     "max_time" => Inf,
#     "verbose"  => false,
#     "scaling"  => "gradient",
#     "logging"  => 0,
#     "visualize" => false,
# )

optimizer = Ipopt.Optimizer()
MOI.set(optimizer, "max_iter" => 10000)

model = Interesso.Optimizer(
    inner=optimizer,
    default_intervals=FixedIntervals(50),
    default_points=LGRPoints(4),
    # default_method=Collocation(),
    # default_method=DAIR(5),
    default_method=DAIROpti(5;tolerance=1e-6),
    # default_method=QPM(5;penalty=0.0001),
    # default_method=SAIR(5),
    # default_bounds=SampledBounds(9)
)

# primal = nothing
@load joinpath(@__DIR__, "race_car_primals.jld2") primal

trackfile = joinpath(@__DIR__, "../tracks", "txt/catalunya_2022_S1.txt")

linear_bicycle(model, trackfile)

MOI.optimize!(model; primal)

# primal = get_primal(model)
# @save joinpath(@__DIR__, "race_car_primals.jld2") primal

# (_, R, _, _) = assess_solution(model);
# res = Float64[]
# push!(res, R)
# map = 1:50
# for i in map
#     refine!(model)
#     (_, R, _, _) = assess_solution(model);
#     push!(res, R)
# end
# map = append!([0], map)
# plt = Plots.plot(map, res)
# savefig(plt, "race_car.svg")

p = Interesso.plot(model)
display(p)

# p1 = plot_trajectory(model, trackfile)
# p2 = plot_solution(model)
# p3 = plot_curvature(trackfile)
# display(p1)
# display(p2)
# display(p3)
