import MathOptInterface as MOI
import DynOptInterface as DOI
using Interesso
using Plots
using SLOW
using Ipopt
using JLD2


include(joinpath(@__DIR__, "aly_chan.jl"))
include(joinpath(@__DIR__, "bang_bang.jl"))
include(joinpath(@__DIR__, "bang_bang_spatial.jl"))
include(joinpath(@__DIR__, "cart_pole.jl"))
include(joinpath(@__DIR__, "cart_pole_implicit.jl"))
include(joinpath(@__DIR__, "double_integrator.jl"))
include(joinpath(@__DIR__, "fuller.jl"))
include(joinpath(@__DIR__, "hyper_sensitive.jl"))
include(joinpath(@__DIR__, "lqr.jl"))
include(joinpath(@__DIR__, "orbit_raising.jl"))
include(joinpath(@__DIR__, "space_shuttle_reentry.jl"))
include(joinpath(@__DIR__, "two_link_robot_arm.jl"))
include(joinpath(@__DIR__, "van_der_pol.jl"))


const NDF = DOI.NonlinearDynamicFunction

# Warm-starts
struct LinearInterpolant <: DOI.AbstractDynamicSolution
    y_a::Float64
    y_b::Float64
    t_0::Float64
    t_f::Float64
end
(li::LinearInterpolant)(t::Real) = li.y_a + (t - li.t_0) * (li.y_b - li.y_a) / (li.t_f - li.t_0)


optimizer = SLOW.Optimizer()
MOI.set(optimizer,
    "dual"     => true,
    "ρ0"       => 10,
    "h_norm"   => 2,
    "γ"        => 1.0,
    "solver"   => "Clarabel",
    "max_iter" => 1000,
    "max_time" => Inf,
    "verbose"  => false,
    "scaling"  => "none",
    "logging"  => 2,
)

optimizer = Ipopt.Optimizer()
MOI.set(optimizer, "max_iter" => 2000)

primal = nothing

model = Interesso.Optimizer(
    inner=optimizer,
    # default_intervals=FlexibleIntervals(10, 0.1),
    default_intervals=FixedIntervals(50),
    # default_points=LGRPoints(3),
    # default_method=Collocation(),
    # default_method=DAIROpti(5),
    default_method=Galerkin(5)
    # default_method=QPM(5;penalty=0.0001),
    # default_method=SAIR(5),
    # default_bounds=SampledBounds(9)
)

van_der_pol(model)

MOI.optimize!(model; primal)

# (_, RR, _, _) = assess_solution(model);

# res = Float64[]

# push!(res, RR)

# map = 1:1

# for i in map
#     refine!(model)
#     (_, RR, _, _) = assess_solution(model);
#     push!(res, RR)
# end

# map = append!([0], map)

# plt = Plots.plot(map, res)

# savefig(plt, "cart_pole.svg")

Interesso.plot(model)
