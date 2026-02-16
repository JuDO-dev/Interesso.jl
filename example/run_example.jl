import MathOptInterface as MOI
import DynOptInterface as DOI
using Interesso
using Plots
using SLOW
using Ipopt


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


# optimizer = SLOW.Optimizer()
# MOI.set(optimizer,
#     "dual"     => false,
#     "ρ0"       => 10.0,
#     "h_norm"   => 1,
#     "solver"   => "Clarabel",
#     "max_iter" => 1000,
#     "max_time" => Inf,
#     "verbose"  => false,
#     "scaling"  => "gradient",
#     "logging"  => 2,
# )

optimizer = Ipopt.Optimizer()
MOI.set(optimizer, "max_iter" => 100_000)

model = Interesso.Optimizer(
    inner=optimizer,
    # default_intervals=FlexibleIntervals(50, 0.1),
    default_intervals=FixedIntervals(20),
    default_points=LGLPoints(3),
    # default_method=Collocation(),
    default_method=DAIR(5),
    # default_method=QPM(5;pen_param = 10),
    # default_method=SAIR(5),
    # default_bounds=SampledBounds(9)
)

two_link_robot_arm(model)

assess_solution(model; q=20)

sols = get_solutions(model)
sol = sols[model.phases[1]]["u1"]
display(plot(tau -> sol(tau), sol.initial, sol.final))
