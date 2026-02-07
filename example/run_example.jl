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

# include(joinpath(@__DIR__, "..", "example", "aly_chan.jl"))
# include(joinpath(@__DIR__, "..", "example", "bang_bang.jl"))
# include(joinpath(@__DIR__, "..", "example", "cart_pole.jl"))
# include(joinpath(@__DIR__, "..", "example", "cart_pole_implicit.jl"))
# include(joinpath(@__DIR__, "..", "example", "double_integrator.jl"))
# include(joinpath(@__DIR__, "..", "example", "fuller.jl"))
# include(joinpath(@__DIR__, "..", "example", "hyper_sensitive.jl"))
# include(joinpath(@__DIR__, "..", "example", "lqr.jl"))
# include(joinpath(@__DIR__, "..", "example", "orbit_raising.jl"))
include(joinpath(@__DIR__, "..", "example", "van_der_pol.jl"))
# include(joinpath(@__DIR__, "..", "example", "vehicle", "vehicle.jl"))
# include(joinpath(@__DIR__, "..", "example", "vehicle_2.jl"))
# include(joinpath(@__DIR__, "..", "example", "vehicle", "vehicle_simple.jl"))

optimizer = SLOW.Optimizer()
MOI.set(optimizer,
    "dual"     => true,
    "ρ0"       => 10,
    "h_norm"   => 1,
    "solver"   => "Clarabel",
    "max_iter" => 3000,
    "max_time" => Inf,
    "verbose"  => false,
    "logging"  => 2,
)



model = Interesso.Optimizer(
    inner=optimizer,
    # default_intervals=FlexibleIntervals(50, 0.1),
    default_intervals=FixedIntervals(50),
    default_points=LGLPoints(3),
    default_method=Collocation(),
    # default_method=QPM(5;pen_param = 1),
    # default_method=SAIR(5),
    # default_bounds=SampledBounds(9)
)
# cart_pole(model)
# hyper_sensitive(model)
# orbit_raising(model)
van_der_pol(model)
# lqr(model)
# bang_bang(model)

assess_solution(model; q=20)

sols = get_solutions(model)
sol = sols["u"]
display(plot(tau -> sol(tau), sol.initial, sol.final))

# # ws = perturb_solutions(Interesso.get_solutions(model), 0.3)

# model2 = Interesso.Optimizer(
#     inner=optimizer,
#     # default_intervals=FlexibleIntervals(50, 0.01),
#     default_intervals=FixedIntervals(40),
#     default_points=LGRPoints(3),
#     default_method=Collocation(),
#     # default_method=QPM(5;pen_param = 1),
#     # default_method=SAIR(5),
#     # default_bounds=SampledBounds(9)
# )
# cart_pole(model2;starts=sols)
# # lqr(model2; starts = ws);