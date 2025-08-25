using Interesso
import MathOptInterface as MOI
using Test
using SLOW

include(joinpath(@__DIR__, "..", "example", "cart_pole.jl"))

@testset "Interesso.jl" begin

    optimizer = SLOW.Optimizer()
    MOI.set(optimizer, MOI.RawOptimizerAttribute("λ0"), nothing)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("ρ0"), 10.0)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("h_norm"), 1)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("max_iter"), 1000)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("max_time"), 60.0)

    model = Interesso.Optimizer(
        inner=optimizer,
        default_intervals=FlexibleIntervals(20, 0.0),
        default_points=LGRPoints(3; order_control=2),
        default_method=IntResidual(5),
        default_bounds=SampledBounds(10)
    )
    ~, ~, ~, model = cart_pole(model)

    status = MOI.get(model.inner, MOI.TerminationStatus())
    @test status in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)

    _eval = eval_funcs(model.inner, model.res_funcs)
    @test all(_eval .≤ 1e-1)
end
