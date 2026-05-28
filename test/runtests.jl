using Interesso
import MathOptInterface as MOI
using Test
using SLOW

include(joinpath(@__DIR__, "..", "example", "cart_pole.jl"))

@testset "Interesso.jl" begin

    optimizer = SLOW.Optimizer()
    MOI.set(optimizer, MOI.RawOptimizerAttribute("dual"), true)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("ρ0"), 10.0)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("h_norm"), 1)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("max_iter"), 1000)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("max_time"), 60.0)

    model = Interesso.Optimizer(
        inner=optimizer,
        default_intervals=FlexibleIntervals(20, 0.0),
        default_points=LGRPoints(3, 2),
        default_method=DAIR(5),
        default_bounds=SampledBounds(10)
    )
    ~, ~, ~, model = cart_pole(model)

    status = MOI.get(model.inner, MOI.TerminationStatus())
    @test status in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)

    @test (assess_solution(model)[2] < 1e-4)

end
