@testset "RigidIntervals" begin
    
    t_0 = -0.2;
    t_f = 1.4;
    n_i = 8;
    intervals = RigidIntervals(-0.2, 1.4, 8);

    @test intervals.n == n_i;
    @test intervals.t ≈ [-0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4];
end