@testset "Context" begin
    
    n_i = 4; n_x = 6; n_u = 2; n_p = 3; n_r = 5;
    n_τ_x = 4; n_τ_u = 3; n_τ_q = 5;

    ctx = I.Context(n_i, n_x, n_u, n_p, n_r, n_τ_x, n_τ_u, n_τ_q);

    @test ctx.x_offset == 24;
    @test ctx.u_offset == 18;
    @test ctx.p_offset == 102;

    @test ctx.z_offsets == [0, 24, 48, 72];
    @test ctx.x_offsets == [5, 7, 9, 11, 13, 15];
    @test ctx.u_offsets == [18, 21];
end

@testset "PrimalDiscretization" begin
    
    n_i = 4; n_x = 6; n_u = 2; n_p = 3; n_r = 5;
    n_τ_x = 4; n_τ_u = 3; n_τ_q = 5;
    ctx = I.Context(n_i, n_x, n_u, n_p, n_r, n_τ_x, n_τ_u, n_τ_q);

    pd = I.PrimalDiscretization{Float64}(undef, ctx);

    @test length(pd.z) == 105;
    @test length(pd.x) == n_i;
    @test length(pd.x[1]) == n_τ_q;
    @test length(pd.x[1][1]) == n_x;
end