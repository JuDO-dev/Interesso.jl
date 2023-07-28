@testset "LegendreLobatto" begin
    
    ll3 = LegendreLobatto(2);
    ll4 = LegendreLobatto(3);
    ll5 = LegendreLobatto(4);

    @test ll3.n == 3;
    @test ll4.n == 4;
    @test ll5.n == 5;

    @test ll3.τ ≈ [-1.0, 0.0, 1.0];
    @test ll4.τ ≈ [-1.0, -sqrt(1/5), sqrt(1/5), 1.0];
    @test ll5.τ ≈ [-1.0, -sqrt(3/7), 0.0, sqrt(3/7), 1.0];

    @test ll3.w_q ≈ [1/3, 4/3, 1/3];
    @test ll4.w_q ≈ [1/6, 5/6, 5/6, 1/6];
    @test ll5.w_q ≈ [0.1, 49/90, 32/45, 49/90, 0.1];

    @test ll3.w_b ≈ [0.5, -1.0, 0.5];
end