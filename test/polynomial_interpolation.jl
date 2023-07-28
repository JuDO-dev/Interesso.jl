@testset "LegendreLobatto Interpolation" begin
    
    p = PolynomialInterpolation(2, LegendreLobatto);

    @test p.n == 3;
    @test p.τ ≈ [-1.0, 0.0, 1.0];
    @test p.w_b ≈ [0.5, -1.0, 0.5];

end

@testset "interpolation_vector" begin
    
    p = PolynomialInterpolation(2, LegendreLobatto);
    τ = [-1.0, -0.5, 0.0, 0.5, 1.0];
    v = [I.interpolation_vector(p, τ_j) for τ_j in τ];
    y = [1.0, 2.0, 5.0];

    @test I.Progradio.dot(v[1], y) == 1.0;
    @test I.Progradio.dot(v[2], y) ≈ 1.25;
    @test I.Progradio.dot(v[3], y) == 2.0;
    @test I.Progradio.dot(v[4], y) ≈ 3.25;
    @test I.Progradio.dot(v[5], y) == 5.0;
end

@testset "differentiation_matrix" begin
    
    p = PolynomialInterpolation(2, LegendreLobatto);
    Dᵀ = I.transposed_differentiation_matrix(p);
    τ = [-1.0, -0.5, 0.0, 0.5, 1.0];
    v = [Dᵀ * I.interpolation_vector(p, τ_j) for τ_j in τ];    
    y = [2.0, 1.0, 2.0];
    
    @test I.Progradio.dot(v[1], y) ≈ -2.0;
    @test I.Progradio.dot(v[2], y) ≈ -1.0;
    @test I.Progradio.dot(v[3], y) ≈ 0.0;
    @test I.Progradio.dot(v[4], y) ≈ 1.0;
    @test I.Progradio.dot(v[5], y) ≈ 2.0;
end