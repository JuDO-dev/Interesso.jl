@testset "Barycentric weights" begin
    
    τ = [-1.0, 0.0, 1.0];
    w_b = I.barycentric_weights(τ);

    @test w_b ≈ [0.5, -1.0, 0.5];
end

@testset "Barycentric interpolation" begin

    τ = [-1.0, 0.0, 1.0];
    w_b = I.barycentric_weights(τ);

    v1 = I.barycentric_interpolation(τ, w_b, 0.5);
    @test v1[1] * 1.0 + v1[2] * 2.0 + v1[3] * 3.0 ≈ 2.5;
    
    v2 = I.barycentric_interpolation(τ, w_b, 1.0);
    @test v2[1] * 1.0 + v2[2] * 2.0 + v2[3] * 3.0 ≈ 3.0;
end

@testset "Barycentric differentiation" begin
    
    τ = [-1.0, 0.0, 1.0];
    w_b = I.barycentric_weights(τ);

    D = I.barycentric_differentiation(τ, w_b);

end