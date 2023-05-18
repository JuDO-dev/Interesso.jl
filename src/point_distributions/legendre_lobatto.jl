struct LegendreLobatto{T} <: PointDistribution{T}
    n::Int
    τ::Vector{T}
    w_q::Vector{T}
    w_b::Vector{T}

    function LegendreLobatto{T}(degree) where {T}

        n = degree + 1;

        τ_float64, w_q_float64 = FastGaussQuadrature.gausslobatto(n);
        τ =   convert(Vector{T}, τ_float64);
        w_q = convert(Vector{T}, w_q_float64);

        w_b = barycentric_weights(τ);
        return new{T}(n, τ, w_q, w_b)
    end
end