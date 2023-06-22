struct LegendreLobatto{T} <: Polynomials{T}
    n::Int
    τ::Vector{T}
    w_q::Vector{T}
    w_b::Vector{T}

    function LegendreLobatto(degree::Integer; type::Type{T}=Float64) where {T}

        n = degree + 1;
        τ, w_q = FGQ.gausslobatto(n);
        w_b = barycentric_weights(τ);

        return new{type}(n, τ, w_q, w_b)
    end
end