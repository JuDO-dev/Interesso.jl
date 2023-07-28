struct QuadratureRule{T<:Real, R<:Rule}
    n::Int
    τ::Vector{T}
    w_q::Vector{T}
    w_b::Vector{T}

    function QuadratureRule(τ::V, w_q::V, w_b::V, R::Type{<:Rule}) where {T, V<:Vector{T}}

        n = length(τ);
        length(w_q) == n  ? nothing : throw(ArgumentError("Ensure "));
        for j in eachindex(τ, w_q)
            -1 ≤ τ[j] ≤ 1 ? nothing : throw(DomainError(τ[j], "Ensure "));
            0 ≤ w_q[j]      ? nothing : throw(DomainError(w_q[j], "Ensure "))  
        end

        return new{T, R}(n, τ, w_q, w_b)
    end
end

function QuadratureRule(n::Integer, R::Type{LegendreLobatto}; T::Type=Float64)

    τ_float64, w_q_float64 = FastGaussQuadrature.gausslobatto(n);
    τ = convert(Vector{T}, τ_float64);
    w_q = convert(Vector{T}, w_q_float64);
    w_b = barycentric_weights(τ); 

    return QuadratureRule(τ, w_q, w_b, R)
end