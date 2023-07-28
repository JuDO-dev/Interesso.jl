struct PolynomialInterpolation{T, R<:Rule}
    n::Int
    τ::Vector{T}
    w_b::Vector{T}

    function PolynomialInterpolation(τ::Vector{T}, w_b::Vector{T}, R::Type{<:Rule}) where {T}

        length(τ) == length(w_b) ? nothing : throw(ArgumentError("Ensure"));
        for τ_j in τ
            -1 ≤ τ_j ≤ 1         ? nothing : throw(DomainError(τ_j, "Ensure"));
        end

        return new{T, R}(length(τ), τ, w_b)
    end
end

function PolynomialInterpolation(degree::Integer, R::Type{LegendreLobatto}; T::Type=Float64)

    τ_float64, _ = FastGaussQuadrature.gausslobatto(degree + 1);
    τ = convert(Vector{T}, τ_float64);

    w_b = barycentric_weights(τ);

    return PolynomialInterpolation(τ, w_b, R)
end

function barycentric_weights(τ::Vector{T}) where {T}
    
    n = length(τ);
    w_b = Vector{T}(undef, n);
    for l in eachindex(τ, w_b)
        denominator_l = one(T);
        for k in 1:l-1
            denominator_l *= τ[l] - τ[k];
        end
        for k in l+1:n
            denominator_l *= τ[l] - τ[k];
        end
        w_b[l] = 1 / denominator_l;
    end
    return w_b
end
    


#=function PolynomialInterpolation(degree::Integer, R::Type{ChebyshevSecond}; T::Type{<:Real}=Float64)
    
    n = degree + 1;
    τ = Vector{T}(undef, n);
    w_b = Vector{T}(undef, n);

    for l in eachindex(τ, w_b)
        τ[l] = cospi((n - l) / degree);
        w_b[l] = ((-1)^(n - l)) * (1 - 0.5 * (l == 1 || l == n));
    end

    return PolynomialInterpolation(τ, w_b, R)
end=#

function interpolation_vector(p::PolynomialInterpolation, τ::T) where {T}

    v = Vector{T}(undef, p.n);
    
    exact = indexin(τ, p.τ);
    if exact[1] isa Integer
        for l in eachindex(p.τ)
            v[l] = l == exact[1] ? one(T) : zero(T);
        end
    else
        denominator = zero(T);
        for l in eachindex(p.τ)
            common_l = p.w_b[l] / (τ - p.τ[l]);
            v[l] = common_l;
            denominator += common_l;
        end

        one_over_denominator = 1 / denominator;
        for l in eachindex(p.τ)
            v[l] *= one_over_denominator;
        end
    end

    return v
end

function transposed_differentiation_matrix(p::PolynomialInterpolation{T}) where {T}

    D = Matrix{T}(undef, length(p.τ), length(p.τ));

    for l in eachindex(p.τ)
        Σ = zero(T);
        for k in eachindex(p.τ)
            if k != l
                D[k,l] = (p.w_b[k] / p.w_b[l]) / (p.τ[l] - p.τ[k]);
                Σ += D[k,l];
            end
        end
        D[l,l] = -Σ;
    end

    return D
end