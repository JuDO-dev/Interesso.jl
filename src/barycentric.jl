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

function interpolation_vector(p::PointDistribution{T}, τ::T) where {T}

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

function differentiation_matrix(p::PointDistribution{T}) where {T}

    D = Matrix{T}(undef, length(p.τ), length(p.τ));

    for l in eachindex(p.τ)
        Σ = zero(T);
        for k in eachindex(p.τ)
            if k != l
                D[l,k] = (p.w_b[k] / p.w_b[l]) / (p.τ[l] - p.τ[k]);
                Σ += D[l,k];
            end
        end
        D[l,l] = -Σ;
    end

    return D
end