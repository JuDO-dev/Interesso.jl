function barycentric_weights(τ::Vector{T})::Vector{T} where {T}
    
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

function barycentric_interpolation(τ::Vector{T}, w_b::Vector{T}, τ_hat::T)::Vector{T} where {T}

    n = length(τ);
    v = Vector{T}(undef, n);

    exact = indexin(τ_hat, τ);
    if exact[1] isa Integer
        for l in eachindex(τ)
            v[l] = l == exact[1] ? one(T) : zero(T);
        end
    else
        denominator = zero(T);
        for l in eachindex(τ, w_b)
            common_l = w_b[l] / (τ_hat - τ[l]);
            v[l] = common_l;
            denominator += common_l;
        end

        one_over_denominator = 1 / denominator;
        @. v *= one_over_denominator;
        #for l in eachindex(τ)
        #    v[l] *= one_over_denominator;
        #end
    end

    return v
end

function barycentric_differentiation(τ::Vector{T}, w_b::Vector{T})::Matrix{T} where {T}

    n = length(τ);
    D = Matrix{T}(undef, n, n);

    for l in eachindex(τ)
        Σ = zero(T);
        for k in eachindex(τ)
            if k != l
                D[l,k] = (w_b[k] / w_b[l]) / (τ[l] - τ[k]);
                Σ += D[l,k];
            end
        end
        D[l,l] = -Σ;
    end

    return D
end