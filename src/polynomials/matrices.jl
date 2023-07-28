function differentiation_matrix(p::BarycentricPolynomial{F}) where {F<:AbstractFloat}
    D = Matrix{F}(undef, p.degree + 1, p.degree + 1);
    for i in 1:(p.degree + 1)
        sum = zero(F);
        for j in 1:(p.degree + 1)
            if j != i
                D[i,j] = (p.weights[j] / p.weights[i]) / (p.nodes[i] - p.nodes[j]);
                sum += D[i,j];
            end
        end
        D[i,i] = -sum;
    end
    return D
end

function interpolation_matrix(p::BarycentricPolynomial{F}, τ::Vector{F}) where {F<:AbstractFloat}
    n_p = p.degree + 1;
    n_τ = length(τ);
    A = Matrix{F}(undef, n_τ, n_p);
    for i in 1:n_τ
        exact = 0;
        sum = zero(F);
        for j in 1:n_p
            difference = τ[i] - p.nodes[j];
            exact = ifelse(difference == 0, j, exact);
            A[i,j] = p.weights[j] / difference;
            sum += A[i,j];
        end
        if sum == 0
            for j in 1:n_p
                A[i,j] = zero(F);
            end
        elseif exact > 0
            for j in 1:n_p
                A[i,j] = j == exact ? one(F) : zero(F);
            end
        else
            for j in 1:n_p
                A[i,j] /= sum;
            end
        end
    end
    return A
end