struct Legendre{F} <: InteressoPolynomial{F}
    degree::Integer
    nodes::Vector{F}
    weights::Vector{F}
    function Legendre{F}(degree::Integer) where {F<:AbstractFloat} 
        (nodes, weights) = gauss(degree + 1);
        for d in 1:(degree+1)
            weights[d] = ((-1)^(d - 1)) * sqrt((1 - nodes[d]^2) * weights[d]);
        end
        return new(degree, nodes, weights)
    end
end

Legendre(::Type{F}, degree::Integer) where {F<:AbstractFloat} = Legendre{F}(degree);

Legendre(degree::Integer) = Legendre(Float64, degree);