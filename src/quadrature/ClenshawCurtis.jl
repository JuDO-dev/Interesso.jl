struct ClenshawCurtis{F} <: InteressoQuadrature{F}
    order::Integer
    nodes::Vector{F}
    weights::Vector{F}
    function ClenshawCurtis{F}(order::Integer) where {F<:AbstractFloat}
        nodes = clenshawcurtisnodes(F, order);
        weights = clenshawcurtisweights(chebyshevmoments1(F, order));
        return new(order, [nodes[i] for i in order:-1:1], weights)
    end
end

ClenshawCurtis(order::Integer) = ClenshawCurtis{Float64}(order);