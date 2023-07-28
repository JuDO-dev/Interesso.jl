struct GaussKronrod{F} <: NestedQuadrature{F}
    order::Integer
    nodes::Vector{F}
    weights::Vector{F}
    higher_nodes::Vector{F}
    higher_weights::Vector{F}
    function GaussKronrod{F}(order::Integer) where {F<:AbstractFloat}
        (GK_nodes, K_weights, G_weights) = kronrod(order);
        higher_nodes = vcat(GK_nodes[1:end], -GK_nodes[end-1:-1:1]);
        higher_weights = vcat(K_weights[1:end], K_weights[end-1:-1:1]);
        nodes = higher_nodes[2:2:end];
        if iseven(order)
            weights = vcat(G_weights[1:end], G_weights[end:-1:1]);
        else 
            weights = vcat(G_weights[1:end], G_weights[end-1:-1:1]);
        end
        return new(order, nodes, weights, higher_nodes, higher_weights)
    end
end

GaussKronrod(::Type{F}, order::Integer) where {F<:AbstractFloat} = GaussKronrod{F}(order);

GaussKronrod(order::Integer) = GaussKronrod(Float64, order);