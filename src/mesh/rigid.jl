struct RigidMesh{F} <: InteressoMesh{F}
    N::Integer
    h::Vector{F}
end

function RigidMesh(::Type{F}, N::Integer) where {F<:AbstractFloat}
    return RigidMesh{F}(N, repeat([one(F) / N], N))
end

RigidMesh(N::Integer) = RigidMesh(Float64, N);

function indexing(mesh::RigidMesh, n_x, n_u, l_x, l_u)
    mid = n_x * (l_x - 2) + n_u * l_u;
    return [(i-1)*(n_x+mid)+1 : i*(n_x+mid)+n_x for i in 1:mesh.N]
end