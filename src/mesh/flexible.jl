struct FlexibleMesh{F} <: InteressoMesh{F}
    N::Integer
    h_min::F
    h_guess::Vector{F}
end

function FlexibleMesh(::Type{F}, N::Integer, Δτ_min::AbstractFloat) where {F<:AbstractFloat}
    return FlexibleMesh{F}(N, convert(F, Δτ_min), repeat([one(F) / N], N))
end

FlexibleMesh(N::Integer, h_min::AbstractFloat) = FlexibleMesh(Float64, N, h_min);