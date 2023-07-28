mutable struct Cache{T}
    ẋ::Vector{T}
    x::Vector{T}
    u::Vector{T}
    r::Vector{T}
    c::T
end

mutable struct IntervalCache{T}
    z_i::Vector{T}
    caches::Cache{T}
    r_i::Vector{T}
    





mutable struct Cache{T}
    z_i_j::Vector{T}
    
    ẋ_i::Vector{Vector{T}}
    x_i::Vector{Vector{T}}
    u_i::Vector{Vector{T}}
    t::Vector{T}
    
    r_i::Vector{Vector{T}}
    
    c_0::T
    c_f::T
    c_i::Vector{T}

    objective_i::T

    Cache(z_i_length::Integer, n_q::Integer, n_x::Integer, n_u::Integer, n_r::Integer; T::Type=Float64) = new{T}(
        Vector{T}(undef, z_i_length),
        [Vector{T}(undef, n_x) for _ in 1:n_q],
        [Vector{T}(undef, n_x) for _ in 1:n_q],
        [Vector{T}(undef, n_u) for _ in 1:n_q],
        [Vector{T}(undef, n_r) for _ in 1:n_q],
        Vector{T}(undef, n_q)
    );
end

function Base.fill!(cache::Cache{T}, a::T) where {T}
    
    fill!(cache.z_i, a);
    
    for ẋ_i_q in cache.ẋ_i fill!(ẋ_i_q, a); end
    for x_i_q in cache.x_i fill!(x_i_q, a); end
    for u_i_q in cache.u_i fill!(u_i_q, a); end

    for r_i_q in cache.r_i fill!(r_i_q, a); end
    fill!(cache.c_i, a);

    return nothing
end