struct Interpolation{F<:AbstractFloat}
    n_x::Integer
    n_u::Integer
    l_x::Integer
    l_u::Integer
    I_ẋ::Matrix{F}
    I_x::Matrix{F}
    I_u::Matrix{F}
end

function Interpolation(n_x::Integer, n_u::Integer, I_ẋ::Matrix{F}, I_x::Matrix{F}, I_u::Matrix{F}) where {F<:AbstractFloat}
    l_x = size(I_x, 2);
    l_u = size(I_u, 2);
    return Interpolation{F}(n_x, n_u, l_x, l_u, I_ẋ, I_x, I_u)
end
#=
function Ẋ!(ẋ::AbstractVector{F}, I::Interpolation{F}, z_i::AbstractVector{F}, j::Integer) where {F<:AbstractFloat}
    n_x = I.n_x;
    l_x = I.l_x;
    offset = n_x * (l_x - 1) + I.n_u * I.l_u;
    for n in 1:n_x
        ẋ[n] =      I.I_ẋ[j,1]      * z_i[n];
        ẋ[n] += dot(I.I_ẋ[j,2:end-1], z_i[n_x+(n-1)*(l_x-2)+1 : n_x+n*(l_x-2)]);
        ẋ[n] +=     I.I_ẋ[j,end]    * z_i[offset+n];
    end
    return nothing
end

function Ẋ(I::Interpolation{F}, z_i::AbstractVector{F}, j::Integer) where {F<:AbstractFloat}
    ẋ = Vector{F}(undef, I.n_x);
    X!(ẋ, I, z_i, j);
    return ẋ
end
=#
function Ẋ(I::Interpolation{F}, z_i::AbstractVector{F}, j::Integer) where {F<:AbstractFloat}
    ẋ = Vector{F}(undef, I.n_x);
    n_x = I.n_x;
    l_x = I.l_x;
    offset = n_x * (l_x - 1) + I.n_u * I.l_u;
    for n in 1:n_x
        ẋ[n] =      I.I_ẋ[j,1]      * z_i[n];
        ẋ[n] += dot(I.I_ẋ[j,2:end-1], z_i[n_x+(n-1)*(l_x-2)+1 : n_x+n*(l_x-2)]);
        ẋ[n] +=     I.I_ẋ[j,end]    * z_i[offset+n];
    end
    return ẋ
end
#=
function X!(x::AbstractVector{F}, I::Interpolation{F}, z_i::AbstractVector{F}, j::Integer) where {F<:AbstractFloat}
    n_x = I.n_x;
    l_x = I.l_x;
    offset = n_x * (l_x - 1) + I.n_u * I.l_u;
    for n in 1:n_x
        x[n] =      I.I_x[j,1]      * z_i[n];
        x[n] += dot(I.I_x[j,2:end-1], z_i[n_x+(n-1)*(l_x-2)+1 : n_x+n*(l_x-2)]);
        x[n] +=     I.I_x[j,end]    * z_i[offset+n];
    end
    return nothing
end

function X(I::Interpolation{F}, z_i::AbstractVector{F}, j::Integer) where {F<:AbstractFloat}
    x = Vector{F}(undef, I.n_x);
    X!(x, I, z_i, j);
    return x
end
=#
function X(I::Interpolation{F}, z_i::AbstractVector{F}, j::Integer) where {F<:AbstractFloat}
    x = Vector{F}(undef, I.n_x);
    n_x = I.n_x;
    l_x = I.l_x;
    offset = n_x * (l_x - 1) + I.n_u * I.l_u;
    for n in 1:n_x
        x[n] =      I.I_x[j,1]      * z_i[n];
        x[n] += dot(I.I_x[j,2:end-1], z_i[n_x+(n-1)*(l_x-2)+1 : n_x+n*(l_x-2)]);
        x[n] +=     I.I_x[j,end]    * z_i[offset+n];
    end
    return x
end
#=
function U!(u::AbstractVector{F}, I::Interpolation{F}, z_i::AbstractVector{F}, j::Integer) where {F<:AbstractFloat}
    n_u = I.n_u;
    l_u = I.l_u;
    offset = I.n_x * (I.l_x - 1);
    for n in 1:n_u
        u[n] = dot(I.I_u[j,:], z_i[offset+(n-1)*(l_u)+1 : offset+n*l_u]);
    end
    return nothing
end

function U(I::Interpolation{F}, z_i::AbstractVector{F}, j::Integer) where {F<:AbstractFloat}
    u = Vector{F}(undef, I.n_u);
    U!(u, I, z_i, j);
    return u
end
=#
function U(I::Interpolation{F}, z_i::AbstractVector{F}, j::Integer) where {F<:AbstractFloat}
    u = Vector{F}(undef, I.n_u);
    n_u = I.n_u;
    l_u = I.l_u;
    offset = I.n_x * (I.l_x - 1);
    for n in 1:n_u
        u[n] = dot(I.I_u[j,:], z_i[offset+(n-1)*(l_u)+1 : offset+n*l_u]);
    end
    return u
end

# Interpolation Matrices [legacy]
