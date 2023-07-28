function interpolate!(d::Discretization{T}, vẋ::VV, vx::VV, vu::VV) where {T<:Real, CTX, VV<:Vector{Vector{T}}}

    for i in 1:d.ctx.n_i
        for q in eachindex(vẋ, vx, vu)
            for j in 1:d.ctx.n_x
                d.ẋ[i][q][j] = vẋ[q][1] * getindex_x_initial(d, i, j) +
                           sum(vẋ[q][k] * getindex_x_inner(d, i, j, k) for k in 2:d.ctx.n_τ_x-1) +
                             vẋ[q][end] * getindex_x_final(d, i, j);
                d.x[i][q][j] = vx[q][1] * getindex_x_initial(d, i, j) +
                           sum(vx[q][k] * getindex_x_inner(d, i, j, k) for k in 2:d.ctx.n_τ_x-1) +
                             vx[q][end] * getindex_x_final(d, i, j);
            end
            for j in 1:d.ctx.n_u
                d.u[i][q][j] = sum(vu[q][k] * getindex_u(d, i, j, k) for k in 1:d.ctx.n_τ_u);
            end
        end
    end
    for j in 1:d.ctx.n_p
        d.p[j] = getindex_p(d, j);
    end
    return nothing
end

function interpolation_vectors(ls::LeastSquares)

    Dᵀ = transposed_differentiation_matrix(ls.x);
    vx = [interpolation_vector(ls.x, τ_j) for τ_j in ls.quadrature.τ];
    vẋ = [Dᵀ * vx_j for vx_j in vx];
    vu = [interpolation_vector(ls.u, τ_j) for τ_j in ls.quadrature.τ];

    return vẋ, vx, vu
end





#=function interpolate_x!(x̃::V, v_q::V, z::Discretization{T}, i::Integer) where {T<:Real, V<:AbstractVector{T}, CTX}

    for j in 1:z.ctx.n_x
        Σ_j = zero(T);
        Σ_j +=     v_q[1] * getindex_x_initial(z, i, j);
        Σ_j += sum(v_q[k] * getindex_x_inner(z, i, j, k) for k in 2:z.ctx.n_τ_x-1);
        Σ_j +=   v_q[end] * getindex_x_final(z, i, j);
        x̃[j] = Σ_j;
    end
    return nothing
end

function interpolate_u!(ũ::V, v_q::V, z::Discretization{T}, i::Integer) where {T<:Real, V<:AbstractVector{T}, CTX}

    for j in 1:z.ctx.n_u
        ũ[j] = sum(v_q[k] * getindex_u(z, i, j, k) for k in 1:z.ctx.n_τ_u);
    end
    return nothing
end=#



#interpolation_vectors(::DirectCollocation) = (nothing, nothing, nothing);


#=function interpolate_x!(xz::Vector{T}, x̃v_j::Vector{T}, z::Vector{T},
    z_offset::Integer, n_x::Integer, n_τ_x::Integer, x_offset::Integer) where {T<:Real}

    for m in 1:n_x
        Σ_m = zero(T);
        Σ_m += x̃v_j[1] * z[z_offset+m];
        for l in 2:n_τ_x-1
            Σ_m += x̃v_j[l] * z[z_offset+n_x+(m-1)*(n_τ_x-2)+l];
        end
        Σ_m += x̃v_j[end] * z[z_offset+x_offset+m]
        xz[m] = Σ_m;     
    end
    return nothing
end=#

#=function interpolate_u!(uz::Vector{T}, ũv_j::Vector{T}, z::Vector{T},
    z_offset::Integer, n_u::Integer, n_τ_u::Integer, u_offset::Integer) where {T<:Real}
    
    for m in 1:n_u
        Σ_m = zero(T);
        for l in 1:n_τ_u
            Σ_m += ũv_j[l] * z[z_offset+u_offset+(m-1)*(n_τ_u)+l];
        end
        uz[m] = Σ_m;
    end
    return nothing
end=#

#=
function X!(xz::Vector{T}, interpolation::Vector{T}, z_i::Vector{T},
    n_x::Integer, n_τ_x::Integer, x_offset::Integer) where {T<:Real}

    for m in 1:n_x
        xz[m] =      interpolation[1] * z_i[m];
        xz[m] += sum(interpolation[l] * z_i[n_x+(m-1)*(n_τ_x-2)+l] for l in 2:n_τ_x-1);
        xz[m] +=   interpolation[end] * z_i[x_offset+m];
    end
    return nothing
end

function U!(uz::Vector{T}, interpolation::Vector{T}, z_i::Vector{T},
    n_u::Integer, n_τ_u::Integer, u_offset::Integer) where {T<:Real}
    
    for m in 1:n_u
        uz[m] = sum(interpolation[l] * z_i[u_offset+(m-1)*(n_τ_u)+l] for l in 1:n_τ_u);
    end
    return nothing
end


struct Interpolation{F<:AbstractFloat}
    n_x::Integer
    n_u::Integer
    n_τ_x::Integer
    n_τ_u::Integer
    I_ẋ::Matrix{F}
    I_x::Matrix{F}
    I_u::Matrix{F}
end

function Interpolation(n_x::Integer, n_u::Integer, I_ẋ::Matrix{F}, I_x::Matrix{F}, I_u::Matrix{F}) where {F<:AbstractFloat}
    n_τ_x = size(I_x, 2);
    n_τ_u = size(I_u, 2);
    return Interpolation{F}(n_x, n_u, n_τ_x, n_τ_u, I_ẋ, I_x, I_u)
end
#=
function Ẋ!(ẋ::AbstractVector{F}, I::Interpolation{F}, z_i::AbstractVector{F}, j::Integer) where {F<:AbstractFloat}
    n_x = I.n_x;
    n_τ_x = I.n_τ_x;
    offset = n_x * (n_τ_x - 1) + I.n_u * I.n_τ_u;
    for n in 1:n_x
        ẋ[n] =      I.I_ẋ[j,1]      * z_i[n];
        ẋ[n] += dot(I.I_ẋ[j,2:end-1], z_i[n_x+(n-1)*(n_τ_x-2)+1 : n_x+n*(n_τ_x-2)]);
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
    n_τ_x = I.n_τ_x;
    offset = n_x * (n_τ_x - 1) + I.n_u * I.n_τ_u;
    for n in 1:n_x
        ẋ[n] =      I.I_ẋ[j,1]      * z_i[n];
        ẋ[n] += dot(I.I_ẋ[j,2:end-1], z_i[n_x+(n-1)*(n_τ_x-2)+1 : n_x+n*(n_τ_x-2)]);
        ẋ[n] +=     I.I_ẋ[j,end]    * z_i[offset+n];
    end
    return ẋ
end

function rrule(::typeof(Ẋ), I, z_i, j)
    function Ẋ_pullback(ẋ̄)
        n_x, n_τ_x, n_u, n_τ_u = I.n_x, I.n_τ_x, I.n_u, I.n_τ_u;
        z̄_i = zeros(n_x * n_τ_x + n_u * n_τ_u);
        offset = n_x * (n_τ_x - 1) + n_u * n_τ_u;
        for n in 1:n_x
            z̄_i[n]                                   =  ẋ̄[n] *  I.I_ẋ[j,1];
            z̄_i[n_x+(n-1)*(n_τ_x-2)+1 : n_x+n*(n_τ_x-2)] .= ẋ̄[n] .* I.I_ẋ[j,2:end-1];
            z̄_i[offset+n]                            =  ẋ̄[n] *  I.I_ẋ[j,end];
        end
        return ZeroTangent(), ZeroTangent(), z̄_i, ZeroTangent()
    end
    return Ẋ(I, z_i, j), Ẋ_pullback
end

#=
function X!(x::AbstractVector{F}, I::Interpolation{F}, z_i::AbstractVector{F}, j::Integer) where {F<:AbstractFloat}
    n_x = I.n_x;
    n_τ_x = I.n_τ_x;
    offset = n_x * (n_τ_x - 1) + I.n_u * I.n_τ_u;
    for n in 1:n_x
        x[n] =      I.I_x[j,1]      * z_i[n];
        x[n] += dot(I.I_x[j,2:end-1], z_i[n_x+(n-1)*(n_τ_x-2)+1 : n_x+n*(n_τ_x-2)]);
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
    n_τ_x = I.n_τ_x;
    offset = n_x * (n_τ_x - 1) + I.n_u * I.n_τ_u;
    for n in 1:n_x
        x[n] =      I.I_x[j,1]      * z_i[n];
        x[n] += dot(I.I_x[j,2:end-1], z_i[n_x+(n-1)*(n_τ_x-2)+1 : n_x+n*(n_τ_x-2)]);
        x[n] +=     I.I_x[j,end]    * z_i[offset+n];
    end
    return x
end

function rrule(::typeof(X), I, z_i, j)
    function X_pullback(x̄)
        n_x, n_τ_x, n_u, n_τ_u = I.n_x, I.n_τ_x, I.n_u, I.n_τ_u;
        z̄_i = zeros(n_x * n_τ_x + n_u * n_τ_u);
        offset = n_x * (n_τ_x - 1) + n_u * n_τ_u;
        for n in 1:n_x
            z̄_i[n]                                   =  x̄[n] *  I.I_x[j,1];
            z̄_i[n_x+(n-1)*(n_τ_x-2)+1 : n_x+n*(n_τ_x-2)] .= x̄[n] .* I.I_x[j,2:end-1];
            z̄_i[offset+n]                            =  x̄[n] *  I.I_x[j,end];
        end
        return ZeroTangent(), ZeroTangent(), z̄_i, ZeroTangent()
    end
    return X(I, z_i, j), X_pullback
end

#=
function U!(u::AbstractVector{F}, I::Interpolation{F}, z_i::AbstractVector{F}, j::Integer) where {F<:AbstractFloat}
    n_u = I.n_u;
    n_τ_u = I.n_τ_u;
    offset = I.n_x * (I.n_τ_x - 1);
    for n in 1:n_u
        u[n] = dot(I.I_u[j,:], z_i[offset+(n-1)*(n_τ_u)+1 : offset+n*n_τ_u]);
    end
    return nothing
end

function U(I::Interpolation{F}, z_i::AbstractVector{F}, j::Integer) where {F<:AbstractFloat}
    u = Vector{F}(undef, I.n_u);
    U!(u, I, z_i, j);
    return u
end
=#


function rrule(::typeof(U), I, z_i, j)
    function U_pullback(ū)
        n_x, n_τ_x, n_u, n_τ_u = I.n_x, I.n_τ_x, I.n_u, I.n_τ_u;
        z̄_i = zeros(n_x * n_τ_x + n_u * n_τ_u);
        offset = n_x * (n_τ_x - 1);
        for n in 1:n_u
            z̄_i[offset+(n-1)*(n_τ_u)+1 : offset+n*n_τ_u] .= ū[n] .* I.I_u[j,:];
        end
        return ZeroTangent(), ZeroTangent(), z̄_i, ZeroTangent()
    end
    return U(I, z_i, j), U_pullback
end
=#