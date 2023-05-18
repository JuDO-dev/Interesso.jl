struct Context
    n_i::Int
    n_x::Int
    n_u::Int
    n_p::Int
    n_r::Int
    n_τ_x::Int
    n_τ_u::Int
    n_τ_q::Int
    x_offset::Int
    u_offset::Int
    p_offset::Int
    z_offsets::Vector{Int}
    x_offsets::Vector{Int}
    u_offsets::Vector{Int}
end

mutable struct ResidualsDiscretization{T}
    r::Vector{Vector{Vector{T}}}
    norm::Vector{Vector{T}}
    norm_integrated::Vector{T}
    norm_integrated_sum::T
end

mutable struct CostDiscretization{T}
    c::Vector{Vector{T}}
    integrated::Vector{T}
    integrated_sum::T
    c_initial::T
    c_final::T
    total::T
end

mutable struct Discretization{T} <: AbstractVector{T}
    z::Vector{T}
    ẋ::Vector{Vector{Vector{T}}}
    x::Vector{Vector{Vector{T}}}
    u::Vector{Vector{Vector{T}}}
    p::Vector{T}
    r::ResidualsDiscretization{T}
    c::CostDiscretization{T}
    ctx::Context
end

Base.size(d::Discretization) = size(d.z);
Base.setindex!(d::Discretization, X, inds...) = setindex!(d.z, X, inds...);
Base.getindex(d::Discretization, inds...) = getindex(d.z, inds...);

Discretization{T}(::UndefInitializer, ctx::Context) where {T} = Discretization(
    Vector{T}(undef, ctx.p_offset + ctx.n_p),
    [[Vector{T}(undef, ctx.n_x) for _ in 1:ctx.n_τ_q] for _ in 1:ctx.n_i],
    [[Vector{T}(undef, ctx.n_x) for _ in 1:ctx.n_τ_q] for _ in 1:ctx.n_i],
    [[Vector{T}(undef, ctx.n_u) for _ in 1:ctx.n_τ_q] for _ in 1:ctx.n_i],
    Vector{T}(undef, ctx.n_p),
    ResidualsDiscretization{T}(undef, ctx),
    CostDiscretization{T}(undef, ctx),
    ctx
);

ResidualsDiscretization{T}(::UndefInitializer, ctx::Context) where {T} = ResidualsDiscretization(
    [[Vector{T}(undef, ctx.n_r) for _ in 1:ctx.n_τ_q] for _ in 1:ctx.n_i],
    [Vector{T}(undef, ctx.n_τ_q) for _ in 1:ctx.n_i],
    Vector{T}(undef, ctx.n_i), zero(T)
);

CostDiscretization{T}(::UndefInitializer, ctx::Context) where {T} = CostDiscretization(
    [Vector{T}(undef, ctx.n_τ_q) for _ in 1:ctx.n_i],
    Vector{T}(undef, ctx.n_i), zero(T), zero(T), zero(T), zero(T)
);

function Context(n_i::I, n_x::I, n_u::I, n_p::I, n_r::I, n_τ_x::I, n_τ_u::I, n_τ_q::I) where {I<:Integer}

    u_offset = n_x * (n_τ_x - 1);
    x_offset = u_offset + n_u * n_τ_u;
    p_offset = x_offset * n_i + n_x;

    z_offsets = [(i - 1) * x_offset for i in 1:n_i];
    x_offsets = [n_x + (j - 1) * (n_τ_x - 2) - 1 for j in 1:n_x];
    u_offsets = [u_offset + (j - 1) * n_τ_u for j in 1:n_u];

    return Context(n_i, n_x, n_u, n_p, n_r, n_τ_x, n_τ_u, n_τ_q,
        x_offset, u_offset, p_offset, z_offsets, x_offsets, u_offsets
    )
end

function JuDOInterface.DOSolution(dop::FixedTimeDOProblem, transcription::Transcription, d::Discretization{T}) where {T} 
    
    x̂ = [[Vector{T}(undef, d.ctx.n_τ_x) for _ in 1:d.ctx.n_x] for _ in 1:d.ctx.n_i];
    û = [[Vector{T}(undef, d.ctx.n_τ_u) for _ in 1:d.ctx.n_u] for _ in 1:d.ctx.n_i];
    p =   Vector{T}(undef, d.ctx.n_p);
    r̂ = [[Vector{T}(undef, d.ctx.n_τ_q) for _ in 1:d.ctx.n_r] for _ in 1:d.ctx.n_i];
    
    decode_x!(x̂, d); decode_u!(û, d); decode_p!(p, d); decode_r!(r̂, d);
    
    t = transcription.intervals.t;
    n_i = transcription.intervals.n;
    intervals = [Interval(t[i], t[i+1]) for i in 1:n_i];
    t̄ =         [(t[i+1] + t[i]) / 2    for i in 1:n_i];
    dtdτ =      [(t[i+1] - t[i]) / 2    for i in 1:n_i];
    
    t_x = [transcription.x.τ .* dtdτ[i] .+ t̄[i] for i in 1:n_i];
    t_u = [transcription.u.τ .* dtdτ[i] .+ t̄[i] for i in 1:n_i];
    t_q = [transcription.q.τ .* dtdτ[i] .+ t̄[i] for i in 1:n_i];

    x̃ = [PiecewiseInterpolant([LagrangeInterpolant(intervals[i], t_x[i], x̂[i][j], transcription.x.w_b
    ) for i in 1:d.ctx.n_i]) for j in 1:d.ctx.n_x];
    ũ = [PiecewiseInterpolant([LagrangeInterpolant(intervals[i], t_u[i], û[i][j], transcription.u.w_b
    ) for i in 1:d.ctx.n_i]) for j in 1:d.ctx.n_u];
    r̃ = [PiecewiseInterpolant([LagrangeInterpolant(intervals[i], t_q[i], r̂[i][j], transcription.q.w_b
    ) for i in 1:d.ctx.n_i]) for j in 1:d.ctx.n_r];

    return DOSolution(dop.t_0, dop.t_f, x̃, ũ, p, r̃, d.r.norm_integrated_sum, nothing)
end

function decode_x!(x::Vector{Vector{Vector{T}}}, d::Discretization{T}) where {T}

    for i in 1:d.ctx.n_i
        for j in 1:d.ctx.n_x
            x[i][j][1] = getindex_x_initial(d, i, j);
            for k in 2:d.ctx.n_τ_x-1
                x[i][j][k] = getindex_x_inner(d, i, j, k);
            end
            x[i][j][end] = getindex_x_final(d, i, j);
        end
    end
    return nothing
end

function decode_u!(u::Vector{Vector{Vector{T}}}, d::Discretization{T}) where {T}

    for i in 1:d.ctx.n_i
        for j in 1:d.ctx.n_u
            for k in 1:d.ctx.n_τ_u
                u[i][j][k] = getindex_u(d, i, j, k);
            end
        end
    end
    return nothing
end

function decode_p!(p::Vector{T}, d::Discretization{T}) where {T}
    
    for j in 1:d.ctx.n_p
        p[j] = getindex_p(d, j);
    end
    return nothing
end

function decode_r!(r::Vector{Vector{Vector{T}}}, d::Discretization{T}) where {T}
    
    for i in 1:d.ctx.n_i
        for j in 1:d.ctx.n_r
            for k in 1:d.ctx.n_τ_q
                r[i][j][k] = d.r.r[i][k][j];
            end
        end
    end
    return nothing
end

function setindex_x_initial!(d::Discretization, X, i::Integer, j)

    return setindex!(d.z, X, d.ctx.z_offsets[i] + j)
end

function getindex_x_initial(d::Discretization, i::Integer, j)

    return getindex(d.z, d.ctx.z_offsets[i] + j)
end

function setindex_x_inner!(d::Discretization, X, i::Integer, j::Integer, k)

    return setindex!(d.z, X, d.ctx.z_offsets[i] + d.ctx.x_offsets[j] + k)
end

function getindex_x_inner(d::Discretization, i::Integer, j::Integer, k)

    return getindex(d.z, d.ctx.z_offsets[i] + d.ctx.x_offsets[j] + k)
end

function setindex_x_final!(d::Discretization, X, i::Integer, j)

    return setindex!(d.z, X, d.ctx.z_offsets[i] + d.ctx.x_offset + j)
end

function getindex_x_final(d::Discretization, i::Integer, j)

    return getindex(d.z, d.ctx.z_offsets[i] + d.ctx.x_offset + j)
end

function setindex_u!(d::Discretization, X, i::Integer, j::Integer, k)

    return setindex!(d.z, X, d.ctx.z_offsets[i] + d.ctx.u_offsets[j] + k)
end

function getindex_u(d::Discretization, i::Integer, j::Integer, k)

    return getindex(d.z, d.ctx.z_offsets[i] + d.ctx.u_offsets[j] + k)
end

function setindex_p!(d::Discretization, X, j)

    return setindex!(d.z, X, d.ctx.z_offsets[i] + j)
end

function getindex_p(d::Discretization, j)
    
    return getindex(d.z, d.ctx.p_offset + j)
end