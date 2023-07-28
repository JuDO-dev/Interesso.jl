struct Context
    n_i::Int
    n_x::Int
    n_u::Int
    n_τ_x::Int
    n_τ_u::Int
    n_τ_q::Int
    z_length::Int
    z_i_length::Int
    z_offsets::Vector{Int}
    x_offsets::Vector{Int}
    u_offsets::Vector{Int}
end

function Context(dfp::FixedTimeDFProblem, transcription::LeastSquares)

    n_i = transcription.intervals.n;
    n_x = length(dfp.x);
    n_u = length(dfp.u);
    n_τ_x = transcription.x_polynomials.n;
    n_τ_u = transcription.u_polynomials.n;
    n_τ_q = transcription.q_polynomials.n;

    z_offsets = [(i - 1) * x_offset for i in 1:n_i];

    u_offsets = [u_offset + (j - 1) * n_τ_u for j in 1:n_u];

    return Context(
        
        
        
        )
end





setindex_x!(z_i::AbstractVector{T}, X::T, ctx::Context, j::Integer, k::Integer) where {T} = setindex!(z_i, X, ctx.x_offsets[j] + k);
setindex_u!(z_i::AbstractVector{T}, X::T, ctx::Context, j::Integer, k::Integer) where {T} = setindex!(z_i, X, ctx.u_offsets[j] + k);
setindex_λ!(z_i::AbstractVector{T}, X::T, ctx::Context, j::Integer, k::Integer) where {T} = setindex!(z_i, X, ctx.λ_offsets[j] + k);
setindex_v!(z::AbstractVector{T}, X::T, ctx::Context, j::Integer) where {T} = setindex!(z, X, ctx.p_offset + j);
setindex_t!(z::AbstractVector{T}, X::T, ctx::Context, i::Integer) where {T} = setindex!(z, X, ctx.t_indices[i]);

getindex_x(z_i::AbstractVector{T}, ctx::Context, j::Integer, k::Integer) where {T} = getindex(z_i, ctx.x_offsets[j] + k);
getindex_u(z_i::AbstractVector{T}, ctx::Context, j::Integer, k::Integer) where {T} = getindex(z_i, ctx.u_offsets[j] + k);
getindex_λ(z_i::AbstractVector{T}, ctx::Context, j::Integer, k::Integer) where {T} = getindex(z_i, ctx.λ_offsets[j] + k);
getindex_v(z::AbstractVector{T}, ctx::Context, j::Integer) where {T} = getindex(z, ctx.p_offset + j);
getindex_t(z::AbstractVector{T}, ctx::Context, i::Integer) where {T} = getindex(z, ctx.t_indices[i]);


function Context(dfp::FixedTimeDOProblem, transcription::Transcription)

    n_i = transcription.intervals.n;
    n_x = length(dfp.x);
    n_u = length(dfp.u);
    n_λ = 0;
    n_p = length(dfp.p);
    n_τ_x = transcription.x_rule.n;
    n_τ_u = transcription.u_rule.n;
    n_τ_q = transcription.q_rule.n;
    
    u_offset = n_x * (n_τ_x - 1);
    x_offset = u_offset + n_u * n_τ_u;
    p_offset = x_offset * n_i + n_x;
    t_offset = 0;

    z_offsets = [(i - 1) * x_offset for i in 1:n_i];
    x_offsets = [n_x - 1 + (j - 1) * (n_τ_x - 2) for j in 1:n_x];
    u_offsets = [u_offset + (j - 1) * n_τ_u for j in 1:n_u];
    λ_offsets = Int[];

    z_length = p_offset + n_p;

    return Context(z_length, n_i, n_x, n_u, n_λ, n_p, n_τ_x, n_τ_u, n_τ_q,
        x_offset, p_offset, t_offset,
        z_offsets, x_offsets, u_offsets, λ_offsets
    );
end

setindex_x_initial!(z::Vector{T}, X::T, ctx::Context, i::Integer, j::Integer) where {T} = setindex!(
    z, X, ctx.z_offsets[i] + j
);

getindex_x_initial(z::Vector{T}, ctx::Context, i::Integer, j::Integer) where {T} = getindex(
    z, ctx.z_offsets[i] + j
);

setindex_x_inner!(z::Vector{T}, X::T, ctx::Context, i::Integer, j::Integer, k::Integer) where {T} = setindex!(
    z, X, ctx.z_offsets[i] + ctx.x_offsets[j] + k
);

getindex_x_inner(z::Vector{T}, ctx::Context, i::Integer, j::Integer, k::Integer) where {T} = getindex(
    z, ctx.z_offsets[i] + ctx.x_offsets[j] + k
);

setindex_x_final!(z::Vector{T}, X::T, ctx::Context, i::Integer, j::Integer) where {T} = setindex!(
    z, X, ctx.z_offsets[i] + ctx.x_offset + j
);

getindex_x_final(z::Vector{T}, ctx::Context, i::Integer, j::Integer) where {T} = getindex(
    z, ctx.z_offsets[i] + ctx.x_offset + j
);

setindex_u!(z::Vector{T}, X::T, ctx::Context, i::Integer, j::Integer, k::Integer) where {T} = setindex!(
    z, X, ctx.z_offsets[i] + ctx.u_offsets[j] + k
);

getindex_u(z::Vector{T}, ctx::Context, i::Integer, j::Integer, k::Integer) where {T} = getindex(
    z, ctx.z_offsets[i] + ctx.u_offsets[j] + k
);

setindex_λ!(z::Vector{T}, X::T, ctx::Context, i::Integer, j::Integer, k::Integer) where {T} = setindex!(
    z, X, ctx.z_offsets[i] + ctx.λ_offsets[j] + k
);

getindex_λ(z::Vector{T}, ctx::Context, i::Integer, j::Integer, k::Integer) where {T} = getindex(
    z, ctx.z_offsets[i] + ctx.λ_offsets[j] + k
);

setindex_p!(z::Vector{T}, X::T, ctx::Context, j::Integer) where {T} = setindex!(
    z, X, ctx.p_offset + j
);

getindex_p(z::Vector{T}, ctx::Context, j::Integer) where {T} = getindex(
    z, ctx.p_offset + j
);

setindex_t!(z::Vector{T}, X::T, ctx::Context, i::Integer) where {T} = setindex!(
    z, X, ctx.z_offsets[i] + ctx.t_offset
);

getindex_t(z::Vector{T}, ctx::Context, i::Integer) where {T} = getindex(
    z, ctx.z_offsets[i] + ctx.t_offset
);