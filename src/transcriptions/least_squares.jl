struct LeastSquares{T, I, P, B} <: Transcription{T, I, P, B}
    intervals::I
    x_polynomials::P
    u_polynomials::P
    q_polynomials::P
    x_bounds::B
    u_bounds::B
    
    a_ẋ::Vector{Vector{T}}
    a_x::Vector{Vector{T}}
    a_u::Vector{Vector{T}}

    function LeastSquares(intervals::I, x_polynomials::P, u_polynomials::P, q_polynomials::P, x_bounds::B, u_bounds::B) where {T,
        I<:Intervals{T}, P<:Polynomials{T}, B<:Bounds}

        a_x = [barycentric_interpolation(x_polynomials.τ, x_polynomials.w_b, τ_j) for τ_j in q_polynomials.τ];
        D = barycentric_differentiation(x_polynomials.τ, x_polynomials.w_b);
        Dᵀ = permutedims(D, (2, 1));
        a_ẋ = [Dᵀ * a_x_j for a_x_j in a_x];
        a_u = [barycentric_interpolation(u_polynomials.τ, u_polynomials.w_b, τ_j) for τ_j in q_polynomials.τ];

        return new{T, I, P, B}(intervals, x_polynomials, u_polynomials, q_polynomials, x_bounds, u_bounds, a_ẋ, a_x, a_u)
    end
end

z_length(dop::DOProblem, transcription::Transcription) = transcription.intervals.n * (length(dop.x) * transcription.x_polynomials.n + length(dop.u) * transcription.u_polynomials.n);

mutable struct LeastSquaresEvaluator{T} <: MOI.AbstractNLPEvaluator



    function LeastSquaresEvaluator()

        return new{T}()
    end
end

function MOI.eval_objective(lse::LeastSquaresEvaluator, z::AbstractVector{T})::T where {T}

    @views begin
        z_views = [z[range] for range in lse.z_ranges];
    end

    for i in eachindex(z_views)
        lse.r_tape.tape.func(z_views[i]);
    end

    return 
end

function MOI.eval_objective_gradient(lse, grz, z)




    return nothing
end
#


















#=nlp_block_data(e::NLPEvaluator, ::LeastSquares) = MOI.NLPBlockData(MOI.NLPBoundsPair[], e, true);

function MOI.eval_objective(e::NLPEvaluator, z::AbstractVector{T}) where {T}

    for i in eachindex(e.caches)
        # Copy z slice
        copyto!(e.caches[i].z_i, 1, z, e.z_offsets[i], e.z_i_length);

        # Evaluate interval objective
        objective!(e.caches[i], e);
    end

    return sum(s.objective for s in e.caches);
end

function MOI.eval_objective_gradient(e::NLPEvaluator, gz::AbstractVector{T}, z::AbstractVector{T}) where {T}
    
    for i in eachindex(e.caches, e.adjoint_caches)
        # Seed storage for AD
        fill!(e.caches[i], zero(T));
        copyto!(e.caches[i].z_i, 1, z, e.ctx.z_offsets[i], e.ctx.z_i_length);
        e.caches[i].objective = 1.0;

        # Automatic differentiation
        Enzyme.autodiff(Enzyme.Reverse, objective!,
            Duplicated(e.caches[i], e.adjoint_caches[i]),
            Const(e)
        );

        # Paste back into gz
        copyto!(gz, e.ctx.z_offsets[i], e.adjoint_caches[i].z_i, 1, e.ctx.z_i_length);
    end
    return nothing
end

MOI.eval_constraint(::NLPEvaluator, ::AbstractVector, ::AbstractVector) = nothing;
MOI.eval_constraint_jacobian(::NLPEvaluator, ::AbstractVector, ::AbstractVector) = nothing;


function objective!(cache::Cache{T}, e::NLPEvaluator{T, PT, LS}) where {T, PT, LS<:LeastSquares}

    decode_ẋ_i!(cache.ẋ_i, cache.z_i, e);
    decode_x_i!(cache.x_i, cache.z_i, e);
    decode_u_i!(cache.u_i, cache.z_i, e);

    evaluate_r_i!(cache, e.problem);

    m_i!(cache)

    return nothing
end

function decode_ẋ!(ẋ::Vector{Vector{T}}, z_i::Vector{T}, e::NLPEvaluator{T, PT, LS}) where {T, PT, LS<:LeastSquares}

    for q in eachindex(ẋ)
        for j in eachindex(ẋ[q])
            ẋ[q][j] = sum(e.transcription.a_ẋ[q][k] * z_i[e.x_offsets[j]+k] for k in eachindex(e.transcription.a_ẋ[q]));
        end
    end
    return nothing
end

function decode_x!(x::Vector{Vector{T}}, z_i::Vector{T}, e::NLPEvaluator{T, PT, LS}) where {T, PT, LS<:LeastSquares}
    for q in eachindex(x)
        for j in eachindex(x[q])
            x[q][j] = sum(e.transcription.a_x[q][k] * z_i[e.x_offsets[j]+k] for k in eachindex(e.transcription.a_x[q]));
        end
    end
    return nothing
end

function decode_u!(u::Vector{Vector{T}}, z_i::Vector{T}, e::NLPEvaluator{T, PT, LS}) where {T, PT, LS<:LeastSquares}

    for q in eachindex(u)
        for j in eachindex(u[q])
            u[q][j] = sum(e.transcription.a_u[q][k] * z_i[e.u_offsets[j]+k] for k in eachindex(e.transcription.a_u[q]));
        end
    end
    return nothing
end














function evaluate_r_i!(cache, problem) where {T}

    for q in eachindex(cache.r)
        problem.r!(cache.r_i[q], ẋ_i[q], x_i[q], u_i[q]);
    end

    return nothing
end



=#





















































#=
function interval_objective!(d::Discretization{T}, v::AbstractVector{T}, t_a::T, t_b::T, problem, ls, ctx)
    
    decode_x!(d, ls, ctx);
    decode_u!(d, ls, ctx);

    evaluate_r!(d, v, t_a, t_b, problem);
    integrate_r!(d);

    return nothing
end

function decode_x!(d::Discretization{T}, ls::LeastSquares{T, I, XR, UR, QR}, ctx::Context) where {T, I, XR, UR, QR}
    for q in 1:ctx.n_τ_q
        for j in 1:ctx.n_x
            d.x[q][j] = sum(ls.vx[q][k] * getindex_x(d.z, ctx, j, k) for k in 1:ctx.n_τ_x);
            d.ẋ[q][j] = sum(ls.vẋ[q][k] * getindex_x(d.z, ctx, j, k) for k in 1:ctx.n_τ_x);
        end
    end
    return nothing
end

function decode_u!(d::Discretization{T}, ls::LeastSquares{T, I, XR, UR, QR}, ctx::Context) where {T, I, XR, UR, QR}
    for q in 1:ctx.n_τ_q
        for j in 1:ctx.n_u
            d.u[q][j] = sum(ls.vu[q][k] * getindex_u(d.z, ctx, j, k) for k in 1:ctx.n_τ_u);
        end
    end
    return nothing
end




MOI.NLPBlockData(problem::DOProblem, ls::LeastSquares, ctx::Context) = MOI.NLPBlockData([], InteressoNLPEvaluator(problem, ls, ctx), true);





getindex_x(z_i::AbstractVector{T}, j::Integer, k::Integer, ctx::Context, ::LeastSquares) where {T} = getindex(z_i, ctx.x_offsets[j] + k);
getindex_u(z_i::AbstractVector{T}, j::Integer, k::Integer, ctx::Context, ::LeastSquares) where {T} = getindex(z_i, ctx.u_offsets[j] + k);
getindex_λ(z_i::AbstractVector{T}, j::Integer, k::Integer, ctx::Context, ::LeastSquares) where {T} = getindex(z_i, ctx.λ_offsets[j] + k);

getindex_x(z::AbstractVector{T}, i::Integer, j::Integer, k::Integer, ctx::Context, ::LeastSquares) where {T} = getindex(z, ctx.z_offsets[i] + ctx.x_offsets[j] + k);
getindex_u(z::AbstractVector{T}, i::Integer, j::Integer, k::Integer, ctx::Context, ::LeastSquares) where {T} = getindex(z, ctx.z_offsets[i] + ctx.u_offsets[j] + k);
getindex_λ(z::AbstractVector{T}, i::Integer, j::Integer, k::Integer, ctx::Context, ::LeastSquares) where {T} = getindex(z, ctx.z_offsets[i] + ctx.λ_offsets[j] + k);

getindex_v(z::AbstractVector{T}, j::Integer, ctx::Context, ::LeastSquares) where {T} = getindex(z, ctx.v_offset + j);
getindex_t(z::AbstractVector{T}, i::Integer, ctx::Context, ::LeastSquares) where {T} = getindex(z, ctx.t_indices[i]);
=#











#=
function decode_x!(x::VVV, z::V, ctx::Context, ls::LeastSquares) where {T,
    V<:Vector{T}, VVV<:Vector{Vector{Vector{T}}}}

    for i in 1:ctx.n_i
        for q in 1:ctx.n_τ_q
            for j in 1:ctx.n_x
                x[i][q][j] = ls.vx[q][1] * getindex_x_initial(z, ctx, i, j) +
                         sum(ls.vx[q][k] * getindex_x_inner(z, ctx, i, j, k) for k in 2:ctx.n_τ_x-1) +
                           ls.vx[q][end] * getindex_x_final(z, ctx, i, j);
            end
        end
    end
    return nothing
end

function get_ẋ!(ẋ::VVV, ::VVV, z::V, ctx::Context, ls::LeastSquares) where {T,
    V<:Vector{T}, VVV<:Vector{Vector{Vector{T}}}}

    for i in 1:ctx.n_i
        for q in 1:ctx.n_τ_q
            for j in 1:ctx.n_x
                ẋ[i][q][j] = ls.vẋ[q][1] * getindex_x_initial(z, ctx, i, j) +
                         sum(ls.vẋ[q][k] * getindex_x_inner(z, ctx, i, j, k) for k in 2:ctx.n_τ_x-1) +
                           ls.vẋ[q][end] * getindex_x_final(z, ctx, i, j);
            end
        end
    end
    return nothing
end

function decode_u!(u::VVV, z::V, ctx::Context, ls::LeastSquares) where {T,
    V<:Vector{T}, VVV<:Vector{Vector{Vector{T}}}}

    for i in 1:ctx.n_i
        for q in 1:ctx.n_τ_q
            for j in 1:ctx.n_u
                u[i][q][j] = sum(ls.vu[q][k] * getindex_u(z, ctx, i, j, k) for k in 1:ctx.n_τ_u);
            end
        end
    end
    return nothing
end=#