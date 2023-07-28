mutable struct Discretization{T}
    z::Vector{T}
    
    ẋ::Vector{Vector{T}}
    x::Vector{Vector{T}}
    u::Vector{Vector{T}}
    λ::Vector{Vector{T}}

    r::Vector{Vector{T}}
    Σr2::Vector{T}
    ∫Σr2::T

    c::Vector{T}
    ∫c::T

    objective::T
end

function slice_z!(d::Vector{PrimalDiscretization{T}}, z::Vector{T}, ctx::Context) where {T}

    for i in 1:ctx.n_i

        @. d[i].z = z[]
    end

    return nothing
end

function evaluate_r!(d::Discretization, v, t_a, t_b, problem)


    return nothing
end

function integrate_r!(d::Discretrization)


    return nothing
end



#=

 <: AbstractVector{T} end

Base.size(d::Discretization) = size(d.z);
Base.setindex!(d::Discretization, X, inds...) = setindex!(d.z, X, inds...);
Base.getindex(d::Discretization, inds...) = getindex(d.z, inds...);

mutable struct PrimalDiscretization{T} <: Discretization{T}
    z::Vector{T}
    
    ẋ::Vector{Vector{Vector{T}}}
    x::Vector{Vector{Vector{T}}}
    u::Vector{Vector{Vector{T}}}
    p::Vector{T}
    t::Vector{T}

    r::Vector{Vector{Vector{T}}}
    r2::Vector{Vector{T}}
    ∫r2::Vector{T}
    Σ∫r2::Vector{T}

    ctx::Context
end

PrimalDiscretization{T}(::UndefInitializer, ctx::Context, n_r::Integer;
    z::Vector{T}=Vector{T}(undef, ctx.z_length)) where {T} = PrimalDiscretization(
    
    z,
    
    [[Vector{T}(undef, ctx.n_x) for _ in 1:ctx.n_τ_q] for _ in 1:ctx.n_i],
    [[Vector{T}(undef, ctx.n_x) for _ in 1:ctx.n_τ_q] for _ in 1:ctx.n_i],
    [[Vector{T}(undef, ctx.n_u) for _ in 1:ctx.n_τ_q] for _ in 1:ctx.n_i],
    Vector{T}(undef, ctx.n_p),
    Vector{T}(undef, ctx.n_i + 1),
    
    [[Vector{T}(undef, n_r) for _ in 1:ctx.n_τ_q] for _ in 1:ctx.n_i],
    [Vector{T}(undef, ctx.n_τ_q) for _ in 1:ctx.n_i],
    Vector{T}(undef, ctx.n_i),
    Vector{T}(undef, 1),

    ctx
);

#=mutable struct PrimalDualDiscretization{T} <: Discretization{T}
    z::Vector{T}
    
    ẋ::Vector{Vector{Vector{T}}}
    x::Vector{Vector{Vector{T}}}
    u::Vector{Vector{Vector{T}}}
    λ::Vector{Vector{Vector{T}}}
    p::Vector{T}
    t::Vector{T}

    r::Vector{Vector{Vector{T}}}
    r2::Vector{Vector{T}}
    ∫r2::Vector{T}
    Σ∫r2::Vector{T}

    c::Vector{Vector{T}}
    ∫c::Vector{T}
    Σ∫c::Vector{T}

    ctx::Context
end=#

function JuDOInterface.DOSolution(dfp::FixedTimeDFProblem, transcription::Transcription, d::PrimalDiscretization) 
    
    x̂ = get_x(d); û = get_u(d); p = get_p(d); r̂ = get_r(d);
    
    t = transcription.intervals.t;
    n_i = transcription.intervals.n;
    intervals = [Interval(t[i], t[i+1]) for i in 1:n_i];
    t̄ =         [(t[i+1] + t[i]) / 2    for i in 1:n_i];
    dtdτ =      [(t[i+1] - t[i]) / 2    for i in 1:n_i];
    
    t_x = [transcription.x_rule.τ .* dtdτ[i] .+ t̄[i] for i in 1:n_i];
    t_u = [transcription.u_rule.τ .* dtdτ[i] .+ t̄[i] for i in 1:n_i];
    t_q = [transcription.q_rule.τ .* dtdτ[i] .+ t̄[i] for i in 1:n_i];

    x̃ = [PiecewiseInterpolant([LagrangeInterpolant(intervals[i], t_x[i], x̂[i][j], transcription.x_rule.w_b
    ) for i in 1:d.ctx.n_i]) for j in 1:d.ctx.n_x];
    ũ = [PiecewiseInterpolant([LagrangeInterpolant(intervals[i], t_u[i], û[i][j], transcription.u_rule.w_b
    ) for i in 1:d.ctx.n_i]) for j in 1:d.ctx.n_u];
    r̃ = [PiecewiseInterpolant([LagrangeInterpolant(intervals[i], t_q[i], r̂[i][j], transcription.q_rule.w_b
    ) for i in 1:d.ctx.n_i]) for j in 1:d.ctx.n_r];

    return DOSolution(dfp.t_0, dfp.t_f, x̃, ũ, p, r̃, d.Σ∫r2[1], nothing)
end

function get_x(d::Discretization{T}) where {T}

    x̂ = [[Vector{T}(undef, d.ctx.n_τ_x) for _ in 1:d.ctx.n_x] for _ in 1:d.ctx.n_i];

    for i in 1:d.ctx.n_i
        for j in 1:d.ctx.n_x
            x[i][j][1] = getindex_x_initial(d.z, d.ctx, i, j);
            for k in 2:d.ctx.n_τ_x-1
                x[i][j][k] = getindex_x_inner(d.z, d.ctx, i, j, k);
            end
            x[i][j][end] = getindex_x_final(d.z, d.ctx, i, j);
        end
    end
    return x̂
end

function get_u(d::Discretization{T}) where {T}

    û = [[Vector{T}(undef, d.ctx.n_τ_u) for _ in 1:d.ctx.n_u] for _ in 1:d.ctx.n_i];

    for i in 1:d.ctx.n_i
        for j in 1:d.ctx.n_u
            for k in 1:d.ctx.n_τ_u
                u[i][j][k] = getindex_u(d.z, d.ctx, i, j, k);
            end
        end
    end
    return û
end

function get_p(d::Discretization{T}) where {T}

    p = Vector{T}(undef, d.ctx.n_p);
    
    decode_p!(p, d.z, d.ctx);
    return p
end

function get_r(d::Discretization{T}) where {T}

    r̂ = [[Vector{T}(undef, d.ctx.n_τ_q) for _ in 1:d.ctx.n_r] for _ in 1:d.ctx.n_i];

    for i in 1:d.ctx.n_i
        for j in 1:d.ctx.n_r
            for q in 1:d.ctx.n_τ_q
                r[i][j][q] = d.r[i][q][j];
            end
        end
    end

    return r̂
end

function decode_p!(p::V, z::V, ctx::Context) where {T, V<:Vector{T}}

    for j in 1:ctx.n_p
        p[j] = getindex_p(z, ctx, j);
    end
    return nothing
end

function decode_t!(t::V, ::V, ::Context, transcription::Transcription{T, RI}) where {T,
    V<:Vector{T}, RI<:RigidIntervals{T}}

    @. t = transcription.intervals.t;
    return nothing
end=#