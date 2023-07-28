transcribe(dop::DOProblem, transcription::Transcription) = transcribe(dop, transcription, DOSolution(dop));

function transcribe(dfp::DFProblem, transcription::Transcription{T, I}, start::DOSolution) where {T, I}

    ctx = Context(dfp, transcription);
    
    z_start = discretize_start(start, ctx, transcription);
    d_start = PrimalDiscretization{T}(undef, ctx, length(dfp.equations); z=z_start);

    constraint_set = discretize_bounds(dfp, ctx, transcription);
    
    function _f(z::V, ẋ::VVV, x::VVV, u::VVV, p::V, t::V, r::VVV, r2::VV, ∫r2::V, Σ∫r2::V, ctx::Context)::T where {T,
        V<:Vector{T}, VV<:Vector{Vector{T}}, VVV<:Vector{Vector{Vector{T}}}}

        decode_x!(x, z, ctx, transcription);
        get_ẋ!(ẋ, x, z, ctx, transcription);
        decode_u!(u, z, ctx, transcription);
        decode_p!(p, z, ctx);
        decode_t!(t, z, ctx, transcription);
        
        evaluate_residuals!(r, r2, ∫r2, Σ∫r2, ẋ, x, u, p, t, ctx, dfp, transcription);

        return Σ∫r2[1]
    end

    f(d::PrimalDiscretization{T})::T where {T} = _f(d.z, d.ẋ, d.x, d.u, d.p, d.t, d.r, d.r2, d.∫r2, d.Σ∫r2, d.ctx);

    gd = PrimalDiscretization{T}(undef, ctx, length(dfp.equations));

    function g!(gz::Vector{T}, d::PrimalDiscretization{T}) where {T}

        gz .= 0;
        gd.z = gz;

        Enzyme.autodiff(Enzyme.Reverse, f, Enzyme.Active, Enzyme.Duplicated(d, gd));
        return nothing
    end
 
    return Progradio.NLProblem(d_start, constraint_set, f, g!)
end

function evaluate_residuals!(r::VVV, r2::VV, ∫r2::V, Σ∫r2::V,
    ẋ::VVV, x::VVV, u::VVV, p::V, t::V, ctx::Context, dop::DOProblem, transcription::Transcription) where {T,
    V<:Vector{T}, VV<:Vector{Vector{T}}, VVV<:Vector{Vector{Vector{T}}}}

    r! = get_residuals(dop.equations);

    for i in 1:ctx.n_i
        t̄_i = (t[i+1] + t[i]) / 2;
        Δt_i = t[i+1] - t[i];
        dtdτ = Δt_i / 2;
        dτdt = 2 / Δt_i;

        for q in 1:ctx.n_τ_q
            @. ẋ[i][q] *= dτdt;
            r!(r[i][q], ẋ[i][q], x[i][q], u[i][q], p, transcription.q_rule.τ[q] * dtdτ + t̄_i);
            r2[i][q] = mapreduce(r_j -> r_j^2, +, r[i][q]);
        end
        ∫r2[i] = sum(transcription.q_rule.w_q[q] * r2[i][q] for q in 1:ctx.n_τ_q) * Δt_i;
    end
    Σ∫r2[1] = sum(∫r2) / (2 * length(dop.equations) * (t[end] - t[1]));

    return nothing
end

#function transcribe(dop::DOProblem, transcription::Transcription, start::DOSolution) end

#function evaluate_cost!() end



get_residuals(equations::ResidualEquations) = equations.r!;



#=
    ctx = Context(ls.intervals.n, length(dfp.x), length(dfp.u), length(dfp.p), length(dfp.equations), ls.x_rule.n, ls.u.n, ls.q.n);
    
    d_start = discretize_start(start, ls, ctx);
    constraint_set = discretize_bounds(dfp, ls, ctx);

    r! = get_residuals(dfp.equations);
    integrated_residuals_scale = 1 / (2 * ctx.n_r * (dfp.t_f - dfp.t_0));



    function _f(z::Vector, ẋ, x, u, p,  ctx::Context)

        

        return 
    end




    



    function f(d::Discretization{T}) where {T}

        interpolate!(d, vẋ, vx, vu);
        
        for i in 1:d.ctx.n_i
            t̄_i = (ls.intervals.t[i+1] + ls.intervals.t[i]) / 2;
            Δt_i = ls.intervals.t[i+1] - ls.intervals.t[i];
            dtdτ = Δt_i / 2;
            dτdt = 2 / Δt_i;

            for q in eachindex(vẋ, vx, vu)
                @. d.ẋ[i][q] *= dτdt;
                r!(d.r.r[i][q], d.ẋ[i][q], d.x[i][q], d.u[i][q], d.p, ls.q.τ[q] * dtdτ + t̄_i);
                d.r.norm[i][q] = mapreduce(r_j -> r_j^2, +, d.r.r[i][q]);
            end
            d.r.norm_integrated[i] = sum(ls.q.w_q[q] * d.r.norm[i][q] for q in 1:d.ctx.n_τ_q) * Δt_i;
        end
        d.r.norm_integrated_sum = sum(d.r.norm_integrated) * integrated_residuals_scale;

        return d.r.norm_integrated_sum
    end

    gd = Discretization{T}(undef, d_start.ctx);

    function g!(gz::Vector{T}, d::Discretization{T}) where {T}
        gd.z = gz;
        Enzyme.gradient!(Enzyme.Reverse, gd, f, d);
        return nothing
    end

    return Progradio.NLProblem(d_start, constraint_set, f, g!)
end

#=function transcribe(dfp::DOProblem{T, R, R, GX, X, GU, U, Nothing, Nothing, E, Nothing},
    transcription::Transcription{T, I}, warmstart::Warmstart) where {T<:Real, GX, X, GU, U, E, I<:FlexibleIntervals{T}}
end

function transcribe(dop::DOProblem{T, R, R, XS, US, PS, E, C},
    transcription::Transcription{T, I}, warmstart::Warmstart) where {T<:Real, XS, US, PS, E, C, I<:RigidIntervals{T}}
end

function transcribe(dop::DOProblem{T, R, R, GX, X, GU, U, Nothing, Nothing, E, C},
    transcription::Transcription{T, I}, warmstart::Warmstart) where {T<:Real, GX, X, GU, U, E, C, I<:FlexibleIntervals{T}}
end=#

#=function integrate_residuals_i!(r̃::V, x̃̇::V, x̃::V, ũ::V, p::AbstractVector{T}, t::V,
    vẋ::VV, vx::VV, vu::VV, z::Discretization{T}, i::Integer, quadrature, residuals!) where {T<:Real, V<:Vector{T}, VV<:Vector{V}, CTX}
    
    t̄_i = (t[i] + t[i+1]) / 2;
    Δt_i = t[i+1] - t[i];
    dtdτ = Δt_i / 2;
    dτdt = 2 / Δt_i;
    
    Σ_i = zero(T);
    for q in eachindex(quadrature.τ, quadrature.w)
        
        interpolate_x!(x̃̇, vẋ[q], z, i);
        @. x̃̇ *= dτdt;
        interpolate_x!(x̃, vx[q], z, i);
        
        for j in
            interpolate_u!(ũ, vu[q], z, i);

        residuals!(r̃, x̃̇, x̃, ũ, p, quadrature.τ[q] * dtdτ + t̄_i);
        Σ_i += quadrature.w[q] * sum(r̃)^2;
    end

    return Σ_i * Δt_i
end=#







#=
function transcribe(dfp::DOProblem{F}, M::RigidMesh{F}, warmstart::Vector{IntervalData{F}}) where {F<:AbstractFloat}

    # warmstart and Bounds
    z_0, z_ℓ, z_u = encode(dfp, M, warmstart);
    
    # Indexing
    n_x = length(dfp.x);
    n_u = length(dfp.u);
    indices = indexing(M, n_x, n_u);
    
    # Interpolation
    D = differentiation_matrix(M.P_x);
    I_x = interpolation_matrix(M.P_x, M.Q.τ);
    I_u = interpolation_matrix(M.P_u, M.Q.τ);
    I = Interpolation(n_x, n_u, I_x * D, I_x, I_u);

    # Integrated residuals per interval
    Δt = dfp.t.final - dfp.t.initial;
    n_r = dfp.equations.n_r;
    function f_i(z_i::Vector{F}, h_i::F) where {F<:AbstractFloat}
        return sum(M.Q.w[j] * sum((dfp.equations.r(
            2 * Ẋ(I, z_i, j) ./ (h_i * Δt),
            X(I, z_i, j),
            U(I, z_i, j))).^2) for j in 1:M.Q.order) / (2 * n_r)
    end

    # Objective
    function f(z::Vector{F}) where {F<:AbstractFloat}
        fz = zero(F);
        for i in 1:length(M.H)
            fz += M.H[i] * f_i(z[indices[i]], M.H[i]);
        end
        return fz
    end

    # Gradient
    function g!(gz::Vector{F}, z::Vector{F}) where {F<:AbstractFloat}
        gz[:] = zeros(F, length(z));
        for i in 1:1:length(M.H)
            gz[indices[i]] .+= M.H[i] .* gradient(z_i -> f_i(z_i, M.H[i]), z[indices[i]])[1];
        end
        return nothing
    end

    return z_0, z_ℓ, z_u, f_i, f, g!
end=#=#