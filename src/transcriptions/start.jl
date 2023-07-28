JuDOInterface.DOSolution(dop::DOProblem{T, T0, TF, XS, US, E, C}) where {T, T0, TF, XS, US, E, C} = DOSolution(
    dop.t_0,
    dop.t_f,
    [x_m.start for x_m in dop.x],
    [u_m.start for u_m in dop.u],
    [p_m.start for p_m in dop.p],
    ConstantInterpolant{T}[],
    typemax(T),
    nothing
);

function discretize_start(start::DOSolution, ctx::Context, transcription::Transcription{T, RI}) where {T, RI<:RigidIntervals}

    z_start = Vector{T}(undef, ctx.p_offset + ctx.n_p);

    τ_x_i = Vector{T}(undef, ctx.n_τ_x);
    τ_u_i = Vector{T}(undef, ctx.n_τ_u);
    
    for i in 1:ctx.n_i
        
        dtdτ = (transcription.intervals.t[i+1] - transcription.intervals.t[i]) / 2;
        t̄_i = (transcription.intervals.t[i+1] + transcription.intervals.t[i]) / 2;
        @. τ_x_i = transcription.x_rule.τ * dtdτ + t̄_i;
        @. τ_u_i = transcription.u_rule.τ * dtdτ + t̄_i;

        for j in 1:ctx.n_x
            setindex_x_initial!(z_start, start.x[j](τ_x_i[1]), ctx, i, j);
            for k in 2:ctx.n_τ_x-1
                setindex_x_inner!(z_start, start.x[j](τ_x_i[k]), ctx, i, j, k);
            end
            setindex_x_final!(z_start, start.x[j](τ_x_i[end]), ctx, i, j);
        end
        for j in 1:ctx.n_u
            for k in 1:ctx.n_τ_u
                setindex_u!(z_start, start.u[j](τ_u_i[k]), ctx, i, j, k);
            end
        end
    end

    for j in 1:ctx.n_p
        setindex_p!(z_start, start.p[j], ctx, j);
    end

    return z_start
end