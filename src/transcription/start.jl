function discretize_start(start::DOSolution, transcription::Transcription{T, I}, ctx::Context) where {T, I<:RigidIntervals}

    d_start = Discretization{T}(undef, ctx);

    τ_x_i = Vector{T}(undef, ctx.n_τ_x);
    τ_u_i = Vector{T}(undef, ctx.n_τ_u);
    
    for i in 1:ctx.n_i
        
        dtdτ = (transcription.intervals.t[i+1] - transcription.intervals.t[i]) / 2;
        t̄_i = (transcription.intervals.t[i+1] + transcription.intervals.t[i]) / 2;
        @. τ_x_i = transcription.x.τ * dtdτ + t̄_i;
        @. τ_u_i = transcription.u.τ * dtdτ + t̄_i;

        for j in 1:ctx.n_x
            setindex_x_initial!(d_start, start.x[j](τ_x_i[1]), i, j);
            for k in 2:ctx.n_τ_x-1
                setindex_x_inner!(d_start, start.x[j](τ_x_i[k]), i, j, k);
            end
            setindex_x_final!(d_start, start.x[j](τ_x_i[end]), i, j);
        end
        for j in 1:ctx.n_u
            for k in 1:ctx.n_τ_u
                setindex_u!(d_start, start.u[j](τ_u_i[k]), i, j, k);
            end
        end
    end

    for j in 1:ctx.n_p
        setindex_p!(d_start, start.p[j], j);
    end

    return d_start
end