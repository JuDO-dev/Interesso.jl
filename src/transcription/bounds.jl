function discretize_bounds(dop::DOProblem, ::Transcription{T, I}, ctx::Context) where {T, I<:RigidIntervals}

    d_lower = Discretization{T}(undef, ctx);
    d_upper = Discretization{T}(undef, ctx);

    for i in 1:ctx.n_i
        for j in 1:ctx.n_x
            setindex_x_initial!(d_lower, dop.x[j].bounds.lower, i, j);
            setindex_x_initial!(d_upper, dop.x[j].bounds.upper, i, j);
            for k in 2:ctx.n_τ_x-1
                setindex_x_inner!(d_lower, dop.x[j].bounds.lower, i, j, k);
                setindex_x_inner!(d_upper, dop.x[j].bounds.upper, i, j, k);
            end
            setindex_x_final!(d_lower, dop.x[j].bounds.lower, i, j);
            setindex_x_final!(d_upper, dop.x[j].bounds.upper, i, j);
        end
        for j in 1:ctx.n_u
            for k in 1:ctx.n_τ_u
                setindex_u!(d_lower, dop.u[j].bounds.lower, i, j, k);
                setindex_u!(d_upper, dop.u[j].bounds.upper, i, j, k);
            end
        end
    end
    for j in 1:ctx.n_p
        setindex_p!(d_lower, dop.p[j].bounds.lower, j);
        setindex_p!(d_upper, dop.p[j].bounds.upper, j);
    end

    for j in 1:ctx.n_x
        setindex_x_initial!(d_lower, dop.x[j].initial.lower, 1, j);
        setindex_x_initial!(d_upper, dop.x[j].initial.upper, 1, j);
        setindex_x_final!(d_lower, dop.x[j].final.lower, ctx.n_i, j);
        setindex_x_final!(d_upper, dop.x[j].final.upper, ctx.n_i, j);
    end

    return Progradio.Box(d_lower.z, d_upper.z)
end