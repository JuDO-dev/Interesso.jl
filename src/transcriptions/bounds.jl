function add_bounds!(m::MOI.ModelLike, x::VVV, u::VVV, λ::VVV, v::V, t::V, ctx, dop) where {}

    for i in 1:ctx.n_i
        for j in 1:ctx.n_x
            MOI.add_constraint(m, x[i][j], MOI.Interval(dop.x[j].bounds.lower, dop.x[j].bounds.upper));
        end
        for j in 1:ctx.n_u
            MOI.add_constraint(m, u[i][j], MOI.Interval(dop.u[j].bounds.lower, dop.u[j].bounds.upper));
        end
        for j in 1:ctx.n_λ
            MOI.add_constraint(m, λ[i][j], MOI.Interval(dop.λ[j].bounds.lower, dop.λ[j].bounds.upper));
        end
    end

    for j in 1:ctx.n_v
        MOI.add_constraint(m, v[i][j], MOI.Interval(dop.v[j].bounds.lower, dop.v[j].bounds.upper));
    end

    for j in 1:ctx.n_t
        MOI.add_constraint(m, t[i][j], MOI.Interval(dop.t[j].bounds.lower, dop.t[j].bounds.upper));
    end

    return nothing
end














function discretize_bounds(dop::DOProblem, ctx::Context, ::Transcription{T, RI}) where {T, RI<:RigidIntervals}

    z_lower = Vector{T}(undef, ctx.z_length);
    z_upper = Vector{T}(undef, ctx.z_length);

    for i in 1:ctx.n_i
        for j in 1:ctx.n_x
            setindex_x_initial!(z_lower, dop.x[j].bounds.lower, ctx, i, j);
            setindex_x_initial!(z_upper, dop.x[j].bounds.upper, ctx, i, j);
            for k in 2:ctx.n_τ_x-1
                setindex_x_inner!(z_lower, dop.x[j].bounds.lower, ctx, i, j, k);
                setindex_x_inner!(z_upper, dop.x[j].bounds.upper, ctx, i, j, k);
            end
            setindex_x_final!(z_lower, dop.x[j].bounds.lower, ctx, i, j);
            setindex_x_final!(z_upper, dop.x[j].bounds.upper, ctx, i, j);
        end
        for j in 1:ctx.n_u
            for k in 1:ctx.n_τ_u
                setindex_u!(z_lower, dop.u[j].bounds.lower, ctx, i, j, k);
                setindex_u!(z_upper, dop.u[j].bounds.upper, ctx, i, j, k);
            end
        end
    end
    for j in 1:ctx.n_p
        setindex_p!(z_lower, dop.p[j].bounds.lower, ctx, j);
        setindex_p!(z_upper, dop.p[j].bounds.upper, ctx, j);
    end

    for j in 1:ctx.n_x
        setindex_x_initial!(z_lower, dop.x[j].initial.lower, ctx, 1, j);
        setindex_x_initial!(z_upper, dop.x[j].initial.upper, ctx, 1, j);
        setindex_x_final!(z_lower, dop.x[j].final.lower, ctx, ctx.n_i, j);
        setindex_x_final!(z_upper, dop.x[j].final.upper, ctx, ctx.n_i, j);
    end

    return Progradio.Box(z_lower, z_upper)
end