struct LeastSquares{T, I, XPD<:PointDistribution{T},
    UPD<:PointDistribution{T}, QPD<:PointDistribution{T}} <: Transcription{T, I}
    
    intervals::I
    x::XPD
    u::UPD
    q::QPD
end

function interpolation_vectors(ls::LeastSquares)

    Dᵀ = transpose(differentiation_matrix(ls.x));
    vx = [interpolation_vector(ls.x, τ_j) for τ_j in ls.q.τ];
    vẋ = [Dᵀ * vx_j for vx_j in vx];
    vu = [interpolation_vector(ls.u, τ_j) for τ_j in ls.q.τ];

    return vẋ, vx, vu
end

function interpolate!(d::Discretization{T}, vẋ::VV, vx::VV, vu::VV) where {T, VV<:Vector{Vector{T}}}

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