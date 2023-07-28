struct LeastSquares{F, M} <: InteressoTranscription{F, M}
    x::BarycentricPolynomial{F}
    u::BarycentricPolynomial{F}
    q::NestedQuadrature{F}
    mesh::M
end

function transcribe(dfp::DFProblem{F}, ls::LeastSquares{F, M}) where {F<:AbstractFloat, M<:RigidMesh}

    # Indexing
    n_x = length(dfp.x);
    n_u = length(dfp.u);
    l_x = ls.x.degree + 1;
    l_u = ls.u.degree + 1;
    indices = indexing(ls.mesh, n_x, n_u, l_x, l_u);

    # Quadrature
    order = ls.q.order;
    τ = ls.q.nodes;
    w = ls.q.weights;
    
    # Interpolation
    D = differentiation_matrix(ls.x);
    I_x = interpolation_matrix(ls.x, τ);
    I_u = interpolation_matrix(ls.u, τ);
    I = Interpolation(n_x, n_u, I_x * D, I_x, I_u);

    # Integrate residuals
    function integrate_residuals(z_i::Vector{F}, h_i::F) where {F<:AbstractFloat}
        return sum(w[j] * sum((dfp.residuals(h_i * Ẋ(I, z_i, j), X(I, z_i, j), U(I, z_i, j))).^2) for j in 1:order)
    end

    # Objective
    function f(z::Vector{F}) where {F<:AbstractFloat}
        fz = zero(F);
        for i in 1:ls.mesh.N
            fz += ls.mesh.h[i] * integrate_residuals(z[indices[i]], ls.mesh.h[i])
        end
        return fz
    end

    # Gradient
    function g!(gz::Vector{F}, z::Vector{F}) where {F<:AbstractFloat}
        gz[:] = zeros(F, length(z));
        for i in 1:ls.mesh.N
            gz[indices[i]] .+= ls.mesh.h[i] .* Enzyme.gradient(Reverse, z_i -> integrate_residuals(z_i, ls.mesh.h[i]), z[indices[i]]);
        end
        return nothing
    end

    # Bounds
    (z_ℓ, z_u) = bounds(dfp, ls, indices);

    # Guess
    z_0 = zeros(F, indices[end][end]);

    return BCProblem(f, g!, z_ℓ, z_u, z_0, 0.01)
end

function transcribe(dfp::DFProblem{F}, ls::LeastSquares{F, M}) where {F<:AbstractFloat, M<:FlexibleMesh}
    return #scp
end