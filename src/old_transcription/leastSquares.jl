struct LeastSquares{F, M} <: InteressoTranscription{F, M}
    x::InteressoPolynomial{F}
    u::InteressoPolynomial{F}
    q::InteressoQuadrature{F}
    mesh::M
end

# Rigid mesh
#function transcribe(dfp::DFProblem)




transcribe(dfp, ls::LeastSquares{F, M}) where {F<:AbstractArray, M<:RigidMesh}















function transcribe(dfp::DFProblem{F}, mesh::RigidMesh{F}, ls::LeastSquares) where F<:AbstractFloat

    #Quadrature points and weights
    quadrature_weights
    interval_indexing = indexing(dfp, mesh);

    function f(z::Vector{F})::F 
        fz = zero(F);
        for indices in interval_indexing
            fz[indices] += rigid_residuals_integration(z[indices], A, Q);
        end
        return fz
    end

    function g!(gz::Vector{F}, z::Vector{F})
        gz[:] = zeros()
        for indices in interval_indexing
            gz[indices] += Zygote.gradient(z_i -> rigid_residuals_integration(z_i, A, Q), z[indices]);
        end
        return nothing
    end

    z_ℓ, z_u = bounds(dfp, rm, ls);
    z_0 = initial(dfp, rm, ls);

    return BCProblem{F}()
end

function rigid_residuals_integration(z_i, w_q)
    integrated_residuals = zeros();

    for w_q in w
        residuals[:] = zeros();
        dfp.residuals!(residuals, ẋ, x, u);
        integrated_residuals += w_q * sum(residuals .^ 2);
    end

    return integrated_residuals
end


function transcribe(dfp::DFProblem, fm::FlexibleMesh, ls::LeastSquares)
    return BCProblem()
end



function interval_indexing(problem::JuDOProblem, mesh::RigidMesh)
    return 
end



function interval_indexing(problem::JuDOProblem, mesh::FlexibleMesh)
    return 
end

function bounds(dfp::DFProblem{F}, rm::RigidMesh{F}, ls::LeastSquares) where F<:AbstractFloat
    n_x = 
    n_u =
    z_ℓ = Vector{F}(undef, n);
    z_u = Vector{F}(undef, n);
 
    return (z_ℓ, z_u)
end

function initial()

    z_0 = similar(z_ℓ);

    return z_0
end




r(x_dot, x, u, p) = 0;

z = [];

x = view(z, x_bit_vector)


r(view(z[indices], x_bits), x)
