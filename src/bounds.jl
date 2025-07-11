abstract type AbstractBounds end
abstract type AbstractBoundsMesh end

function build_bounds_mesh end


# Exact Bounds

struct ExactBounds     <: AbstractBounds     end
struct ExactBoundsMesh <: AbstractBoundsMesh end

mesh_type(::Type{ExactBounds}) = ExactBoundsMesh

build_bounds_mesh(::ExactBounds, ::AbstractPoints) = ExactBoundsMesh()


# Sampled Bounds

struct SampledBounds{T<:AbstractPoints} <: AbstractBounds
    samp_points::T
end

SampledBounds(number::Integer) = SampledBounds(CGLPoints(number))

struct SampledBoundsMesh <: AbstractBoundsMesh 

    sampled_dif::Matrix{Float64}
    sampled_alg::Matrix{Float64}

    function SampledBoundsMesh(points::AbstractPoints, samp_points::AbstractPoints)

        bary_weights_dif = _barycentric_weights(points.points_dif_τ)
        bary_weights_alg = _barycentric_weights(points.points_alg_τ)

        sampled_dif = _interpolation_matrix(
            points.points_dif_τ,
            bary_weights_dif,
            samp_points.points_alg_τ,
        )

        sampled_alg = _interpolation_matrix(
            points.points_alg_τ,
            bary_weights_alg,
            samp_points.points_alg_τ,
        )

        return new(sampled_dif, sampled_alg)
    end

end

mesh_type(::Type{SampledBounds}) = SampledBoundsMesh

build_bounds_mesh(bounds::SampledBounds, points::AbstractPoints) = SampledBoundsMesh(points, bounds.samp_points)

get_points_samp_length(mesh::AbstractBoundsMesh) = size(mesh.sampled_alg, 1)


# Bernstein Bounds

struct BernsteinBounds     <: AbstractBounds end
struct BernsteinBoundsMesh <: AbstractBoundsMesh

    bernstein_dif::Matrix{Float64}
    bernstein_alg::Matrix{Float64}

    function BernsteinBoundsMesh(points_dif::Vector{<:Real}, points_alg::Vector{<:Real})

        if !all(0 .≤ points_dif .≤ 1)
            throw(DomainError(points_dif, "Please ensure 0 .≤ points_dif .≤ 1."))

        elseif !all(0 .≤ points_alg .≤ 1)
            throw(DomainError(points_alg, "Please ensure 0 .≤ points_alg .≤ 1."))

        else
            bernstein_dif = _bernstein(length(points_dif)) / _vandermonde(points_dif)
            bernstein_alg = _bernstein(length(points_alg)) / _vandermonde(points_alg)
        end

        return new(bernstein_dif, bernstein_alg)
    end
end

mesh_type(::Type{BernsteinBounds}) = BernsteinBoundsMesh

function build_bounds_mesh(::BernsteinBounds, points::AbstractPoints)
    
    return BernsteinBoundsMesh(
        0.5 * points.points_dif_τ .+ 0.5,
        0.5 * points.points_alg_τ .+ 0.5,
    )
end

_vandermonde(τ::Vector{<:Real}) = [τ_i^k for τ_i in τ, k in 0:(length(τ) - 1)]

_bernstein(n_τ::Integer) = [binomial(k, r) / binomial(n_τ - 1, r) for k in 0:(n_τ - 1), r in 0:(n_τ - 1)]