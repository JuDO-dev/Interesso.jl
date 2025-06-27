abstract type AbstractBounds end
abstract type AbstractBoundsMesh end

function build_bounds_mesh end


# Sampled Bounds

struct SampledBounds     <: AbstractBounds     end
struct SampledBoundsMesh <: AbstractBoundsMesh end

mesh_type(::Type{SampledBounds}) = SampledBoundsMesh

build_bounds_mesh(::SampledBounds, ::AbstractPoints) = SampledBoundsMesh()

# Bernstein Bounds

struct BernsteinBounds     <: AbstractBounds end
struct BernsteinBoundsMesh <: AbstractBoundsMesh

    bernstein_dif::Matrix{Float64}
    bernstein_alg::Matrix{Float64}

    function BernsteinBoundsMesh(points_dif::Vector{<:Real}, points_alg::Vector{<:Real})

        if !all(0 .≤ points_dif .≤ 1) == true
            throw(DomainError(points_dif, "Please ensure 0 .≤ points_dif`` .≤ 1."))

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