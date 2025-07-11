"""
    AbstractPoints

Supertype for point distributions in the ``[-1, 1]`` interval.
"""
abstract type AbstractPoints end

"""
    AbstractPointsMesh

Supertype for point distributions in an interval ``[t_a, t_b]``.
"""
abstract type AbstractPointsMesh end

"""
    build_points_mesh(::AbstractPoints, t_a::Real, t_b::Real)::AbstractPointsMesh

"""
function build_points_mesh(::AbstractPoints, ::Real, ::Real)::AbstractPointsMesh end

get_points_dif_length(mesh::AbstractPointsMesh) = length(mesh.points_dif)
get_points_alg_length(mesh::AbstractPointsMesh) = length(mesh.points_alg)
get_points_quad_length(mesh::AbstractPointsMesh) = length(mesh.points_alg)

## Legendre-Gauss-Radau

"""
    LGRPoints(size::Integer)


"""
struct LGRPoints <: AbstractPoints
    points_dif_τ::Vector{Float64}
    points_alg_τ::Vector{Float64}
    quad_weights_τ::Vector{Float64}

    function LGRPoints(number::Integer)
        
        if !(number ≥ 1)
            throw(DomainError("Please ensure number ≥ 1."))
        end

        points_alg_τ, quad_weights_τ = FGQ.gaussradau(number)
        points_dif_τ = vcat(points_alg_τ, 1.0)

        return new(points_dif_τ, points_alg_τ, quad_weights_τ)
    end
end

"""
    LGRPointsMesh(::LGRPoints, t_a::Real, t_b::Real)


"""
struct LGRPointsMesh <: AbstractPointsMesh

    t_a::Float64
    t_b::Float64

    points_dif::Vector{Float64}
    points_alg::Vector{Float64}

    quad_weights::Vector{Float64}

    bary_weights_dif::Vector{Float64}
    bary_weights_alg::Vector{Float64}

    differentiation::Matrix{Float64}

    function LGRPointsMesh(points::LGRPoints, t_a::Real, t_b::Real)

        _throw_if_invalid_bounds(t_a, t_b)

        Δt = t_b - t_a
        Σt = t_a + t_b

        points_dif = [0.5 * Δt * p_j .+ 0.5 * Σt for p_j in points.points_dif_τ]
        points_alg = [0.5 * Δt * p_j .+ 0.5 * Σt for p_j in points.points_alg_τ]

        quad_weights = [0.5 * Δt * w_j for w_j in points.quad_weights_τ]

        bary_weights_dif = _barycentric_weights(points_dif)
        bary_weights_alg = _barycentric_weights(points_alg)
        
        differentiation = _differentiation_matrix(points_dif, bary_weights_dif)

        return new(t_a, t_b, points_dif, points_alg, quad_weights, bary_weights_dif, 
            bary_weights_alg, differentiation,
        )
    end
end

mesh_type(::Type{LGRPoints}) = LGRPointsMesh

build_points_mesh(points::LGRPoints, t_a::Real, t_b::Real) = LGRPointsMesh(points, t_a, t_b)


## Gauss-Legendre

struct GLPoints <: AbstractPoints
    points_alg_τ::Vector{Float64}
    quad_weights_τ::Vector{Float64}

    function GLPoints(number::Integer)
        
        if !(number ≥ 1)
            throw(DomainError("Please ensure number ≥ 1."))
        end

        points_alg_τ, quad_weights_τ = FGQ.gausslegendre(number)

        return new(points_alg_τ, quad_weights_τ)
    end
end


struct GLPointsMesh <: AbstractPointsMesh

    t_a::Float64
    t_b::Float64

    points_alg::Vector{Float64}

    quad_weights::Vector{Float64}

    function GLPointsMesh(points::GLPoints, t_a::Real, t_b::Real)

        _throw_if_invalid_bounds(t_a, t_b)

        Δt = t_b - t_a
        Σt = t_a + t_b

        points_alg = [0.5 * Δt * p_j .+ 0.5 * Σt for p_j in points.points_alg_τ]

        quad_weights = [0.5 * Δt * w_j for w_j in points.quad_weights_τ]

        return new(t_a, t_b, points_alg, quad_weights)
    end
end

mesh_type(::Type{GLPoints}) = GLPointsMesh

build_method_mesh(points::GLPoints, mesh::AbstractPointsMesh) = GLPointsMesh(points, mesh.t_a, mesh.t_b)


##  Chebyshev-Gauss-Lobatto

struct CGLPoints <: AbstractPoints
    points_alg_τ::Vector{Float64}

    function CGLPoints(number::Integer)
        
        if !(number ≥ 1)
            throw(DomainError("Please ensure number ≥ 1."))
        end

        points_alg_τ = cos.((number-1:-1:0) .* π ./ (number-1))

        return new(points_alg_τ)
    end
end

function _throw_if_invalid_bounds(lower::Real, upper::Real)

    if lower ≥ upper
        return throw(ArgumentError("Please ensure lower < upper."))
    end
end