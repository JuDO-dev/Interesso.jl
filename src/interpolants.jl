"""
    AbstractInterpolant

Supertype for Interesso interpolants.
    
It is a subtype of
[`DOI.AbstractDynamicSolution`](@extref DynOptInterface.AbstractDynamicSolution).
"""
abstract type AbstractInterpolant <: DOI.AbstractDynamicSolution end

struct PiecewiseInterpolant{I<:AbstractInterpolant} <: AbstractInterpolant
    initial::Float64
    final::Float64
    pieces::Vector{I}
    pieces_initials::Vector{Float64}

    function PiecewiseInterpolant(pieces::Vector{I}) where {I<:AbstractInterpolant}
        for i in 2:length(pieces)
            if !(pieces[i].initial ≈ pieces[i-1].final)
                throw(ArgumentError("Please ensure the pieces are continuous."))
            end
        end
        
        return new{I}(
            first(pieces).initial,
            last(pieces).final,
            pieces,
            [piece.initial for piece in pieces],
        )
    end
end

function (piecewise::PiecewiseInterpolant)(t::Real)
    if piecewise.initial ≤ t ≤ piecewise.final
        return piecewise.pieces[searchsortedlast(piecewise.pieces_initials, t)](t)
    else
        return NaN
    end
end

"""
    LagrangeInterpolant(
        initial::T,
        final::T,
        points::Vector{T},
        weights::Vector{T},
        values::Vector{T},
    ) where {T<:Real}

Represent a
[Lagrange interpolating polynomial](https://wikipedia.org/wiki/Lagrange_polynomial) (in the
barycentric form).

It is a subtype of [`AbstractInterpolant`](@ref AbstractInterpolant).
"""
struct LagrangeInterpolant <: AbstractInterpolant
    initial::Float64
    final::Float64
    points::Vector{Float64}
    weights::Vector{Float64}
    values::Vector{Float64}

    function LagrangeInterpolant(
        initial::T,
        final::T,
        points::Vector{T},
        weights::Vector{T},
        values::Vector{T},
    ) where {T<:Real}
        if !(length(points) == length(weights) == length(values))
            throw(ArgumentError("Please ensure length(points) == length(weights) == length(values)."))
        elseif !(issorted(points) == true)
            throw(ArgumentError("Please ensure issorted(points) == true."))
        elseif !(initial < first(points) || initial ≈ first(points))
            throw(DomainError(points, "Please ensure initial .≤ points."))
        elseif !(last(points) < final || last(points) ≈ final)
            throw(DomainError(points, "Please ensure points .≤ final."))
        end

        return new(initial, final, points, weights, values)
    end
end

"""
    (::LagrangeInterpolant)(t::Real)

Evaluate a Lagrange interpolant at `t`.
"""
function (interpolant::LagrangeInterpolant)(t::Real)
    if !(interpolant.initial ≤ t ≤ interpolant.final)
        throw(DomainError(t, "Ensure that initial ≤ t ≤ final."))
    end
    return _interpolate(interpolant.points, interpolant.weights, interpolant.values, t)
end

function _barycentric_weights(points::Vector{T}) where {T<:Real}

    if length(points) == 1
        return [one(T)]
    else
        return 1 ./ [
            prod(points[j] - points[k] for k in eachindex(points) if k != j)
            for j in eachindex(points)
        ]
    end
end

function _interpolate(
    points::Vector{T}, weights::Vector{T}, values::Vector{T}, point::T) where {T<:Real}

    numerator = zero(T)
    denominator = zero(T)

    for j in eachindex(points, weights, values)
        if points[j] == point
            return values[j]
        end

        common_j = weights[j] / (point - points[j])
        numerator += common_j * values[j]
        denominator += common_j
    end

    return numerator / denominator
end

function _differentiation_matrix(points::Vector{T}, weights::Vector{T}) where {T<:Real}

    D = Matrix{T}(undef, length(points), length(points))

    for j in eachindex(points, weights)
        D_jj = zero(T)
        for k in eachindex(points)
            if k != j
                D_jk = (weights[k] / weights[j]) / (points[j] - points[k])
                D[j,k] = D_jk
                D_jj -= D_jk
            end
        end
        D[j,j] = D_jj
    end
    return D
end

function _interpolation_matrix(
    points::Vector{T}, weights::Vector{T}, points_interp::Vector{T}) where {T<:Real}

    matrix = Matrix{T}(undef, length(points_interp), length(points))

    for k in eachindex(points_interp)

        exact = 0
        sum = zero(T)

        for j in eachindex(points, weights)
            point_difference = points_interp[k] - points[j]
            exact = point_difference == 0 ? j : exact
            matrix[k,j] = weights[j] / point_difference
            sum += matrix[k,j]
        end

        if sum == 0
            for j in eachindex(points)
                matrix[k,j] = zero(T)
            end
        elseif exact > 0
            for j in eachindex(points)
                matrix[k,j] = j == exact ? one(T) : zero(T)
            end
        else
            for j in eachindex(points)
                matrix[k,j] /= sum
            end
        end
    end
    return matrix
end

struct ZOHInterpolant <: Interesso.AbstractInterpolant
    initial::Float64
    final::Float64
    points::Vector{Float64}
    values::Vector{Float64}
end

function ZOHInterpolant(points::Vector{Float64}, values::Vector{Float64})
    length(points) == length(values) || error("length mismatch")
    return ZOHInterpolant(points[1], points[end], points, values)
end

function (itp::ZOHInterpolant)(s::Real)
    i = searchsortedlast(itp.points, s)
    i = clamp(i, 1, length(itp.values))
    return itp.values[i]
end

struct CubicInterpolant <: AbstractInterpolant
    initial::Float64
    final::Float64
    points::Vector{Float64}
    a::Vector{Float64}
    b::Vector{Float64}
    c::Vector{Float64}
    d::Vector{Float64}
end

function CubicInterpolant(points::Vector{Float64}, values::Vector{Float64})
    n = length(points) - 1
    h = diff(points)

    # Natural cubic spline (tridiagonal solve)
    α = zeros(n + 1)
    for i in 2:n
        α[i] = 3.0 / h[i] * (values[i+1] - values[i]) -
               3.0 / h[i-1] * (values[i] - values[i-1])
    end

    cvec = zeros(n + 1)
    l = ones(n + 1)
    μ = zeros(n + 1)
    z = zeros(n + 1)
    for i in 2:n
        l[i] = 2.0 * (points[i+1] - points[i-1]) - h[i-1] * μ[i-1]
        μ[i] = h[i] / l[i]
        z[i] = (α[i] - h[i-1] * z[i-1]) / l[i]
    end

    bvec = zeros(n)
    dvec = zeros(n)
    for j in n:-1:1
        cvec[j] = z[j] - μ[j] * cvec[j+1]
        bvec[j] = (values[j+1] - values[j]) / h[j] -
                  h[j] * (cvec[j+1] + 2.0 * cvec[j]) / 3.0
        dvec[j] = (cvec[j+1] - cvec[j]) / (3.0 * h[j])
    end

    return CubicInterpolant(points[1], points[end], points, values[1:n], bvec, cvec[1:n], dvec)
end

function (pc::CubicInterpolant)(s::Real)
    i = searchsortedlast(pc.points, s)
    i = clamp(i, 1, length(pc.a))
    ds = s - pc.points[i]
    return pc.a[i] + ds * (pc.b[i] + ds * (pc.c[i] + ds * pc.d[i]))
end
