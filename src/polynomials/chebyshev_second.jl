struct Chebyshev1{T} <: Polynomials{T}
    degree::Integer
    nodes::Vector{T}
    weights::Vector{T}
    Chebyshev1{T}(degree::Integer) where {T<:Real} = new(degree,
        [cospi((2 * d + 1) / (2 * degree + 2)) for d in degree:-1:0],
        [((-1)^d) * sinpi((2 * d + 1) / (2 * degree + 2))  for d in degree:-1:0]
    )
end

Chebyshev1(::Type{T}, degree::Integer) where {T<:Real} = Chebyshev1{T}(degree);

Chebyshev1(degree::Integer) = Chebyshev1(Float64, degree);

struct Chebyshev2{T} <: Polynomials{T}
    degree::Integer
    nodes::Vector{T}
    weights::Vector{T}
    Chebyshev2{T}(degree::Integer) where {T<:Real} = new(degree,
        [cospi(d / degree) for d in degree:-1:0],
        [((-1)^d) * (1.0 - 0.5 * (d == degree || d == 0)) for d in degree:-1:0]
    )
end

Chebyshev2(::Type{T}, degree::Integer) where {T<:Real} = Chebyshev2{T}(degree);

Chebyshev2(degree::Integer) = Chebyshev2(Float64, degree);