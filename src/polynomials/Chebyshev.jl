struct Chebyshev1{F} <: BarycentricPolynomial{F}
    degree::Integer
    nodes::Vector{F}
    weights::Vector{F}
    Chebyshev1{F}(degree::Integer) where {F<:AbstractFloat} = new(degree,
        [cospi((2 * d + 1) / (2 * degree + 2)) for d in degree:-1:0],
        [((-1)^d) * sinpi((2 * d + 1) / (2 * degree + 2))  for d in degree:-1:0]
    )
end

Chebyshev1(::Type{F}, degree::Integer) where {F<:AbstractFloat} = Chebyshev1{F}(degree);

Chebyshev1(degree::Integer) = Chebyshev1(Float64, degree);

struct Chebyshev2{F} <: BarycentricPolynomial{F}
    degree::Integer
    nodes::Vector{F}
    weights::Vector{F}
    Chebyshev2{F}(degree::Integer) where {F<:AbstractFloat} = new(degree,
        [cospi(d / degree) for d in degree:-1:0],
        [((-1)^d) * (1.0 - 0.5 * (d == degree || d == 0)) for d in degree:-1:0]
    )
end

Chebyshev2(::Type{F}, degree::Integer) where {F<:AbstractFloat} = Chebyshev2{F}(degree);

Chebyshev2(degree::Integer) = Chebyshev2(Float64, degree);