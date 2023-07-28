module BarycentricPolynomials

abstract type BarycentricPolynomial{F<:AbstractFloat} end
export BarycentricPolynomial

# Chebyshev Polynomials
include("Chebyshev.jl")
export Chebyshev1, Chebyshev2

# Differentiation and Interpolation Matrices
include("matrices.jl")
export differentiation_matrix, interpolation_matrix

end