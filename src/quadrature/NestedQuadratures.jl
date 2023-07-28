module NestedQuadratures

using QuadGK: kronrod

abstract type NestedQuadrature{F<:AbstractFloat} end
export NestedQuadrature

# Gauss-Kronrod
include("GaussKronrod.jl")
export GaussKronrod

end