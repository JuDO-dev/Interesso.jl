module Interesso

using Reexport
@reexport using JuDOBase
@reexport using Progradio
#using AbstractDifferentiation
using Enzyme
#using Zygote
#using ReverseDiff
using QuadGK: gauss, kronrod
#using FastTransforms: clenshawcurtisnodes, chebyshevmoments1, clenshawcurtisweights
using LinearAlgebra: dot

# Barycentric Polynomials
include("polynomials/BarycentricPolynomials.jl")
using .BarycentricPolynomials
export Chebyshev1, Chebyshev2

# Nested Quadrature
include("quadrature/NestedQuadratures.jl")
using .NestedQuadratures
export GaussKronrod

# Mesh
abstract type InteressoMesh{F<:AbstractFloat} end
include("mesh/rigid.jl")
include("mesh/flexible.jl")
export RigidMesh, FlexibleMesh

# Transcription
abstract type InteressoTranscription{F<:AbstractFloat, M<:InteressoMesh} end
include("transcription/bounds.jl")
include("transcription/interpolation.jl")
#include("transcription/warmstarting.jl")
include("transcription/leastSquares.jl")
#include("transcription/collocation.jl")
export LeastSquares

# Refinement
#abstract type InteressoRefinement{F<:AbstractFloat, T<:InteressoTranscription} end
#include("refinement/no.jl")
#include("refinement/convergent.jl")
#include("refinement/predictive.jl")

# Solve
#include("solve.jl")

end