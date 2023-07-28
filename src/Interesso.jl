module Interesso

import MathOptInterface    as MOI
import FastGaussQuadrature as FGQ
import LinearAlgebra       as LA
import ReverseDiff         as RD
import Reexport
Reexport.@reexport using JuDOInterface

abstract type Intervals{T} end
abstract type Polynomials{T} end
abstract type Bounds end
abstract type Transcription{T, I<:Intervals, P<:Polynomials, B<:Bounds} end
abstract type Refinement{T, TT<:Transcription} end

include("intervals/rigid.jl")
export RigidIntervals

include("polynomials/legendre_lobatto.jl")
include("polynomials/barycentric_interpolation.jl")
export LegendreLobatto

include("bounds/simple.jl")
export SimpleApproximation

include("transcriptions/direct_collocation.jl")
#include("transcriptions/least_squares.jl")
export DirectCollocation#, LeastSquares

end