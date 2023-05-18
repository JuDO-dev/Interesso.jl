module Interesso

import Reexport
Reexport.@reexport using JuDOInterface
import Progradio
import FastGaussQuadrature
import Enzyme

abstract type PointDistribution{T} end
include("point_distributions/legendre_lobatto.jl")
#include("point_distributions/chebyshev_second.jl")
export LegendreLobatto#, ChebyshevSecond

abstract type Intervals{T} end
include("intervals/rigid.jl")
#include("intervals/flexible.jl")
export RigidIntervals#, FlexibleIntervals

# Transcription
abstract type Transcription{T, I<:Intervals{T}} end
include("transcription/discretization.jl")
include("transcription/direct_collocation.jl")
include("transcription/least_squares.jl")
include("transcription/start.jl")
include("transcription/bounds.jl")
include("transcription/transcribe.jl")
export LeastSquares, DirectCollocation, transcribe

# Refinement
#include("refinement/partition.jl")
#abstract type AbstractRefinement{F<:AbstractFloat} end
#abstract type AbstractRefinementState{F<:AbstractFloat} end
#include("refinement/no.jl")
#include("refinement/convergent.jl")
#include("refinement/predictive.jl")
#export NoRefinement#, ConvergentRefinement, PredictiveRefinement

# Solve
#include("solve.jl")

include("barycentric.jl")

end