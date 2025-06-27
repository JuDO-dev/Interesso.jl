module Interesso

import Ipopt
import MathOptInterface as MOI
import DynOptInterface as DOI
import FastGaussQuadrature as FGQ

using OrderedCollections: OrderedSet, OrderedDict

include("interpolants.jl")
include("points.jl")
include("methods.jl")
include("bounds.jl")
include("intervals.jl")

include("DOI/aliases.jl")
include("DOI/optimizer.jl")
include("DOI/attributes.jl")
include("DOI/ingredients.jl")
include("DOI/solutions.jl")

include("transcription/dyn_funs.jl")
include("transcription/bou_funs.jl")
include("transcription/ingredients.jl")
include("transcription/sol_dyn_var.jl")
include("transcription/sol_derivative.jl")

export AbstractInterpolant, PiecewiseInterpolant, LagrangeInterpolant
export AbstractPoints, AbstractPointsMesh, LGRPoints
export AbstractMethod, AbstractMethodMesh, Collocation
export AbstractBounds, AbstractBoundsMesh, SampledBounds, BernsteinBounds
export AbstractIntervals, AbstractIntervalsMesh, FixedIntervals, FlexibleIntervals

end