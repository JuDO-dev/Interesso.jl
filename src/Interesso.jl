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
include("DOI/macro.jl")

include("transcription/dyn_funs.jl")
include("transcription/bou_funs.jl")
include("transcription/objective.jl")
include("transcription/ingredients.jl")
include("transcription/sol_dyn_var.jl")
include("transcription/sol_derivative.jl")

include("post_solve/post_analyze.jl")
include("post_solve/perturb.jl")

export AbstractInterpolant, PiecewiseInterpolant, LagrangeInterpolant
export AbstractPoints, AbstractPointsMesh, LGRPoints, LGLPoints
export AbstractMethod, AbstractMethodMesh, AbstractIntRes, AbstractIntResMesh, Collocation, DAIR, QPM, SAIR, SAPM
export AbstractBounds, AbstractBoundsMesh, ExactBounds, SampledBounds, BernsteinBounds
export AbstractIntervals, AbstractIntervalsMesh, FixedIntervals, FlexibleIntervals
export get_solutions, warmstart!
export eval_funcs, eval_accuracy, assess_solution
export perturb_solution, perturb_solutions
export @variable

end