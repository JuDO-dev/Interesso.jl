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

include("post_solve/post_analyze.jl")
include("post_solve/mesh_refinement.jl")
include("post_solve/perturb.jl")
include("post_solve/plot.jl")

export AbstractInterpolant, PiecewiseInterpolant, LagrangeInterpolant, ZOHInterpolant, CubicInterpolant
export AbstractPoints, AbstractPointsMesh, LGRPoints, LGLPoints
export AbstractMethod, AbstractMethodMesh, AbstractIntRes, AbstractIntResMesh, Collocation, Galerkin, AbstractDAIR, DAIR, DAIRFeas, DAIROpti, QPM, SAIR, SAPM
export AbstractBounds, AbstractBoundsMesh, ExactBounds, SampledBounds, BernsteinBounds
export AbstractIntervals, AbstractIntervalsMesh, FixedIntervals, FlexibleIntervals
export @variable, @control
export get_solutions, warmstart!
export get_primal, get_dual, get_primal_dual, set_primal_start!, set_dual_start!
export eval_funcs, eval_accuracy, assess_solution, residual_map, summarize_residual
export refine!
export perturb_solution, perturb_solutions
export plot, plot_residual, plot_residual!

end