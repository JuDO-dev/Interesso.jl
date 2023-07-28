using Test
using Interesso
const I = Interesso;

include("barycentric_interpolation.jl")

include("rules/legendre_lobatto.jl")
#include("rules/chebyshev_second.jl")

include("intervals/rigid.jl")
#include("intervals/flexible.jl")

#include("transcription/discretizations.jl")
#include("transcription/warmstart.jl")
#include("transcription/bounds.jl")