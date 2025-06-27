const VAR = MOI.VariableIndex
const PHS = DOI.PhaseIndex
const DYN_VAR = DOI.DynamicVariableIndex
const NDF = DOI.NonlinearDynamicFunction
const NBF = DOI.NonlinearBoundaryFunction

const EQ64 = MOI.EqualTo{Float64}
const IV64 = MOI.Interval{Float64}

const STARTS = OrderedDict{DYN_VAR,DOI.AbstractDynamicSolution}

const DIF_FUN = DOI.ExplicitDifferentialFunction{NDF}

const DIF_CONS = OrderedDict{
    MOI.ConstraintIndex{DIF_FUN,EQ64},
    Tuple{DIF_FUN,EQ64},
}

const ALG_CONS = OrderedDict{
    MOI.ConstraintIndex{NDF,EQ64},
    Tuple{NDF,EQ64},
}

const LINKAGES = OrderedDict{
    MOI.ConstraintIndex{DOI.Linkage{DYN_VAR},EQ64},
    Tuple{DOI.Linkage{DYN_VAR},EQ64}
}

const OBJ = DOI.Bolza{NBF,DOI.MultiPhaseIntegral{NDF}}

const SOLS{F<:DOI.AbstractDynamicFunction} = OrderedDict{
    F,
    PiecewiseInterpolant{LagrangeInterpolant},
}

const MESHES = OrderedDict{PHS,AbstractIntervalsMesh}

const PHS_VARS = OrderedDict{PHS,Vector{VAR}}
const DYN_VAR_VARS = OrderedDict{DYN_VAR,Vector{Vector{VAR}}}