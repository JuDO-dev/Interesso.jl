const VAR = MOI.VariableIndex
const PHS = DOI.PhaseIndex
const TIME_VAR = Union{Float64,VAR}
const DYN_VAR = DOI.DynamicVariableIndex
const NDF = DOI.NonlinearDynamicFunction
const NBF = DOI.NonlinearBoundaryFunction

const EQ64 = MOI.EqualTo{Float64}
const IV64 = MOI.Interval{Float64}
const LE64 = MOI.LessThan{Float64}
const GE64 = MOI.GreaterThan{Float64}
const LC64 = Union{EQ64,IV64,LE64,GE64}

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

const PATH_CONS = OrderedDict{
    Union{
        MOI.ConstraintIndex{NDF,IV64},
        MOI.ConstraintIndex{NDF,LE64},
        MOI.ConstraintIndex{NDF,GE64}
    },
    Tuple{NDF,Union{IV64,LE64,GE64}},
}

const BOU_CONS = OrderedDict{
    Union{
        MOI.ConstraintIndex{NBF,EQ64},
        MOI.ConstraintIndex{NBF,IV64},
        MOI.ConstraintIndex{NBF,LE64},
        MOI.ConstraintIndex{NBF,GE64}
    },
    Tuple{NBF,LC64},
}

const LINKAGES = OrderedDict{
    MOI.ConstraintIndex{DOI.Linkage{DYN_VAR},EQ64},
    Tuple{DOI.Linkage{DYN_VAR},EQ64}
}

const BOLZA = DOI.Bolza{NBF,DOI.MultiPhaseIntegral{NDF}}
const OBJ = Union{NBF, DOI.MultiPhaseIntegral{NDF}, BOLZA}

const SOLS{F<:DOI.AbstractDynamicFunction} = OrderedDict{
    F,
    PiecewiseInterpolant{LagrangeInterpolant},
}

const MESHES = OrderedDict{PHS,AbstractIntervalsMesh}

const PHS_VARS = OrderedDict{PHS,Vector{VAR}}
const TIME_VARS = OrderedDict{PHS,TIME_VAR}
const DYN_VAR_VARS = OrderedDict{DYN_VAR,Vector{Vector{VAR}}}

const WS{S<:DOI.AbstractDynamicSolution} = OrderedDict{String,S}
const WSS{S<:DOI.AbstractDynamicSolution} = OrderedDict{PHS,OrderedDict{String,S}}