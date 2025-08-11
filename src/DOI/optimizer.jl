mutable struct Optimizer <: MOI.AbstractOptimizer

    # Attributes
    default_intervals::AbstractIntervals
    default_points::AbstractPoints
    default_method::AbstractMethod
    default_bounds::AbstractBounds
    
    # Dynamic Optimization Problem
    name::String
    phases::OrderedSet{PHS}
    phase_initials::OrderedDict{PHS,Float64}
    phase_finals::OrderedDict{PHS,Float64}
    dyn_vars::OrderedDict{PHS,OrderedSet{DYN_VAR}}
    dyn_var_bounds::OrderedDict{PHS,OrderedDict{DYN_VAR,IV64}}
    dyn_var_initials::OrderedDict{PHS,OrderedDict{DYN_VAR,EI64}}
    dyn_var_finals::OrderedDict{PHS,OrderedDict{DYN_VAR,EI64}}
    linkages::LINKAGES
    dif_dyn_vars::OrderedSet{DYN_VAR}
    dif_cons::OrderedDict{PHS,DIF_CONS}
    alg_cons::OrderedDict{PHS,ALG_CONS}
    objective_sense::MOI.OptimizationSense
    objective::Union{OBJ,Nothing}

    last_index_phases::Int64
    last_index_dyn_vars::Int64
    last_index_dif_cons::Int64
    last_index_alg_cons::Int64
    last_index_linkages::Int64

    # Start
    start_dyn_vars::OrderedDict{PHS,STARTS}

    # Phase Attributes
    phase_intervals::OrderedDict{PHS,<:AbstractIntervals}
    phase_points::OrderedDict{PHS,<:AbstractPoints}
    phase_method::OrderedDict{PHS,<:AbstractMethod}
    phase_bounds::OrderedDict{PHS,<:AbstractBounds}

    # Transcription
    meshes::MESHES
    inner::MOI.AbstractOptimizer
    phase_vars::PHS_VARS
    dyn_var_vars::DYN_VAR_VARS
    penalty_funs::OrderedDict{PHS,MOI.ScalarNonlinearFunction}

    # Solution
    sol_dyn_vars::OrderedDict{PHS,SOLS{DYN_VAR}}
    sol_derivatives::OrderedDict{PHS,SOLS{DOI.Derivative{DYN_VAR}}}

    # Analyze
    dif_res_funcs::Vector{MOI.AbstractFunction}
    res_funcs::Vector{MOI.AbstractFunction}

    function Optimizer(;
        inner::MOI.ModelLike=Ipopt.Optimizer(),
        default_intervals::AbstractIntervals=FixedIntervals(1),
        default_points::AbstractPoints=LGRPoints(5),
        default_method::AbstractMethod=Collocation(),
        default_bounds::AbstractBounds=ExactBounds(),
    )

        if default_method isa Collocation
            if !(length(default_points.points_dif_τ) == length(default_points.points_alg_τ) + 1)
                throw(DomainError("Collocation method requires states and control to be of same order."))
            end
            if default_bounds isa SampledBounds
                throw(DomainError("Collocation method does not support sampled bounds."))
            end
        end

        inner = MOI.Bridges.full_bridge_optimizer(inner, Float64)

        return new(
            default_intervals,
            default_points,
            default_method,
            default_bounds,
            "",
            OrderedSet{PHS}(),
            OrderedDict{PHS,Float64}(),
            OrderedDict{PHS,Float64}(),
            OrderedDict{PHS,OrderedSet{DYN_VAR}}(),
            OrderedDict{PHS,OrderedDict{DYN_VAR,IV64}}(),
            OrderedDict{PHS,OrderedDict{DYN_VAR,EI64}}(),
            OrderedDict{PHS,OrderedDict{DYN_VAR,EI64}}(),
            LINKAGES(),
            OrderedSet{DYN_VAR}(),
            OrderedDict{PHS,DIF_CONS}(),
            OrderedDict{PHS,ALG_CONS}(),
            MOI.FEASIBILITY_SENSE,
            nothing,
            0,
            0,
            0,
            0,           
            0,
            OrderedDict{PHS,STARTS}(),
            OrderedDict{PHS,AbstractIntervals}(),
            OrderedDict{PHS,AbstractPoints}(),
            OrderedDict{PHS,AbstractMethod}(),
            OrderedDict{PHS,AbstractBounds}(),
            MESHES(),
            inner,
            PHS_VARS(),
            DYN_VAR_VARS(),
            OrderedDict{PHS,MOI.ScalarNonlinearFunction}(),
            OrderedDict{PHS,SOLS{DYN_VAR}}(),
            OrderedDict{PHS,SOLS{DOI.Derivative{DYN_VAR}}}(),
            Vector{MOI.AbstractFunction}(),
            Vector{MOI.AbstractFunction}(),
        )
    end
end

function MOI.empty!(model::Optimizer)

    empty!(model.phases)
    empty!(model.phase_initials)
    empty!(model.phase_finals)
    empty!(model.dyn_vars)
    empty!(model.dyn_var_bounds)
    empty!(model.dyn_var_initials)
    empty!(model.dyn_var_finals)
    empty!(model.linkages)
    empty!(model.dif_dyn_vars)
    empty!(model.dif_cons)
    empty!(model.alg_cons)
    model.objective_sense = MOI.FEASIBILITY_SENSE
    model.objective = nothing
    model.last_index_phases = 0
    model.last_index_dyn_vars = 0
    model.last_index_dif_cons = 0
    model.last_index_alg_cons = 0
    model.last_index_linkages = 0
    empty!(model.start_dyn_vars)
    empty!(model.phase_intervals)
    empty!(model.phase_points)
    empty!(model.phase_method)
    empty!(model.phase_bounds)
    empty!(model.meshes)
    MOI.empty!(model.inner)
    empty!(model.phase_vars)
    empty!(model.dyn_var_vars)
    empty!(model.penalty_funs)
    empty!(model.sol_dyn_vars)
    empty!(model.sol_derivatives)

    return nothing
end

function MOI.is_empty(model::Optimizer)

    return isempty(model.phases)          &&                                     
        isempty(model.phase_initials)     && isempty(model.phase_finals)        &&
        isempty(model.dyn_vars)           && isempty(model.dyn_var_bounds)      &&
        isempty(model.dyn_var_initials)   && isempty(model.dyn_var_finals)      &&
        isempty(model.linkages)           && isempty(model.dif_dyn_vars)        &&
        isempty(model.dif_cons)           && isempty(model.alg_cons)            &&
        model.objective_sense == MOI.FEASIBILITY_SENSE                          &&
        isnothing(model.objective)        &&
        iszero(model.last_index_phases)   && iszero(model.last_index_dyn_vars)  &&
        iszero(model.last_index_dif_cons) && iszero(model.last_index_alg_cons)  &&
        iszero(model.last_index_linkages) && isempty(model.start_dyn_vars)      &&
        isempty(model.phase_intervals)    && isempty(model.phase_points)        &&
        isempty(model.phase_method)       && isempty(model.phase_bounds)        &&
        
        isempty(model.meshes)             && MOI.is_empty(model.inner)  && 
        isempty(model.phase_vars)         && isempty(model.dyn_var_vars)        &&
        isempty(model.penalty_funs)       &&
        isempty(model.sol_dyn_vars)       && isempty(model.sol_derivatives)
end

function MOI.optimize!(model::Optimizer)

    ## Build Mesh

    if isempty(model.meshes)

        for phase in model.phases
            
            if !haskey(model.phase_initials, phase)
                error("Please ensure that all phases have a fixed initial value.")
            elseif !haskey(model.phase_finals, phase)
                error("Please ensure that all phases have a fixed final value.")
            end
            
            model.meshes[phase] = build_intervals_mesh(
                get(model.phase_intervals, phase, model.default_intervals),
                get(model.phase_points, phase, model.default_points),
                get(model.phase_method, phase, model.default_method),
                get(model.phase_bounds, phase, model.default_bounds),
                model.phase_initials[phase],
                model.phase_finals[phase],
            )
        end
    end

    ## Transcribe Problem

    for phase in model.phases

        transcribe_phase!(model, phase, model.meshes[phase])

        transcribe_dyn_vars!(model, phase, model.meshes[phase])

        transcribe_dyn_var_starts!(model, phase, model.meshes[phase])
        
        transcribe_bounds!(model, phase, model.meshes[phase])
    
        transcribe_initials!(model, phase, model.meshes[phase])
    
        transcribe_finals!(model, phase, model.meshes[phase])
    
        transcribe_dif_cons!(model, phase, model.meshes[phase])

        transcribe_alg_cons!(model, phase, model.meshes[phase])
    end

    transcribe_linkages!(model, model.meshes)

    transcribe_objective!(model, model.meshes)

    ## Optimize
    MOI.optimize!(model.inner)

    ## Update Meshes
    for phase in model.phases

        update_mesh!(
            model.meshes[phase],
            model.inner,
            model.phase_vars,
            phase,
            get(model.phase_points, phase, model.default_points),
        )
    end
    return nothing
end