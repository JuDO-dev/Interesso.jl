# MOI Optimizer Attributes

MOI.get(::Optimizer, ::MOI.SolverName) = "Interesso"

MOI.get(::Optimizer, ::MOI.SolverVersion) = "0.0.0"

## Name
MOI.supports(::Optimizer, ::MOI.Name) = true
function MOI.set(model::Optimizer, ::MOI.Name, name::String) 
    model.name = name
    return nothing
end
MOI.get(model::Optimizer, ::MOI.Name) = model.name

## Silent
MOI.supports(model::Optimizer, attr::MOI.Silent) = MOI.supports(model.inner, attr)
MOI.set(model::Optimizer, attr::MOI.Silent, bool::Bool) = MOI.set(model.inner, attr, bool)
MOI.get(model::Optimizer, attr::MOI.Silent)             = MOI.get(model.inner, attr)

## ObjectiveSense
function MOI.supports(model::Optimizer, attr::MOI.ObjectiveSense)
    return MOI.supports(model.inner, attr)
end

function MOI.set(model::Optimizer, attr::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    MOI.set(model.inner, attr, sense)
    return nothing
end

MOI.get(model::Optimizer, attr::MOI.ObjectiveSense) = MOI.get(model.inner, attr)

MOI.supports_incremental_interface(::Optimizer) = true


## Interesso Optimizer Attributes

struct DefaultIntervals <: MOI.AbstractOptimizerAttribute end

function MOI.set(model::Optimizer, ::DefaultIntervals, intervals::AbstractIntervals)
    model.default_intervals = intervals
    return nothing
end

MOI.get(model::Optimizer, ::DefaultIntervals) = model.intervals


struct DefaultMethod <: MOI.AbstractOptimizerAttribute end

function MOI.set(model::Optimizer, ::DefaultMethod, method::AbstractMethod)
    model.default_method = method
    return nothing
end

MOI.get(model::Optimizer, ::DefaultMethod) = model.method


struct DefaultBounds <: MOI.AbstractOptimizerAttribute end

function MOI.set(model::Optimizer, ::DefaultBounds, bounds::AbstractBounds)
    model.default_bounds = bounds
    return nothing
end

MOI.get(model::Optimizer, ::DefaultBounds) = model.bounds


## Interesso Phase Attributes

struct PhaseIntervals <: DOI.AbstractPhaseAttribute end

function MOI.set(
    model::Optimizer,
    ::PhaseIntervals,
    phase::DOI.PhaseIndex,
    intervals::AbstractIntervals,
)
    _throw_if_invalid_index(model, phase)

    model.phase_intervals[phase] = intervals
    return nothing
end

struct PhaseMethod <: DOI.AbstractPhaseAttribute end

function MOI.set(
    model::Optimizer,
    ::PhaseMethod,
    phase::DOI.PhaseIndex,
    method::AbstractMethod,
)
    _throw_if_invalid_index(model, phase)

    model.phase_method[phase] = method
    return nothing
end

struct PhaseBounds <: DOI.AbstractPhaseAttribute end

function MOI.set(
    model::Optimizer,
    ::PhaseBounds,
    phase::DOI.PhaseIndex,
    bounds::AbstractBounds,
)
    _throw_if_invalid_index(model, phase)

    model.phase_bounds[phase] = bounds
    return nothing
end