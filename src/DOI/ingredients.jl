# MOI Variables

MOI.add_variable(model::Optimizer) = MOI.add_variable(model.inner)

MOI.is_valid(model::Optimizer, var::VAR) = MOI.is_valid(model.inner, var)

MOI.supports(::Optimizer, ::MOI.VariablePrimalStart, ::Type{VAR}) = true
function MOI.set(model::Optimizer, attr::MOI.VariablePrimalStart, var::VAR, value::Real) 
    return MOI.set(model.inner, attr, var, value)
end

# MOI Constraints

function MOI.supports_constraint(
    model::Optimizer,
    F::Type{<:MOI.AbstractFunction},
    S::Type{<:MOI.AbstractSet},
)
    return MOI.supports_constraint(model.inner, F, S)
end

function MOI.add_constraint(model::Optimizer, fun::MOI.AbstractFunction, set::MOI.AbstractSet)
    return MOI.add_constraint(model.inner, fun, set)
end

function MOI.is_valid(model::Optimizer, con::MOI.ConstraintIndex)
    return MOI.is_valid(model.inner, con)
end


# Boundary Conditions

MOI.supports_constraint(
    ::Optimizer,
    ::Type{<:DOI.AbstractBoundaryFunction},
    ::Type{<:LC64},
) = true

function MOI.add_constraint(
    model::Optimizer,
    fun::BF,
    set::S,
) where {BF<:DOI.AbstractBoundaryFunction,S<:LC64}

    index = MOI.ConstraintIndex{BF,S}(model.last_index_bou_cons + 1)
    model.bou_cons[index] = (fun, set)
    model.last_index_bou_cons += 1

    return index
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintFunction,
    con::MOI.ConstraintIndex{BF,S},
) where {BF<:DOI.AbstractBoundaryFunction,S<:LC64}
    return model.bou_cons[con][1]
end

function MOI.get(
    model::Optimizer,
    ::MOI.ConstraintSet,
    con::MOI.ConstraintIndex{BF,S},
) where {BF<:DOI.AbstractBoundaryFunction,S<:LC64}
    return model.bou_cons[con][2]
end


# Phases

DOI.supports_phases(::Optimizer) = true

function DOI.add_phase(model::Optimizer)

    phase = PHS(model.last_index_phases + 1)
    push!(model.phases, phase)

    model.dyn_vars[phase]         = OrderedSet{DYN_VAR}()
    model.dyn_var_bounds[phase]   = OrderedDict{DYN_VAR,LC64}()
    model.dyn_var_initials[phase] = OrderedDict{DYN_VAR,LC64}()
    model.dyn_var_finals[phase]   = OrderedDict{DYN_VAR,LC64}()
    model.dif_cons[phase]         = DIF_CONS()
    model.alg_cons[phase]         = ALG_CONS()
    model.path_cons[phase]        = PATH_CONS()
    model.start_dyn_vars[phase]   = STARTS()
    model.sol_dyn_vars[phase]     = SOLS{DYN_VAR}()

    model.last_index_phases += 1
    
    return phase
end

MOI.is_valid(model::Optimizer, phase::PHS) = phase in model.phases

function _throw_if_invalid_index(model::Optimizer, phase::PHS)
    if !MOI.is_valid(model, phase)
        throw(DOI.InvalidPhaseIndex(phase))
    end
    return nothing
end


# Phase boundaries

MOI.supports_constraint(::Optimizer, ::Type{DOI.Initial{PHS}}, ::Type{<:LC64}) = true
MOI.supports_constraint(::Optimizer, ::Type{DOI.Final{PHS}},   ::Type{<:LC64}) = true

function MOI.add_constraint(model::Optimizer, phase_initial::DOI.Initial{PHS}, set::S) where {S<:LC64}

    phase = phase_initial.dyn_fun

    _throw_if_invalid_index(model, phase)

    if haskey(model.phase_initials, phase)
        throw(MOI.AddConstraintNotAllowed{typeof(phase_initial),typeof(set)}(
            "Initial phase value already set."
        ))
    end

    model.phase_initials[phase] = set

    return MOI.ConstraintIndex{DOI.Initial{PHS},S}(phase.value)
end

function MOI.add_constraint(model::Optimizer, phase_final::DOI.Final{PHS}, set::S) where {S<:LC64}

    phase = phase_final.dyn_fun

    _throw_if_invalid_index(model, phase)

    if haskey(model.phase_finals, phase)
        throw(MOI.AddConstraintNotAllowed{typeof(phase_final),typeof(set)}(
            "Final phase value already set."
        ))
    end

    if set isa MOI.LessThan
        if set.upper <= 0.0
            throw(DomainError("Final time should have a positive upper bound."))
        end
    end

    if set isa MOI.Interval
        if set.lower <= 0.0
            throw(DomainError("Final time should have a positive lower bound."))
        end
    end

    model.phase_finals[phase] = set

    return MOI.ConstraintIndex{DOI.Final{PHS},S}(phase.value)
end


# Dynamic variables

DOI.supports_dynamic_variables(::Optimizer) = true

function DOI.add_dynamic_variable(model::Optimizer, phase::PHS)

    _throw_if_invalid_index(model, phase)

    dyn_var = DYN_VAR(model.last_index_dyn_vars + 1, phase)
    push!(model.dyn_vars[phase], dyn_var)
    model.last_index_dyn_vars += 1
    return dyn_var
end

function MOI.is_valid(model::Optimizer, dyn_var::DYN_VAR)
    phase = DOI.phase_index(dyn_var)
    return dyn_var in model.dyn_vars[phase]
end

function _throw_if_invalid_index(model::Optimizer, dyn_var::DYN_VAR)
    if !MOI.is_valid(model, dyn_var)
        throw(DOI.InvalidDynamicVariableIndex(dyn_var))
    end
    return nothing
end


# Dynamic variable bounds

MOI.supports_constraint(::Optimizer, ::Type{DYN_VAR}, ::Type{<:LC64}) = true

function MOI.add_constraint(model::Optimizer, dyn_var::DYN_VAR, set::S) where {S<:LC64}

    _throw_if_invalid_index(model, dyn_var)

    phase = DOI.phase_index(dyn_var)
    
    if haskey(model.dyn_var_bounds[phase], dyn_var)
        throw(MOI.AddConstraintNotAllowed{DYN_VAR,typeof(set)}("Bound already set."))
    end

    model.dyn_var_bounds[phase][dyn_var] = set

    return MOI.ConstraintIndex{DYN_VAR,typeof(set)}(dyn_var.value)
end


# Dynamic variable boundaries

MOI.supports_constraint(::Optimizer, ::DOI.Initial{DYN_VAR}, ::Type{<:LC64}) = true
MOI.supports_constraint(::Optimizer, ::DOI.Final{DYN_VAR},   ::Type{<:LC64}) = true

function MOI.add_constraint(model::Optimizer, dyn_var_initial::DOI.Initial{DYN_VAR}, set::S) where {S<:LC64}

    dyn_var = dyn_var_initial.dyn_fun

    _throw_if_invalid_index(model, dyn_var)
    
    phase = DOI.phase_index(dyn_var_initial.dyn_fun)
    
    if haskey(model.dyn_var_initials[phase], dyn_var)
        throw(MOI.AddConstraintNotAllowed{typeof(dyn_var_initial),S}("Initial value already set."))
    end

    model.dyn_var_initials[phase][dyn_var] = set

    return MOI.ConstraintIndex{DOI.Initial{DYN_VAR},S}(dyn_var.value)
end

function MOI.add_constraint(model::Optimizer, dyn_var_final::DOI.Final{DYN_VAR}, set::S) where {S<:LC64}

    dyn_var = dyn_var_final.dyn_fun

    _throw_if_invalid_index(model, dyn_var)
    
    phase = DOI.phase_index(dyn_var_final.dyn_fun)
    
    if haskey(model.dyn_var_finals[phase], dyn_var)
        throw(MOI.AddConstraintNotAllowed{typeof(dyn_var_final),S}("Final value already set."))
    end

    model.dyn_var_finals[phase][dyn_var] = set

    return MOI.ConstraintIndex{DOI.Final{DYN_VAR},S}(dyn_var.value)
end


# Dynamic variable names

MOI.supports(::Optimizer, ::DOI.DynamicVariableName) = true

function MOI.set(
    model::Optimizer,
    ::DOI.DynamicVariableName,
    dyn_var::DYN_VAR,
    name::String,
)
    _throw_if_invalid_index(model, dyn_var)

    model.dyn_var_names[dyn_var] = name

    return nothing
end

function MOI.get(
    model::Optimizer,
    ::DOI.DynamicVariableName,
    dyn_var::DYN_VAR,
)
    _throw_if_invalid_index(model, dyn_var)

    return get(model.dyn_var_names, dyn_var, nothing)
end


# Dynamic variable starts

MOI.supports(::Optimizer, ::DOI.DynamicVariableStart) = true

function MOI.set(
    model::Optimizer,
    ::DOI.DynamicVariableStart,
    dyn_var::DYN_VAR,
    start::DOI.AbstractDynamicSolution,
)
    _throw_if_invalid_index(model, dyn_var)

    phase = DOI.phase_index(dyn_var)

    model.start_dyn_vars[phase][dyn_var] = start

    return nothing
end

# Linkages

MOI.supports_constraint(::Optimizer, ::Type{DOI.Linkage{DYN_VAR}}, ::EQ64) = true

function MOI.add_constraint(model::Optimizer, linkage::DOI.Linkage{DYN_VAR}, set::EQ64)

    _throw_if_invalid_index(model, linkage.dyn_fun_final)
    _throw_if_invalid_index(model, linkage.dyn_fun_initial)

    index = MOI.ConstraintIndex{typeof(linkage),typeof(set)}(model.last_index_linkages)
    model.linkages[index] = (linkage, set)
    model.last_index_linkages += 1
    return nothing
end

# Differential Constraints

function MOI.supports_constraint(
    ::Optimizer,
    ::DIF_FUN,
    ::EQ64,
)
    return true
end

function MOI.add_constraint(
    model::Optimizer,
    dif_fun::DIF_FUN,
    set::EQ64,
)
    return MOI.add_constraint(model, dif_fun, set, 1.0)
end

function MOI.add_constraint(
    model::Optimizer,
    dif_fun::DIF_FUN,
    set::EQ64,
    scaling::Real,
)
    phase = DOI.phase_index(dif_fun)
    _throw_if_invalid_index(model, phase)
    index = MOI.ConstraintIndex{DIF_FUN,EQ64}(model.last_index_dif_cons + 1)
    model.dif_cons[phase][index] = (dif_fun, set, Float64(scaling))
    model.last_index_dif_cons += 1

    if !(dif_fun.dyn_var in model.dif_dyn_vars)
        push!(model.dif_dyn_vars, dif_fun.dyn_var)
    end

    return index
end


# Algebraic Equations

function MOI.supports_constraint(
    ::Optimizer,
    ::DOI.NonlinearDynamicFunction,
    ::Type{<:LC64},
)
    return true
end

function MOI.add_constraint(
    model::Optimizer,
    alg_fun::DOI.NonlinearDynamicFunction,
    set::EQ64,
)
    return MOI.add_constraint(model, alg_fun, set, 1.0)
end

function MOI.add_constraint(
    model::Optimizer,
    alg_fun::DOI.NonlinearDynamicFunction,
    set::EQ64,
    scaling::Real,
)
    phase = DOI.phase_index(alg_fun)
    _throw_if_invalid_index(model, phase)
    index = MOI.ConstraintIndex{DOI.NonlinearDynamicFunction,EQ64}(model.last_index_alg_cons + 1)
    model.alg_cons[phase][index] = (alg_fun, set, Float64(scaling))
    model.last_index_alg_cons += 1
    _push_dif_vars!(model, alg_fun)
    return index
end

function MOI.add_constraint(
    model::Optimizer,
    path_fun::DOI.NonlinearDynamicFunction,
    set::S,
) where {S<:Union{IV64,LE64,GE64}}
    phase = DOI.phase_index(path_fun)
    _throw_if_invalid_index(model, phase)
    index = MOI.ConstraintIndex{DOI.NonlinearDynamicFunction,S}(model.last_index_path_cons + 1)
    model.path_cons[phase][index] = (path_fun, set)
    model.last_index_path_cons += 1
    _push_dif_vars!(model, path_fun)
    return index
end

function _push_dif_vars!(::Optimizer, ::Any)
    return nothing
end

function _push_dif_vars!(model::Optimizer, fun::DOI.NonlinearDynamicFunction)
    for arg in fun.args
        _push_dif_vars!(model, arg)
    end
    return nothing
end

function _push_dif_vars!(model::Optimizer, fun::DOI.Derivative{DYN_VAR})
    if !(fun.dyn_fun in model.dif_dyn_vars)
        push!(model.dif_dyn_vars, fun.dyn_fun)
    end
    return nothing
end


# Objective

function MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{T}) where {T<:OBJ}
    return true
end

function MOI.set(
    model::Optimizer,
    ::MOI.ObjectiveFunction{T},
    obj_fun::T,
) where {T<:OBJ}
    model.objective = obj_fun
    return nothing
end

function MOI.get(model::Optimizer, ::MOI.ObjectiveFunction{T}) where {T<:OBJ}
    return model.objective
end