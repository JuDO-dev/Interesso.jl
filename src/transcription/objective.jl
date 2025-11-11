function transcribe_objective!(model::Optimizer, meshes::MESHES)

    MOI.set(
        model.inner,
        MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction}(),
        transcribe_objective(model.objective, model, meshes),
    )
    
    return nothing
end

transcribe_objective(::Nothing, ::Optimizer, ::MESHES) = nothing

function transcribe_objective(obj::OBJ, model::Optimizer, meshes::MESHES)
    terms = Any[]
    push!(terms, transcribe_bou_fun(obj, model, meshes))
    push!(terms, transcribe_penalty(model, meshes))
    return _add_terms(terms)
end

function transcribe_penalty(model::Optimizer, meshes::MESHES)
    terms = Any[]
    for (phase, mesh) in meshes
        push!(terms, transcribe_penalty(model, phase, mesh))
    end
    return _add_terms(terms)
end

function transcribe_penalty(
    ::Optimizer,
    ::PHS,
    ::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:CollocationMesh,BM}
    return nothing
end

function transcribe_penalty(
    ::Optimizer,
    ::PHS,
    ::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:DAIRMesh,BM}
    return nothing
end

function transcribe_penalty(
    ::Optimizer,
    ::PHS,
    ::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:ASIRMesh,BM}
    return nothing
end

function transcribe_penalty(
    model::Optimizer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:QPMMesh,BM}

    pen_fun = transcribe_dyn_least_square(model, phase, mesh)

    return MOI.ScalarNonlinearFunction(:*, [1.0, pen_fun])
end

function _add_terms(terms::Vector{Any})
    filtered = Any[]
    for term in terms
        if !_is_trivial_term(term)
            push!(filtered, term)
        end
    end

    if isempty(filtered)
        return nothing
    elseif length(filtered) == 1 && filtered[1] isa MOI.ScalarNonlinearFunction
        return filtered[1]
    else
        return MOI.ScalarNonlinearFunction(:+, filtered)
    end
end

_is_trivial_term(::Nothing) = true
_is_trivial_term(term::Number) = iszero(term)
_is_trivial_term(::Any) = false