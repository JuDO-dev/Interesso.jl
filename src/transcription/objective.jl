function transcribe_objective!(model::Optimizer, meshes::MESHES)
    terms = MOI.ScalarNonlinearFunction[]

    penalty_only = all(
        mesh -> mesh isa AbstractIntervalsMesh{<:Any,<:DAIRFeasMesh,<:Any},
        values(meshes),
    )

    if model.objective !== nothing && !penalty_only
        push!(terms, transcribe_bou_fun(model.objective, model, meshes))
    end

    for (phase, mesh) in meshes
        push!(terms, transcribe_penalty_terms(model, phase, mesh))
    end

    MOI.set(model.inner, MOI.ObjectiveSense(), model.objective_sense)
    MOI.set(
        model.inner,
        MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction}(),
        MOI.ScalarNonlinearFunction(:+, terms),
    )

    return nothing
end

function transcribe_penalty_terms(
    ::Optimizer,
    ::PHS,
    ::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:Union{CollocationMesh,DAIROptiMesh,SAIRMesh,GalerkinMesh},BM}
    return MOI.ScalarNonlinearFunction(:+, [])
end

function transcribe_penalty_terms(
    model::Optimizer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:DAIRFeasMesh,BM}
    return transcribe_dyn_least_square(model, phase, mesh)
end

function transcribe_penalty_terms(
    model::Optimizer,
    phase::PHS,
    mesh::AbstractIntervalsMesh{PM,MM,BM},
) where {PM,MM<:Union{QPMMesh,SAPMMesh},BM}
    method_mesh = get_method_mesh(mesh, 1)
    pen_fun = transcribe_dyn_least_square(model, phase, mesh)
    return MOI.ScalarNonlinearFunction(:*, [method_mesh.penalty, pen_fun])
end
