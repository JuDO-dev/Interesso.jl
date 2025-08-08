function transcribe_objective!(model::Optimizer, meshes::MESHES)

    if !(model.objective isa Nothing)

        MOI.set(
            model.inner,
            MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction}(),
            transcribe_bou_fun(model.objective, model, meshes),
        )
    end
    
    return nothing
end