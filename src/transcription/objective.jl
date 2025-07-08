function transcribe_objective!(model::Optimizer, meshes::MESHES)
    
    if (model.objective isa Nothing)

        MOI.set(
            model.inner,
            MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction}(),
            transcribe_dyn_least_square(model, meshes),
        )

    else

        MOI.set(
            model.inner,
            MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction}(),
            MOI.ScalarNonlinearFunction(
                :+,
                [
                    transcribe_bou_fun(model.objective, model, meshes),
                    MOI.ScalarNonlinearFunction(
                        :*,
                        [
                            1.0,   # some penalty parameter, intended to be dynamic
                            transcribe_dyn_least_square(model, meshes),
                        ]
                    )
                ]
            )
        )

    end
    return nothing
end