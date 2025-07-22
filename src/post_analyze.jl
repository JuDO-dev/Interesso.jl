function eval_funcs(optimizer::MOI.ModelLike, funcs::Vector{<:MOI.AbstractFunction})
    var_idxs = MOI.get(optimizer, MOI.ListOfVariableIndices())
    x_vals   = MOI.get(optimizer, MOI.VariablePrimal(), var_idxs)

    nl = MOI.Nonlinear.Model()
    backend = MOI.Nonlinear.SparseReverseMode()

    for f in funcs
        MOI.Nonlinear.add_constraint(nl, f, MOI.EqualTo(0.0))
    end
    
    evaluator = MOI.Nonlinear.Evaluator(nl, backend, var_idxs)
    MOI.initialize(evaluator, Symbol[])
    eval = fill(NaN, length(funcs))
    MOI.eval_constraint(evaluator, eval, x_vals)

    return eval
end