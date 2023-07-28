function add_continuity!(m::MOI.ModelLike, x::Vector{Vector{Vector{VariableIndex}}}, ctx)
    for i in 2:ctx.n_i
        for j in 1:ctx.n_x
            MOI.add_constraint(m, MOI.ScalarAffineFunction(
                [MOI.ScalarAffineTerm(1.0, x[i][j][1]), MOI.ScalarAffineTerm(-1.0, x[i-1][j][end])]
            ), MOI.EqualTo(0.0));
        end
    end
    return nothing
end