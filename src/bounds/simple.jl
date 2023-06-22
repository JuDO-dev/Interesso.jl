struct SimpleApproximation <: Bounds end

function add_bounds!(m::MOI.ModelLike, variable_indices, variables)

    for i in eachindex(variable_indices)
        for j in eachindex(variables)
            for k in axes(variable_indices, 2)
                MOI.add_constraint(m, variable_indices[i][j,k], MOI.GreaterThan(variables[j].bounds.lower));
                MOI.add_constraint(m, variable_indices[i][j,k],    MOI.LessThan(variables[j].bounds.upper));
            end
        end
    end
    return nothing
end

#=function add_bounds!(m::MOI.ModelLike, e::NLPEvaluator) 
    
    for i in eachindex(e.x_idx, e.u_idx)
        for j in eachindex(e.problem.x, e.x_idx[i])
            for k in eachindex(e.x_idx[i][j])
                MOI.add_constraint(m, e.x_idx[i][j][k], MOI.GreaterThan(e.problem.x[j].bounds.lower));
                MOI.add_constraint(m, e.x_idx[i][j][k],    MOI.LessThan(e.problem.x[j].bounds.upper));
            end
        end
        for j in eachindex(e.problem.u, e.u_idx[i])
            for k in eachindex(e.x_idx[i][j])
                MOI.add_constraint(m, e.u_idx[i][j][k], MOI.GreaterThan(e.problem.u[j].bounds.lower));
                MOI.add_constraint(m, e.u_idx[i][j][k],    MOI.LessThan(e.problem.u[j].bounds.upper));
            end
        end
    end
    return nothing
end=#