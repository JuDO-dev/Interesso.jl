function add_variables_and_constraints(m::MOI.ModelLike, problem, transcription, ctx::Context)

    x, u, v = add_variables!(m, ctx);

    set_start!(m, );

    add_bounds!(m, );

    add_continuity!(m, );

    return nothing
end
    
    
    
function add_variables!(m::MOI.ModelLike, ctx::Context)
    
    
    x = [[Vector{MOI.VariableIndex}(undef, ctx.n_τ_x) for _ in 1:ctx.n_x] for _ in 1:ctx.n_i];
    u = [[Vector{MOI.VariableIndex}(undef, ctx.n_τ_u) for _ in 1:ctx.n_u] for _ in 1:ctx.n_i];
    
    v = Vector{MOI.VariableIndex}(undef, ctx.n_v);
    
    return x, u, v
end

function set_start!(m, x, u, v)

    for i in 
        for j in 
            MOI.set(m, MOI.VariablePrimalStart(), x[i][j][k], start)

    return nothing
end
    
function add_bounds!(m::MOI.ModelLike, transcription)

    for i in 
        for j in 

        end
    end

    add
end

function add_continuity!(m)

    MOI.add_constraint(m, )
    
    
    
    
    
    
    
    
    
    
    
    z_variable_indices = MOI.add_variables(m, ctx.z_length);

    z_i_variable_indices = [Vector{MOI.VariableIndex}(undef, ctx.z_i_length) for _ in 1:ctx.n_i];

    for i in 1:ctx.n_i
        z_i_variable_indices[i] = view(z_variable_indices, z_indices);
    end


    v_variable_indices[j] = view()

    
    
    
    
    
    
    
    
    
    
    x, u, = add_variables()
    
    
    
    
    
    
    
    
    
    
    
    

    
    

        for j in 1:ctx.n_x
            for k in 1:ctx.n_τ_x
                x_variable_indices[i][j][k] = getindex_x(z_i_variable_indices, ctx, j, k)


            end

        end

        for j in 1:ctx.n_u
            for k in 1:ctx.n_τ_u
                



