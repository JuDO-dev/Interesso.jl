struct DirectCollocation{T, I, P, B} <: Transcription{T, I, P, B}
    intervals::I
    x_polynomials::P
    u_polynomials::P
    q_polynomials::P
    x_bounds::B
    u_bounds::B

    D::Matrix{T}
    
    function DirectCollocation(intervals::I, polynomials::P, bounds::B) where {T,
        I<:Intervals{T}, P<:Polynomials{T}, B<:Bounds}
        
        D = barycentric_differentiation(polynomials.τ, polynomials.w_b);
        
        return new{T, I, P, B}(intervals, polynomials, polynomials, polynomials, bounds, bounds, D)
    end
end

z_length(dop::DOProblem, transcription::Transcription) = transcription.intervals.n * (length(dop.x) * transcription.x_polynomials.n + length(dop.u) * transcription.u_polynomials.n);

mutable struct DirectCollocationEvaluator <: MOI.AbstractNLPEvaluator

    z_ranges
    r_ranges
    J_ranges
    J_structure

    x_indices
    u_indices
    
    g∫c_tape    
    Jr_tape

    function DirectCollocationEvaluator(z_indices::Vector{MOI.VariableIndex}, dop::DOProblem, dc::DirectCollocation)
 
        n_i = dc.intervals.n;
        n_x = length(dop.x);
        n_u = length(dop.u);
        n_r = length(dop.equations);
        n_τ_x = dc.x_polynomials.n;
        n_τ_u = dc.u_polynomials.n;
        n_τ_q = dc.q_polynomials.n;
        x_i_length = n_x * n_τ_x;
        u_i_length = n_u * n_τ_u;
        z_i_length = x_i_length + u_i_length;
        r_i_length = n_r * n_τ_q;
        J_i_length = z_i_length * r_i_length;
        
        z_ranges = [UnitRange((i-1) * z_i_length + 1, i * z_i_length) for i in 1:n_i];
        r_ranges = [UnitRange((i-1) * r_i_length + 1, i * r_i_length) for i in 1:n_i];
        J_ranges = [UnitRange((i-1) * J_i_length + 1, i * J_i_length) for i in 1:n_i];
        J_structure = [(
            (i-1) * r_i_length + m,
            (i-1) * z_i_length + l
        ) for i in 1:n_i for l in 1:z_i_length for m in 1:r_i_length];

        # Variable indices
        @views begin
            z_indices_views = [z_indices[range] for range in z_ranges];
            x_indices = [reshape(z_i[1:x_i_length], n_x, n_τ_x) for z_i in z_indices_views];
            u_indices = [reshape(z_i[x_i_length + 1:z_i_length], n_u, n_τ_u) for z_i in z_indices_views];
        end;
        
        # Cost
        #x_i_ranges = [UnitRange((k-1) * n_x + 1, k * n_x) for k in 1:n_τ_x];
        #u_i_ranges = [UnitRange((k-1) * n_u + 1, k * n_u) for k in 1:n_τ_u];

        function ∫c_closure(z_i::AbstractArray{<:Real})
            @views begin
                x_i = reshape(z_i[1:x_i_length], n_x, n_τ_x);
                u_i = reshape(z_i[x_i_length + 1:z_i_length], n_u, n_τ_u);
            end
            
            return sum(dc.q_polynomials.w_q[k] * dop.cost.running(x_i[:,k], u_i[:,k])
                for k in eachindex(dc.q_polynomials.τ)
            )
        end

        # Gradient of the cost
        g∫c_tape = RD.compile(RD.GradientTape(∫c_closure, rand(z_i_length)));

        # Residuals
        r! = get_residuals(dop.equations);
        function r_closure!(r_i::AbstractArray{<:Real}, z_i::AbstractArray{<:Real}) 
            @views begin
                r_i = reshape(r_i,                          n_r, n_τ_q);
                x_i = reshape(z_i[1:x_i_length], n_x, n_τ_x);
                u_i = reshape(z_i[x_i_length + 1:z_i_length], n_u, n_τ_u);
            end;
            
            ẋ_i = dc.D * permutedims(x_i, (2,1));
            
            for k in eachindex(dc.q_polynomials.τ)
                r!(r_i[:,k], ẋ_i[k,:], x_i[:,k], u_i[:,k]);
            end
            return nothing
        end

        # Jacobian of the residuals
        Jr_tape = RD.compile(RD.JacobianTape(r_closure!, rand(r_i_length), rand(z_i_length)));

        return new(z_ranges, r_ranges, J_ranges, J_structure, x_indices, u_indices, g∫c_tape, Jr_tape)
    end
end

function add_continuity!(m::MOI.ModelLike, variable_indices)

    (j_range, k_range) = axes(variable_indices[begin]);
    
    for i in 2:length(variable_indices)
        for j in j_range, k in k_range
            MOI.add_constraint(m, MOI.ScalarAffineFunction([
                MOI.ScalarAffineTerm( 1.0, variable_indices[i][j,begin]),
                MOI.ScalarAffineTerm(-1.0, variable_indices[i-1][j,end])
            ], 0.0), MOI.EqualTo(0.0));
        end
    end
    return nothing
end

function add_boundary_conditions!(m::MOI.ModelLike, variable_indices, variables)

    for j in eachindex(variables)
        MOI.add_constraint(m, MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, variable_indices[begin][j,begin])], 0.0), MOI.GreaterThan(variables[j].initial.lower));
        MOI.add_constraint(m, MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, variable_indices[begin][j,begin])], 0.0),    MOI.LessThan(variables[j].initial.upper));
        MOI.add_constraint(m, MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, variable_indices[end][j,end])],     0.0), MOI.GreaterThan(variables[j].final.lower));
        MOI.add_constraint(m, MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, variable_indices[end][j,end])],     0.0),    MOI.LessThan(variables[j].final.upper));
    end
    return nothing
end

function MOI.eval_objective(dce::DirectCollocationEvaluator, z::AbstractVector{T})::T where {T}

    @views begin
        z_views = [z[range] for range in dce.z_ranges];
    end;

    return sum(dce.g∫c_tape.tape.func(z_views[i]) for i in eachindex(z_views))
end

function MOI.eval_objective_gradient(dce::DirectCollocationEvaluator, gcz::AbstractVector{T}, z::AbstractVector{T})::Nothing where {T}

    @views begin
        gcz_views = [gcz[range] for range in dce.z_ranges];
          z_views =   [z[range] for range in dce.z_ranges];
    end;
    
    for i in eachindex(gcz_views, z_views)
        RD.gradient!(gcz_views[i], dce.g∫c_tape, z_views[i]);
    end
    return nothing
end

function MOI.eval_constraint(dce::DirectCollocationEvaluator, rz::AbstractVector{T}, z::AbstractVector{T})::Nothing where {T}

    @views begin
        rz_views = [rz[range] for range in dce.r_ranges];
         z_views =  [z[range] for range in dce.z_ranges];
    end;
    
    for i in eachindex(rz_views, z_views)
        dce.Jr_tape.tape.func(rz_views[i], z_views[i]);
    end
    return nothing
end

MOI.jacobian_structure(dce::DirectCollocationEvaluator)::Vector{Tuple{Int64,Int64}} = dce.J_structure;

function MOI.eval_constraint_jacobian(dce::DirectCollocationEvaluator, J::AbstractVector{T}, z::AbstractVector{T})::Nothing where {T}
    
    @views begin
        J_views = [J[range] for range in dce.J_ranges];
        z_views = [z[range] for range in dce.z_ranges];
    end;

    for i in eachindex(J_views, z_views)
        RD.jacobian!(J_views[i], dce.Jr_tape, z_views[i]);
    end
    return nothing
end

nlp_block_data(dce::DirectCollocationEvaluator) = MOI.NLPBlockData([MOI.NLPBoundsPair(0.0, 0.0) for _ in Base.OneTo(length(dce.J_structure))] , dce, true);

get_residuals(equations::ResidualEquations) = equations.r!;