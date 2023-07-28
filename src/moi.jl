mutable struct NLPEvaluator{T, PT<:DOProblem, TT<:Transcription} <: MOI.AbstractNLPEvaluator
    problem::PT
    transcription::TT
    




    gradient_caches::Vector{Vector{Cache{T}}}
    gradient_adjoint_caches::Vector{Vector{Cache{T}}}

    jacobian_caches::Vector{Vector{Cache{T}}}
    jacobian_adjoint_caches::Vector{Vector{NTuple{N, Cache{T}}}}







    z_i_length::Int
    z_length::Int
    z_offsets::Vector{Int}
    x_offsets::Vector{Int}
    u_offsets::Vector{Int}
    
    x_idx::Vector{Vector{Vector{MOI.VariableIndex}}}
    u_idx::Vector{Vector{Vector{MOI.VariableIndex}}}
    #v::Vector{MOI.VariableIndex}
    #t::Vector{MOI.VariableIndex}

    caches::Vector{Cache{T}}
    adjoint_caches::Vector{Cache{T}}

    function NLPEvaluator(z_idx::Vector{MOI.VariableIndex}, problem::PT, transcription::TT; T::Type=Float64) where {PT<:DOProblem, TT<:Transcription}

        n_x = length(problem.x);
        n_u = length(problem.u);
        n_r = length(problem.equations);
        n_i = transcription.intervals.n;
        n_τ_x = transcription.x_polynomials.n;
        n_τ_u = transcription.u_polynomials.n;
        n_τ_q = transcription.q_polynomials.n;

        z_i_length = n_x * n_τ_x + n_u * n_τ_u;
        z_length = n_i * z_i_length;
        z_offsets = [(i-1) * z_i_length          for i in 1:n_i];
        x_offsets = [(j-1) * n_τ_x               for j in 1:n_x];
        u_offsets = [n_x * n_τ_x + (j-1) * n_τ_u for j in 1:n_u];
        
        x_idx = [[Vector{MOI.VariableIndex}(undef, n_τ_x) for _ in 1:n_x] for _ in 1:n_i];
        u_idx = [[Vector{MOI.VariableIndex}(undef, n_τ_u) for _ in 1:n_u] for _ in 1:n_i];
        for i in eachindex(x_idx, u_idx)
            for j in eachindex(x_idx[i])
                for k in eachindex(x_idx[i][j]) x_idx[i][j][k] = z_idx[z_offsets[i]+x_offsets[j]+k]; end
            end
            for j in eachindex(u_idx[i])
                for k in eachindex(u_idx[i][j]) u_idx[i][j][k] = z_idx[z_offsets[i]+u_offsets[j]+k]; end
            end
        end

        caches =         [Cache(z_i_length, n_τ_q, n_x, n_u, n_r) for _ in 1:n_i];
        adjoint_caches = [Cache(z_i_length, n_τ_q, n_x, n_u, n_r) for _ in 1:n_i];

        return new{T, PT, TT}(problem, transcription, z_i_length, z_length, z_offsets, x_offsets, u_offsets, x_idx, u_idx, caches, adjoint_caches)
    end
end

function add_continuity!(m::MOI.ModelLike, x_idx::Vector{Vector{Vector{MOI.VariableIndex}}})

    for i in 2:length(x_idx)
        for j in eachindex(x_idx[i])
            MOI.add_constraint(m, MOI.ScalarAffineFunction([
                MOI.ScalarAffineTerm(1.0, x_idx[i][j][1]),
                MOI.ScalarAffineTerm(-1.0, x_idx[i-1][j][end])
            ], 0.0), MOI.EqualTo(0.0));
        end
    end
    return nothing
end

function add_boundary_conditions!(m::MOI.ModelLike, x::Vector{DifferentialVariable{T, IT}}, x_idx::Vector{Vector{Vector{MOI.VariableIndex}}}) where {T, IT}

    for j in eachindex(x, x_idx[1], x_idx[end])
        MOI.add_constraint(m, MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, x_idx[1][j][1])],     0.0), MOI.GreaterThan(x[j].initial.lower));
        MOI.add_constraint(m, MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, x_idx[1][j][1])],     0.0),    MOI.LessThan(x[j].initial.upper));
        MOI.add_constraint(m, MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, x_idx[end][j][end])], 0.0), MOI.GreaterThan(x[j].final.lower));
        MOI.add_constraint(m, MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, x_idx[end][j][end])], 0.0),    MOI.LessThan(x[j].final.upper));
    end
    return nothing
end