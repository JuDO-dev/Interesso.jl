struct NoRefinement{P,T,O} <: Refinement{T}
    problem::P
    transcription::T
    optimizer::O
end

function Base.iterate(nr::NoRefinement)

    MOI.empty!(nr.optimizer);
    z_idx = MOI.add_variables(z_length(nr.problem, nr.transcription));
    nlp_evaluator = NLPEvaluator(z_idx, nr.problem, nr.transcription);
    
    # Add constraints
    add_continuity!(nr.optimizer, nlp_evaluator.x_idx);
    add_boundary_conditions!(nr.optimizer, nr.problem.x, nlp_evaluator.x_idx);
    add_bounds!(nr.optimizer, nlp_evaluator, nr.transcription.bounds);

    MOI.set(nr.optimizer, MOI.NLPBlock(), nlp_block_data(nlp_evaluator, nr.transcription));
    MOI.set(nr.optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE);
    #MOI.optimize!(nr.optimizer);

    return (_, state)
end

function Base.iterate(nr::NoRefinement, state::NoRefinementState)



    return (_, state)
end





#=
struct NoRefinement{F} <: AbstractRefinement{F}
    problem::DOProblem{F}
    mesh::AbstractMesh{F}
    direction::Progradio.ProgradioDirection{F}
    search::Progradio.ProgradioSearch{F}
    maxIterations::Int
    ϵ_f::F
end

struct NoRefinementState{F} <: AbstractRefinementState{F}
    i::Int
    data::Vector{IntervalData{F}}
    optimality::F
    #trace::Vector{F}
end

function Base.iterate(refinement::NoRefinement{F}) where {F<:AbstractFloat}
    #Initialization
    data = coldstart(refinement.problem, refinement.mesh);
        
    # Transcription
    z_0, z_ℓ, z_u, _, f, g! = transcribe(refinement.problem, refinement.mesh, data);
    bci = Iterator(BCProblem(z_0, z_ℓ, z_u, 0.01, f, g!), refinement.direction, refinement.search; i_max = refinement.maxIterations);

    # Optimization
    fz, bcs = iterate(bci);
    #println(fz)
    i = 0;
    converged = false
    #trace = Vector{F}(undef, 0);
    while !converged
        next = iterate(bci, bcs);
        if next isa Nothing
            break
        else
            i += 1;
            fz, bcs = next;
                if fz <= refinement.ϵ_f
                    #println("Converged!")
                    converged = true
                    break
                end 
            #println("fz = ", fz)
            #append!(trace, fz);
        end
    end

    data = decode(refinement.problem, refinement.mesh, bcs.x);
    optimality = optimality_function(bci.bcp, bcs);
    return fz, NoRefinementState(i, data, optimality)#, trace)
end=#