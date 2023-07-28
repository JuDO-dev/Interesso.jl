struct PredictiveRefinement{F} <: InteressoRefinement{F}
    problem::DOProblem{F}
    M_0::AbstractMesh{F}
    optimizer::Progradio.ProgradioOptimizer
    maxIterations::Integer
    maxRefinements::Integer
    maxPartitions::Integer
    ϵ_f::F
    ϵ_θ::F
end

mutable struct PredictiveRefinementState{F} <: InteressoRefinementState{F}
    r::Integer
    i::Integer
    M::AbstractMesh{F}
    data::Vector{IntervalData{F}}
    optimality::F
    trace::Vector{F}
end

function Base.iterate(refinement::PredictiveRefinement{F}) where {F<:AbstractFloat}
    #Initialization
    return Inf, PredictiveRefinementState(0, 0, deepcopy(refinement.M_0), coldstart(refinement.problem, refinement.M_0), -Inf, Vector{F}(undef, 0))
end

function Base.iterate(refinement::PredictiveRefinement{F}, state::PredictiveRefinementState{F}) where {F<:AbstractFloat}
    # Check maximum refinements
    if state.r >= refinement.maxRefinements
        return nothing
    else
        # Transcription
        z_0, z_ℓ, z_u, f_i, f, g! = transcribe(refinement.problem, state.M, state.data);
        bci = iterator(BCProblem(f, g!, z_ℓ, z_u, z_0, 0.01), refinement.optimizer; i_max = refinement.maxIterations - state.i);

        # Optimization
        fz, bcs = iterate(bci);
        θ = optimality_function(bci.bcp, bcs);
        #println("fz = ", fz)
        #println("θz = ", θ)
        B = deepcopy(bcs.B);
        fine = true
        converged = false
        while fine && !converged
            next = iterate(bci, bcs);
            if next isa Nothing
                break
            else
                state.i += 1;
                fz_next, bcs = next;
                if fz_next <= refinement.ϵ_f
                    converged = true
                    #println("Converged!")
                    break
                end
                θ_next = optimality_function(bci.bcp, bcs);

                #Print stuff:
                #println("fz = ", fz_next)
                #println("θz = ", θ_next)
                #println("ΣΔB = ", sum(bcs.B .!= B))

                # Convergence
                if θ_next > refinement.ϵ_θ

                    # Binding Set changes
                    if sum(bcs.B .!= B) == 0
                        fine = false;
                        #println("Refine!")

                        # Partition H
                        partitions = repeat([1], length(state.M.H));
                        indices = indexing(state.M, length(refinement.problem.x), length(refinement.problem.u));
                        data = decode(refinement.problem, state.M, bcs.x);

                        for i in 1:length(state.M.H)
                            fz_i = f_i(bcs.x[indices[i]], state.M.H[i]);
                            ϵ_fz_i = fz_i / state.M.H[i];
                            if ϵ_fz_i > refinement.ϵ_f
                                partitions[i] = min(ceil((ϵ_fz_i / refinement.ϵ_f)^(1 / (state.M.P_x.degree + 1))), refinement.maxPartitions);
                            end
                        end

                        state.M.H, state.data = partition(state.M, data, partitions);
                        state.r += 1;
                        state.optimality = -Inf;
                        return fz, state
                    end
                end
                #append!(state.trace, fz_next);
                fz = fz_next;
                θ = θ_next;
                B = deepcopy(bcs.B);
            end
        end
                
        state.data = decode(refinement.problem, state.M, bcs.x);
        state.optimality = θ;
        return fz, state
    end
end