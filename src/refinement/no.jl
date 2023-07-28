struct NoRefinement{F, T} <: InteressoRefinement{F, T}
    problem::DFProblem
    transcription::T
    optimiser::ProgradioOptimiser
    maxIterations::Integer
end

mutable struct NoRefinementState
    i::Integer
end

function Base.iterate(refinement<:NoRefinement)
    return (f0, NoRefinementState(0))
end

function Base.iterate(refinement::NoRefinement, state::NoRefinementState)
    # Check maximum iterations
    if state.i >= refinement.maxIterations
        return nothing
    else
        # Transcription
        bcp = transcribe(refinement.problem, refinement.transcription);
        bci = iterator(bcp, refinement.optimiser; i_max = refinement.maxIterations);

        # Optimization
        (fz, bcs) = iterate(bci);
        while state.i < refinement.maxIterations
            
            next = iterate(bci, bcs);
            if next isa Nothing
                break
            else
                state.i += 1;
        end

        return fz
    end
end