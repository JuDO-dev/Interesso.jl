struct FilteringIterator{T<:InteressoTranscription, F} <: InteressoIterator{F}
    problem::DFProblem
    optimiser::ProgradioOptimiser
    maxIterations::Integer
    maxRefinements::Integer
    #constant
end

struct FilteringIteratorState{F} <: InteressoIteratorState{F}
    i::Integer
    r::Integer
    mesh::Vector{F}
    solution::JuDOBase.DFSolution{F}
    θ::F
end

function Base.iterate(I::FilteringIterator)
    

    return (Fz, state)
end

function Base.iterate(I::FilteringIterator{T,F}, state::FilteringIteratorState{F}) where {T<:InteressoTranscription, F<:AbstractFloat}
    # Check maximum iterations/refinements
    if state.i >= I.maxIterations || state.r >= I.maxRefinements
        return nothing

    else
        # Mesh Refinement
        mesh = refine(state.mesh);


        # Transcription
        sbi = transcribe(T, I.problem, I.optimiser, mesh, state.solution);
        fz, sbs = iterate(sbi);

        # Optimisation
        θ = optimality(sbi, sbs);
        while θ### WIP
            (fz, sbs) = iterate(sbi, sbs);
            i += 1;
            θ = Progradio.optimality(sbi, sbs);
        end

        # Solution
        sol = solution(mesh, sbs);

        return (fz, FilteringIteratorState(state.i, state.r - 1, mesh, sol, θ))
    end
end