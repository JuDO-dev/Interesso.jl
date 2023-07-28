abstract type InteressoIterator{F<:AbstractFloat} end
abstract type InteressoIteratorState{F<:AbstractFloat} end
Base.eltype(::InteressoIterator{F}) where F<:AbstractFloat = F;
Base.length(iterator::InteressoIterator) = iterator.r_max;

# Dynamic Feasibility
struct DFIterator{F, R<:InteressoRefinement} <: InteressoIterator{F}
    # Main
    dfp::DFProblem{F}
    refinement::R
    # Termination
    r_max
    i_max
end

function iterator(dfp::DFProblem, refinement::InteressoRefinement; i_max)
    return DFIterator{T}(dfp, refinement, r_max, i_max)
end

function Base.iterate(dfi::DFIterator{T, O, R}) 

    #agnostic to refinement, transcription, optimizer



    return f0, refinement_state(dfi)
end

function Base.iterate(dfi::DFIterator, state::InteressoIteratorState)
    if state.r >= dfi.r_max || state.i >= dfi.i_max
        return nothing
    else
        iterate!(dfi, state);
    end
end

# Dynamic Optimization
struct DOIterator{T<:InteressoTranscription, C<:ProgradioConstrainer, R<:InteressoRefinement, F} <: InteressoIterator{F}
    # Main
    dop::DOProblem{F}
    transcription::T
    constrainer::C
    refinement::R
    # Termination
    r_max
    i_max
end