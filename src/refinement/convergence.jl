struct ConvergenceRefinement{F} <: InteressoRefinement{F}
    



mutable struct ConvergenceIteratorState{F} <: InteressoIteratorState{F}

    integrated_residuals::F

end

function refinement_state(dfi::DFIterator{T, O, R}) where {T, O, R<:ConvergenceIteratorState}


function iterate!(dfi::DFIterator{T, F}, state::ConvergenceIteratorState{F}) 
    
    bcp = transcribe(dfi.dfp, dfi.transcription)


    return nothing
end

function iterate!(doi::DOIterator{T, F}, state::ConvergenceIteratorState{F})
    return nothing
end