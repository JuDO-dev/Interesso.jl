function coldstart(problem::JuDOProblem{F}, transcription::InteressoTranscription{F}, indices::Vector{UnitRange}) where {F<:AbstractFloat}
    return zeros(F, indices[end][end])
end

function warmstart()



    return nothing
end