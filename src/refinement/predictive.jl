struct PredictiveRefinement{F, T} <: InteressoRefinement{F, T}
    problem::DFProblem
    transcription::T
    optimiser::ProgradioOptimiser
    maxIterations::Integer
    maxRefinements::Integer
end

function Base.iterate()

end

function Base.iterate()
    
end