struct Collocation <: InteressoTranscription
    P_x
    P_u
end

function transcribe(dfp::DFProblem, rm::RigidMesh, c::Collocation)

    return BCProblem()
end