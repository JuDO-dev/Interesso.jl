function perturb_solution(
    solution::Interesso.PiecewiseInterpolant{I},
    sigma::Real
) where {I<:Interesso.LagrangeInterpolant}
    perturbed_pieces = [
        Interesso.LagrangeInterpolant(
            piece.initial,
            piece.final,
            copy(piece.points),
            copy(piece.weights),
            piece.values .* (1.0 .+ sigma .* randn(length(piece.values))),
        )
        for piece in solution.pieces
    ]
    return Interesso.PiecewiseInterpolant(perturbed_pieces)
end

function perturb_solutions(
    solutions::Interesso.WSS{T},
    sigma::Real,
)::Interesso.WSS{DOI.AbstractDynamicSolution} where {T<:DOI.AbstractDynamicSolution}
    out = Interesso.WSS{DOI.AbstractDynamicSolution}()

    for (phase, phase_sols) in solutions
        out_phase = Interesso.WS{DOI.AbstractDynamicSolution}()
        for (name, sol) in phase_sols
            out_phase[name] = perturb_solution(sol, sigma)
        end
        out[phase] = out_phase
    end

    return out
end