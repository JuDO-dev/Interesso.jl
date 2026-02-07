function perturb_solution(
    solution::Interesso.PiecewiseInterpolant{I},
    sigma::Real
) where{I<:Interesso.LagrangeInterpolant}
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
    solutions::Dict{String,DOI.AbstractDynamicSolution},
    sigma::Real;
)::Dict{String,DOI.AbstractDynamicSolution}
    return Dict(name => perturb_solution(sol, sigma) for (name, sol) in solutions)
end