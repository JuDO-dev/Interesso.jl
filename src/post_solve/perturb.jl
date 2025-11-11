import DynOptInterface as DOI

function perturb_solution(
    solution::Interesso.PiecewiseInterpolant,
    sigma::Real
)
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
    solutions::Dict{String, DOI.AbstractDynamicSolution},
    sigma::Real;
)
    return Dict(name => perturb_solution(sol, sigma) for (name, sol) in solutions)
end