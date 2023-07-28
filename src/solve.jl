function solve(problem::JuDOProblem, iterator::InteressoIterator)

    # First iteration
    next = iterate(iterator);

    # Subsequent iterations
    while next !== nothing
        (fz, state) = next;
        next = iterate(iterator, state);
    end

    return solution(problem, state)
end

function solution(dfp::DFProblem, state::InteressoIteratorState)



    return dfsolution()
end

function solution(dop::DOProblem, state::InteressoIteratorState)


    return dosolution()
end