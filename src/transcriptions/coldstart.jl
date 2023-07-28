function coldstart(problem::DOProblem{F}, M::RigidMesh{F}) where {F<:AbstractFloat}
    n_x = length(problem.x);
    n_u = length(problem.u);

    data = [IntervalData(undef, n_x, n_u, M) for _ in 1:length(M.H)];

    for n in 1:length(problem.x)
        x_initial = problem.x[n].initial_guess;#(problem.x[n].initial[1] + problem.x[n].initial[2]) / 2;
        x_terminal = problem.x[n].final_guess;#(problem.x[n].terminal[1] + problem.x[n].terminal[1]) / 2;

        for i in 1:length(M.H)
            data[i].x[n] = collect(range(
                start = x_initial + (i - 1) * (x_terminal - x_initial) / length(M.H),
                stop = x_initial + i * (x_terminal - x_initial) / length(M.H),
                length = M.P_x.degree + 1));
        end
    end
    for n in 1:length(problem.u)
        for i in 1:length(M.H)
            data[i].u[n] = zeros(F, M.P_u.degree + 1)
        end
    end
    return data
end