"""
    find_split_point(res::IntervalResidual)

Return the time at which the pointwise residual is largest within the interval.
"""
function find_split_point(model::Optimizer, res::IntervalResidual)
    interval = argmax(R.interval_residual)
    q = div(length(res.nodes), length(res.interval_residual))

    first_node = (interval - 1) * q + 1
    last_node = interval * q

    node = first_node + argmax(res.residual[first_node:last_node]) - 1
    t = res.nodes[node]

    if model.phase_finals[res.phase] isa MOI.EqualTo
        t_0 = model.phase_initials[res.phase].value
        t_f = model.phase_finals[res.phase].value
        return (t - t_0) / (t_f - t_0)
    else
        return t
    end
end


"""
    refine!(model::Optimizer; q::Integer=10)

In-place mesh refinement:

1. Calls `eval_accuracy` to get per-interval residuals (with pointwise data).
2. Identifies the interval with the largest integrated residual.
3. Reads the split point from the stored pointwise residuals.
4. Stores the current solution as interpolant warm-starts.
5. Reads breakpoints from the mesh, inserts the split point.
6. Resets the transcription and installs the new `FixedIntervals`.

After this returns, call `MOI.optimize!(model)` to re-solve.
"""
function refine!(model::Optimizer; q::Integer=10)

    residuals = eval_accuracy(model; q)
    aggregates = aggregate_residuals(residuals)

    R = argmax(r -> maximum(r.interval_residual), aggregates)
    interval = argmax(R.interval_residual)
    τ = find_split_point(model, R)

    model.phase_finals[R.phase] isa MOI.EqualTo && normalize_solutions!(model)
    solutions = get_solutions(model)
    warmstart!(model, solutions)

    phase_interval = get(model.phase_intervals, R.phase, model.default_intervals)
    points = phase_interval.points
    insert!(points, interval + 1, τ)

    reset!(model)
    model.phase_intervals[R.phase] = FixedIntervals(points)

    MOI.optimize!(model)

    return model
end
