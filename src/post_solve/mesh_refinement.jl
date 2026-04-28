"""
    find_split_point(res::IntervalResidual)

Return the time at which the pointwise residual is largest within the interval.
"""
function find_split_point(model::Optimizer, res::IntervalResidual)
    idx = argmax(res.point_res)
    t = res.point_mesh[idx]
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

    # 1. per-interval residuals (with pointwise data)
    residuals = eval_accuracy(model; q)
    R = argmax(r -> r.residual, residuals)

    # 2. split point — just a lookup, no re-evaluation
    τ = find_split_point(model, R)

    # 3. warm-start from current solution
    model.phase_finals[R.phase] isa MOI.EqualTo && normalize_solutions!(model)
    solutions = get_solutions(model)
    warmstart!(model, solutions)

    # 4. read breakpoints from mesh before reset clears it
    interval = get(model.phase_intervals, R.phase, model.default_intervals)
    points = interval.points
    insert!(points, R.interval + 1, τ)

    # 5. reset and install
    reset!(model)
    model.phase_intervals[R.phase] = FixedIntervals(points)

    MOI.optimize!(model)

    return model
end
