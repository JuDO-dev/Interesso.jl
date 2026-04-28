using Plots

include(joinpath(@__DIR__, "tracks", "get_track.jl"))

function plot_trajectory(model::Interesso.Optimizer, trackfile::AbstractString)
    plt = plot_track(trackfile)
    plot_trajectory!(plt, model, trackfile)
    return plt
end

function plot_track(trackfile::AbstractString)
    _, xref, yref, ψref, _, nlref, nrref, _ = getTrack(trackfile)

    xleft = xref .- nlref .* sin.(ψref)
    yleft = yref .+ nlref .* cos.(ψref)
    xright = xref .+ nrref .* sin.(ψref)
    yright = yref .- nrref .* cos.(ψref)

    plt = Plots.plot(
        xref,
        yref;
        linestyle = :dash,
        color = :black,
        linewidth = 0.2,
        label = "",
        xlabel = "x [m]",
        ylabel = "y [m]",
        aspect_ratio = :equal,
        legend = false,
    )
    Plots.plot!(plt, xleft, yleft; color = :black, linewidth = 0.5)
    Plots.plot!(plt, xright, yright; color = :black, linewidth = 0.5)

    return plt
end

function plot_trajectory!(plt::Plots.Plot, model::Interesso.Optimizer, trackfile::AbstractString)

    solutions = Interesso.get_solutions(model)
    sref, xref, yref, ψref, _, _, _, _ = getTrack(trackfile)

    for phase in model.phases
        phase_sols = get(solutions, phase, nothing)
        n_sol = get(phase_sols, "n", nothing)
        v_sol = get(phase_sols, "v", nothing)

        n_eval = [n_sol(s) for s in sref]
        v_eval = [v_sol(s) for s in sref]

        x_traj = xref .- n_eval .* sin.(ψref)
        y_traj = yref .+ n_eval .* cos.(ψref)

        Plots.plot!(plt, x_traj, y_traj; color = :grey14, linewidth = 0.4, label = "")

        scatter!(plt, x_traj, y_traj;
            marker_z = v_eval, color = :imola,
            markersize = 0.5, markerstrokewidth = 0,
            colorbar = true, colorbar_title = "\nv [m/s]", right_margin = 10Plots.mm, label = "")
    end

    return plt
end

# ────────────────────────────────────────────────────────────────────
# plot_solution — state and control time-histories vs arc-length
# ────────────────────────────────────────────────────────────────────
function plot_solution(model::Interesso.Optimizer; N::Int = 500)

    solutions = Interesso.get_solutions(model)

    var_specs = [
        ("t", "t [s]",   :x, nothing),
        ("n", "n [m]",   :x, [-0.12, 0.12]),
        ("α", "α [rad]", :x, nothing),
        ("v", "v [m/s]", :x, nothing),
        ("D", "D",       :u, [-1.0, 1.0]),
        ("δ", "δ [rad]", :u, [-0.40, 0.40]),
    ]

    subplots = Plots.Plot[]

    for (name, ylabel, tag, bounds) in var_specs
        sp = Plots.plot(; ylabel = ylabel, label = "", linewidth = 1.5)
        if bounds !== nothing
            hline!(sp, bounds; linestyle = :dash, color = :red, label = "")
        end

        for phase in model.phases
            phase_sols = get(solutions, phase, nothing)
            phase_sols === nothing && continue
            sol = get(phase_sols, name, nothing)
            sol === nothing && continue

            s_eval = range(sol.initial, sol.final; length = N)
            y_eval = [sol(s) for s in s_eval]
            if tag == :x
                Plots.plot!(sp, collect(s_eval), y_eval; color = :dodgerblue, label = "", linewidth = 1.5)
            elseif tag == :u
                Plots.plot!(sp, collect(s_eval), y_eval; color = :orange, label = "", linewidth = 1.5)
            end
        end

        push!(subplots, sp)
    end

    return Plots.plot(subplots...;
        layout = (length(subplots), 1), size = (800, 200*length(subplots)), link = :x,
        left_margin = 8Plots.mm)
end

# ────────────────────────────────────────────────────────────────────
# plot_curvature — curvature spline vs track data
# ────────────────────────────────────────────────────────────────────

function plot_curvature!(plt::Plots.Plot, trackfile::AbstractString)
    sref, _, _, _, κref, _, _, _ = getTrack(trackfile)

    Plots.plot!(plt, sref, κref;
        linewidth = 2,
        label = "interpolated κ(s)",
        xlabel = "s [m]",
        ylabel = "κ(s) [1/m]",
        size = (1000, 400),
        left_margin = 5Plots.mm,
        bottom_margin = 5Plots.mm,
    )
    scatter!(plt, sref, κref; markersize = 2, label = "data")
    return plt
end

function plot_curvature(trackfile::AbstractString)
    plt = Plots.plot()
    plot_curvature!(plt, trackfile)
    return plt
end