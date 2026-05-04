using Plots

include(joinpath(@__DIR__, "../tracks/get_track.jl"))

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
        n_sol  = get(phase_sols, "e_y", nothing)
        vx_sol = get(phase_sols, "v_x", nothing)
        vy_sol = get(phase_sols, "v_y", nothing)

        n_eval = [n_sol(s) for s in sref]
        v_eval = [sqrt(vx_sol(s)^2 + vy_sol(s)^2) for s in sref]

        x_traj = xref .- n_eval .* sin.(ψref)
        y_traj = yref .+ n_eval .* cos.(ψref)

        Plots.plot!(plt, x_traj, y_traj; color = :grey14, linewidth = 0.4, label = "")

        scatter!(plt, x_traj, y_traj;
            marker_z = v_eval, color = :imola,
            markersize = 1.0, markerstrokewidth = 0,
            colorbar = true, colorbar_title = "\nv [m/s]", right_margin = 10Plots.mm, label = "")
    end

    return plt
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