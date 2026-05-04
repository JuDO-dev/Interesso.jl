module InteressoPlots

using Interesso
import Plots

function ylims(sol; n=100, minspan=1e-3)

    ts = range(sol.initial, sol.final; length=n)
    ys = sol.(ts)

    ymin = minimum(ys)
    ymax = maximum(ys)
    center = (ymin + ymax) / 2

    span = max(ymax - ymin, minspan)

    return (center - span/2, center + span/2)
end


function Interesso.plot(model::Interesso.Optimizer, var_name::String)
    
    solutions = get_solutions(model)

    plt = Plots.plot(title=var_name, legend=:topright)

    has_solution = false

    for (p, phase) in enumerate(model.phases)
        phase_solutions = get(solutions, phase, nothing)
        phase_solutions === nothing && continue

        sol = get(phase_solutions, var_name, nothing)
        sol === nothing && continue

        Plots.plot!(
            plt,
            τ -> sol(τ),
            sol.initial,
            sol.final,
            ylims=ylims(sol);
            label=false,
            # label="phase $(p)",
        )
        has_solution = true
    end

    has_solution || throw(ArgumentError("No dynamic variable with name '$var_name' was found."))

    return plt
end

function Interesso.plot(model::Interesso.Optimizer)

    solutions = get_solutions(model)
    var_names = unique(values(model.dyn_var_names))

    n_vars = length(var_names)
    n_cols = n_vars > 6 ? 2 : 1
    n_rows = cld(n_vars, n_cols)

    plt = Plots.plot(
        layout=(n_rows, n_cols),
        legend=:topright,
        size=(1000 * n_cols, max(300 * n_rows, 400)),
    )

    for (i, var_name) in enumerate(var_names)
        
        Plots.plot!(plt; title=var_name, subplot=i)

        has_solution = false

        for (p, phase) in enumerate(model.phases)
            phase_solutions = get(solutions, phase, nothing)
            phase_solutions === nothing && continue

            sol = get(phase_solutions, var_name, nothing)
            sol === nothing && continue

            Plots.plot!(
                plt,
                τ -> sol(τ),
                sol.initial,
                sol.final,
                ylims=ylims(sol);
                label=false,
                # label="phase $(p)",
                subplot=i,
            )
            has_solution = true
        end
    end

    return plt
end

function Interesso.plot_residual!(
    plt::Plots.Plot,
    residuals::Vector{Interesso.IntervalResidual};
    aggregate::Bool=false,
    label::String="aggregate",
)
    plot_residuals = aggregate ?
        Interesso.aggregate_residuals(residuals; label) :
        residuals

    for res in plot_residuals
        Plots.plot!(
            plt,
            res.nodes,
            res.residual;
            xlabel = "domain",
            ylabel = "residual",
            label = res.label,
            grid = true,
        )
        # Plots.bar!(
        #     plt,
        #     res.nodes,
        #     res.residual;
        #     xlabel = "domain",
        #     ylabel = "residual",
        #     label = res.label,
        #     grid = true,
        #     linealpha = 0,
        # )
        # Plots.scatter!(
        #     plt,
        #     res.nodes,
        #     res.residual;
        #     xlabel = "domain",
        #     ylabel = "residual",
        #     label = res.label,
        #     grid = true,
        #     markersize = 1,
        #     markeralpha = 0.9,
        #     markerstrokewidth = 0,
        # )
    end

    return plt
end

function Interesso.plot_residual!(
    plt::Plots.Plot,
    model::Interesso.Optimizer;
    q::Integer=10,
    aggregate::Bool=false,
    label::String="aggregate",
)
    residuals = Interesso.residual_map(model; q)
    return Interesso.plot_residual!(plt, residuals; aggregate, label)
end

function Interesso.plot_residual(
    residuals::Vector{Interesso.IntervalResidual};
    aggregate::Bool=false,
    label::String="aggregate",
)
    plt = Plots.plot()
    Interesso.plot_residual!(plt, residuals; aggregate, label)
    Plots.plot!(plt; title = "Residual at quadrature points")
    return plt
end

function Interesso.plot_residual(
    model::Interesso.Optimizer;
    q::Integer=10,
    aggregate::Bool=false,
    label::String="aggregate",
)
    residuals = Interesso.residual_map(model; q)
    plt = Interesso.plot_residual(residuals; aggregate, label)
    Plots.plot!(plt; title = "Residual at interpolated $(q) quadrature points")
    return plt
end

end