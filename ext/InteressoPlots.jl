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

# Per-bar widths so neighbouring bars touch instead of collapsing to the
# minimum node gap (LGR nodes cluster, which makes auto-width bars slivers).
function _bar_widths(nodes::AbstractVector)
    n = length(nodes)
    n < 2 && return fill(1.0, n)
    d = diff(nodes)
    return [d; d[end]]
end

function Interesso.plot_residual!(
    plt::Plots.Plot,
    residuals::Vector{Interesso.IntervalResidual};
    aggregate::Bool=false,
    label::String="",
    style::Symbol=:line,
)
    style in (:line, :bar) ||
        throw(ArgumentError("style must be :line or :bar, got :$style."))

    plot_residuals = aggregate ?
        Interesso.aggregate_residuals(residuals; label) :
        residuals

    for res in plot_residuals
        if style === :line
            Plots.plot!(
                plt,
                res.nodes,
                res.residual;
                xlabel = "domain",
                ylabel = "residual",
                label = res.label,
                grid = true,
            )
        else
            Plots.bar!(
                plt,
                res.nodes,
                res.residual;
                bar_width = _bar_widths(res.nodes),
                xlabel = "domain",
                ylabel = "residual",
                label = res.label,
                grid = true,
                linealpha = 0,
                fillalpha = 0.5,
            )
        end
    end

    return plt
end

function Interesso.plot_residual!(
    plt::Plots.Plot,
    model::Interesso.Optimizer;
    q::Integer=10,
    aggregate::Bool=false,
    label::String="",
    style::Symbol=:line,
)
    residuals = Interesso.residual_map(model; q)
    return Interesso.plot_residual!(plt, residuals; aggregate, label, style)
end

function Interesso.plot_residual(
    residuals::Vector{Interesso.IntervalResidual};
    aggregate::Bool=false,
    label::String="",
    style::Symbol=:line,
)
    plt = Plots.plot()
    Interesso.plot_residual!(plt, residuals; aggregate, label, style)
    Plots.plot!(plt; title = "Residual at quadrature points")
    return plt
end

function Interesso.plot_residual(
    model::Interesso.Optimizer;
    q::Integer=10,
    aggregate::Bool=false,
    label::String="",
    style::Symbol=:line,
)
    residuals = Interesso.residual_map(model; q)
    plt = Interesso.plot_residual(residuals; aggregate, label, style)
    Plots.plot!(plt; title = "Residual at interpolated $(q) quadrature points")
    return plt
end

end