using CSV
using DataFrames
using Statistics

const PLOT_COLORS = [
    "#1f77b4",
    "#d62728",
    "#2ca02c",
    "#9467bd",
    "#8c564b",
    "#17becf",
]
const PLOT_DASHES = ["", "8,4", "2,4", "10,3,2,3", "4,2,1,2", "12,4"]
const PLOT_MARKERS = [:circle, :square, :triangle, :diamond, :cross, :x]

_svg_escape(x) = replace(string(x), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;")

function _numeric_or_missing(x)
    if x isa Missing || isnothing(x)
        return missing
    elseif x isa Number
        return Float64(x)
    end
    sx = strip(String(x))
    isempty(sx) && return missing
    lower = lowercase(sx)
    lower in ("missing", "na", "nan", "nothing") && return missing
    try
        return Float64(sx)
    catch
        return missing
    end
end

function _bool_or_missing(x)
    if x isa Missing || isnothing(x)
        return missing
    elseif x isa Bool
        return x
    end
    sx = lowercase(strip(String(x)))
    sx in ("true", "1", "yes") && return true
    sx in ("false", "0", "no") && return false
    return missing
end

function _first_existing(df::DataFrame, candidates::Vector{String})
    for name in candidates
        if name in names(df)
            return name
        end
    end
    return nothing
end

function _require_column(df::DataFrame, candidates::Vector{String}, label::String)
    found = _first_existing(df, candidates)
    isnothing(found) && error("The input CSV must contain a column for $(label). Tried: $(join(candidates, ", ")).")
    return found
end

function _solver_config(df::DataFrame)
    if "solver_config" in names(df)
        return [String(v) for v in df[!, "solver_config"]]
    end
    solver_col = _require_column(df, ["solver", "Solver", "solver_family"], "solver name")
    warm_col = _first_existing(df, ["warmstart", "Warmstart"])
    configs = String[]
    for i in 1:nrow(df)
        solver = String(df[i, solver_col])
        if lowercase(solver) == "blclag" && !isnothing(warm_col)
            warm = _bool_or_missing(df[i, warm_col])
            if warm === true
                push!(configs, "BlCLag+BlC ws=on")
            elseif warm === false
                push!(configs, "BlCLag+BlC ws=off")
            else
                push!(configs, "BlCLag+BlC")
            end
        elseif lowercase(solver) in ("blc", "blcsolver")
            push!(configs, "BlC")
        else
            push!(configs, solver)
        end
    end
    return configs
end

function _normalize_benchmark_df(df::DataFrame)
    runtime_col = _require_column(df, ["runtime", "Runtime"], "runtime")
    status_col = _first_existing(df, ["status", "Opt_status", "reported_status"])
    opt_col = _first_existing(df, ["opt", "Opt"])
    bnodes_col = _first_existing(df, ["bnodes", "BNodes"])
    nusers_col = _require_column(df, ["nusers", "instance_nusers"], "number of users")
    alpha_col = _require_column(df, ["alpha", "instance_alpha"], "alpha")
    case_col = _first_existing(df, ["case_id", "instance_name", "experiment_name", "run_id"])
    runtime_limit_col = _first_existing(df, ["runtime_limit", "time_limit"])
    first_var_col = _first_existing(df, ["first_level_vars", "subproblem_first_level_variable_count_total"])
    second_var_col = _first_existing(df, ["second_level_vars", "subproblem_second_level_variable_count_total"])
    linking_var_col = _first_existing(df, ["linking_vars", "subproblem_binary_linking_variable_count_total"])

    out = DataFrame()
    out.solver_config = _solver_config(df)
    out.runtime = [_numeric_or_missing(v) for v in df[!, runtime_col]]
    out.status = isnothing(status_col) ? fill("Unknown", nrow(df)) : [String(v) for v in df[!, status_col]]
    out.opt = isnothing(opt_col) ? fill(missing, nrow(df)) : [_numeric_or_missing(v) for v in df[!, opt_col]]
    out.bnodes = isnothing(bnodes_col) ? fill(missing, nrow(df)) : [_numeric_or_missing(v) for v in df[!, bnodes_col]]
    out.nusers = [Int(round(_numeric_or_missing(v))) for v in df[!, nusers_col]]
    out.alpha = [_numeric_or_missing(v) for v in df[!, alpha_col]]
    out.case_id = isnothing(case_col) ? [string("case_", i) for i in 1:nrow(df)] : [String(v) for v in df[!, case_col]]
    out.runtime_limit = isnothing(runtime_limit_col) ? fill(missing, nrow(df)) : [_numeric_or_missing(v) for v in df[!, runtime_limit_col]]
    out.first_level_vars = isnothing(first_var_col) ? fill(missing, nrow(df)) : [_numeric_or_missing(v) for v in df[!, first_var_col]]
    out.second_level_vars = isnothing(second_var_col) ? fill(missing, nrow(df)) : [_numeric_or_missing(v) for v in df[!, second_var_col]]
    out.linking_vars = isnothing(linking_var_col) ? fill(missing, nrow(df)) : [_numeric_or_missing(v) for v in df[!, linking_var_col]]
    return out
end

function _plot_style_map(labels::Vector{String})
    unique_labels = sort(unique(labels))
    return Dict(
        label => (
            color=PLOT_COLORS[mod1(i, length(PLOT_COLORS))],
            dash=PLOT_DASHES[mod1(i, length(PLOT_DASHES))],
            marker=PLOT_MARKERS[mod1(i, length(PLOT_MARKERS))]
        ) for (i, label) in enumerate(unique_labels)
    )
end

function _svg_header(io, width, height)
    println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$width" height="$height" viewBox="0 0 $width $height">""")
    println(io, """<rect x="0" y="0" width="$width" height="$height" fill="white"/>""")
    println(io, """<style>
        text { font-family: Arial, Helvetica, sans-serif; fill: #111; }
        .title { font-size: 20px; font-weight: bold; }
        .axislabel { font-size: 14px; font-weight: bold; }
        .tick { font-size: 11px; }
        .legend { font-size: 12px; }
    </style>""")
end

function _svg_footer(io)
    println(io, "</svg>")
end

function _draw_axes(io, x0, y0, x1, y1; xlabel="", ylabel="", title="")
    println(io, """<line x1="$x0" y1="$y1" x2="$x1" y2="$y1" stroke="#222" stroke-width="1.5"/>""")
    println(io, """<line x1="$x0" y1="$y0" x2="$x0" y2="$y1" stroke="#222" stroke-width="1.5"/>""")
    println(io, """<text x="$(x0 + (x1-x0)/2)" y="$(y0 - 22)" text-anchor="middle" class="title">$(_svg_escape(title))</text>""")
    println(io, """<text x="$(x0 + (x1-x0)/2)" y="$(y1 + 38)" text-anchor="middle" class="axislabel">$(_svg_escape(xlabel))</text>""")
    println(io, """<text x="$(x0 - 10)" y="$(y0 - 10)" text-anchor="start" class="axislabel">$(_svg_escape(ylabel))</text>""")
end

function _draw_y_ticks(io, x0, y0, y1, ymin, ymax; nticks=5, fmt=x -> round(x; digits=2))
    ymax == ymin && (ymax = ymin + 1)
    for i in 0:nticks
        v = ymin + (ymax - ymin) * i / nticks
        y = y1 - (y1 - y0) * i / nticks
        println(io, """<line x1="$x0" y1="$y" x2="$(x0-5)" y2="$y" stroke="#222" stroke-width="1"/>""")
        println(io, """<line x1="$x0" y1="$y" x2="$(x0 + (y1-y0)*1.2)" y2="$y" stroke="#e5e5e5" stroke-width="1"/>""")
        println(io, """<text x="$(x0-8)" y="$(y+4)" text-anchor="end" class="tick">$(_svg_escape(fmt(v)))</text>""")
    end
end

function _marker_svg(marker, x, y, color)
    if marker == :circle
        return """<circle cx="$x" cy="$y" r="4" fill="$color" stroke="#111" stroke-width="1"/>"""
    elseif marker == :square
        return """<rect x="$(x-4)" y="$(y-4)" width="8" height="8" fill="$color" stroke="#111" stroke-width="1"/>"""
    elseif marker == :triangle
        pts = "$(x),$(y-5) $(x-5),$(y+4) $(x+5),$(y+4)"
        return """<polygon points="$pts" fill="$color" stroke="#111" stroke-width="1"/>"""
    elseif marker == :diamond
        pts = "$(x),$(y-5) $(x-5),$y $(x),$(y+5) $(x+5),$y"
        return """<polygon points="$pts" fill="$color" stroke="#111" stroke-width="1"/>"""
    elseif marker == :cross
        return """<g stroke="$color" stroke-width="2"><line x1="$(x-4)" y1="$y" x2="$(x+4)" y2="$y"/><line x1="$x" y1="$(y-4)" x2="$x" y2="$(y+4)"/></g>"""
    else
        return """<g stroke="$color" stroke-width="2"><line x1="$(x-4)" y1="$(y-4)" x2="$(x+4)" y2="$(y+4)"/><line x1="$(x-4)" y1="$(y+4)" x2="$(x+4)" y2="$(y-4)"/></g>"""
    end
end

function _draw_legend(io, styles::AbstractDict{String,<:NamedTuple}, labels::Vector{String}, x, y)
    for (idx, label) in enumerate(labels)
        style = styles[label]
        yy = y + 22 * (idx - 1)
        println(io, """<line x1="$x" y1="$yy" x2="$(x+28)" y2="$yy" stroke="$(style.color)" stroke-width="2.5" stroke-dasharray="$(style.dash)"/>""")
        println(io, _marker_svg(style.marker, x + 14, yy, style.color))
        println(io, """<text x="$(x+38)" y="$(yy+4)" class="legend">$(_svg_escape(label))</text>""")
    end
end

function _performance_profile(df::DataFrame, output_path::String)
    solved_labels = sort(unique(df.solver_config))
    styles = _plot_style_map(solved_labels)
    runtime_limit_default = maximum(skipmissing(coalesce.(df.runtime_limit, df.runtime, 1.0)))::Float64
    instances = sort(unique(df.case_id))
    ratios = Dict(label => Float64[] for label in solved_labels)

    for case_id in instances
        sub = df[df.case_id .== case_id, :]
        effective = Dict{String,Float64}()
        for row in eachrow(sub)
            status = lowercase(String(row.status))
            runtime = ismissing(row.runtime) ? runtime_limit_default : Float64(row.runtime)
            limit = ismissing(row.runtime_limit) ? runtime_limit_default : Float64(row.runtime_limit)
            effective[row.solver_config] = status == "optimal" ? runtime : 2.0 * limit
        end
        best = minimum(values(effective))
        for label in solved_labels
            push!(ratios[label], effective[label] / best)
        end
    end

    for label in solved_labels
        sort!(ratios[label])
    end

    width, height = 980, 620
    x0, y0, x1, y1 = 90, 80, 730, 520
    max_ratio = maximum(maximum(vals) for vals in values(ratios))

    open(output_path, "w") do io
        _svg_header(io, width, height)
        _draw_axes(io, x0, y0, x1, y1; xlabel="Performance Ratio (runtime / best runtime on instance)", ylabel="Solved Fraction", title="Performance Profile")
        _draw_y_ticks(io, x0, y0, y1, 0.0, 1.0; nticks=5, fmt=x -> string(round(x; digits=2)))
        for i in 0:5
            r = 1.0 + (max_ratio - 1.0) * i / 5
            x = x0 + (x1 - x0) * (r - 1.0) / max(1e-9, (max_ratio - 1.0))
            println(io, """<line x1="$x" y1="$y1" x2="$x" y2="$(y1+5)" stroke="#222" stroke-width="1"/>""")
            println(io, """<text x="$x" y="$(y1+20)" text-anchor="middle" class="tick">$(_svg_escape(round(r; digits=2)))</text>""")
        end
        for label in solved_labels
            style = styles[label]
            vals = ratios[label]
            n = length(vals)
            pts = Tuple{Float64,Float64}[]
            push!(pts, (1.0, 0.0))
            for (i, v) in enumerate(vals)
                frac_prev = (i - 1) / n
                frac = i / n
                push!(pts, (v, frac_prev))
                push!(pts, (v, frac))
            end
            path = IOBuffer()
            first = true
            for (rx, fy) in pts
                x = x0 + (x1 - x0) * (rx - 1.0) / max(1e-9, (max_ratio - 1.0))
                y = y1 - (y1 - y0) * fy
                print(path, first ? "M $x $y " : "L $x $y ")
                first = false
            end
            println(io, """<path d="$(String(take!(path)))" fill="none" stroke="$(style.color)" stroke-width="2.5" stroke-dasharray="$(style.dash)"/>""")
        end
        _draw_legend(io, styles, solved_labels, 770, 120)
        _svg_footer(io)
    end
end

function _avg_bnodes_plot(df::DataFrame, output_path::String)
    usable = df[.!ismissing.(df.bnodes), :]
    labels = sort(unique(usable.solver_config))
    styles = _plot_style_map(labels)
    avg_nodes = [mean(skipmissing(usable[usable.solver_config .== label, :bnodes])) for label in labels]
    ymax = maximum(avg_nodes) * 1.15

    width, height = 920, 560
    x0, y0, x1, y1 = 90, 80, 820, 470
    barw = 0.65 * (x1 - x0) / max(1, length(labels))

    open(output_path, "w") do io
        _svg_header(io, width, height)
        _draw_axes(io, x0, y0, x1, y1; xlabel="Solver / Parameter Combination", ylabel="Average B&B Nodes", title="Average Branch-and-Bound Nodes")
        _draw_y_ticks(io, x0, y0, y1, 0.0, ymax; nticks=5, fmt=x -> string(round(Int, x)))
        for (i, label) in enumerate(labels)
            style = styles[label]
            cx = x0 + (i - 0.5) * (x1 - x0) / length(labels)
            h = (y1 - y0) * avg_nodes[i] / max(1e-9, ymax)
            println(io, """<rect x="$(cx - barw/2)" y="$(y1 - h)" width="$barw" height="$h" fill="$(style.color)" fill-opacity="0.65" stroke="#111" stroke-width="1.2"/>""")
            println(io, """<text x="$cx" y="$(y1 + 18)" text-anchor="middle" class="tick">$(_svg_escape(label))</text>""")
            println(io, """<text x="$cx" y="$(y1 - h - 8)" text-anchor="middle" class="tick">$(_svg_escape(round(Int, avg_nodes[i])))</text>""")
        end
        _svg_footer(io)
    end
end

function _quartiles(values::Vector{Float64})
    vals = sort(values)
    return (
        minimum(vals),
        quantile(vals, 0.25),
        quantile(vals, 0.5),
        quantile(vals, 0.75),
        maximum(vals),
    )
end

function _boxplot_panel(io, x0, y0, x1, y1, labels, values_by_label, styles, title, ylabel)
    all_vals = reduce(vcat, values_by_label)
    ymin = minimum(all_vals)
    ymax = maximum(all_vals)
    ymax == ymin && (ymax = ymin + 1.0)
    _draw_axes(io, x0, y0, x1, y1; xlabel="Solver / Parameter Combination", ylabel=ylabel, title=title)
    _draw_y_ticks(io, x0, y0, y1, ymin, ymax; nticks=5, fmt=x -> string(round(Int, x)))
    boxw = 0.55 * (x1 - x0) / max(1, length(labels))
    for (i, label) in enumerate(labels)
        style = styles[label]
        vals = values_by_label[i]
        qmin, q1, q2, q3, qmax = _quartiles(vals)
        cx = x0 + (i - 0.5) * (x1 - x0) / length(labels)
        function ymap(v)
            y1 - (y1 - y0) * (v - ymin) / max(1e-9, (ymax - ymin))
        end
        y_q1, y_q2, y_q3 = ymap(q1), ymap(q2), ymap(q3)
        y_min, y_max = ymap(qmin), ymap(qmax)
        println(io, """<line x1="$cx" y1="$y_max" x2="$cx" y2="$y_q3" stroke="#111" stroke-width="1.2"/>""")
        println(io, """<line x1="$cx" y1="$y_q1" x2="$cx" y2="$y_min" stroke="#111" stroke-width="1.2"/>""")
        println(io, """<line x1="$(cx - boxw/4)" y1="$y_max" x2="$(cx + boxw/4)" y2="$y_max" stroke="#111" stroke-width="1.2"/>""")
        println(io, """<line x1="$(cx - boxw/4)" y1="$y_min" x2="$(cx + boxw/4)" y2="$y_min" stroke="#111" stroke-width="1.2"/>""")
        println(io, """<rect x="$(cx - boxw/2)" y="$y_q3" width="$boxw" height="$(y_q1 - y_q3)" fill="$(style.color)" fill-opacity="0.55" stroke="#111" stroke-width="1.2"/>""")
        println(io, """<line x1="$(cx - boxw/2)" y1="$y_q2" x2="$(cx + boxw/2)" y2="$y_q2" stroke="#111" stroke-width="1.6" stroke-dasharray="$(style.dash)"/>""")
        println(io, """<text x="$cx" y="$(y1 + 18)" text-anchor="middle" class="tick">$(_svg_escape(label))</text>""")
    end
end

function _variable_boxplots(df::DataFrame, output_path::String)
    labels = sort(unique(df.solver_config))
    styles = _plot_style_map(labels)
    metrics = [
        ("First-Level Variables", :first_level_vars),
        ("Second-Level Variables", :second_level_vars),
        ("Linking Variables", :linking_vars),
    ]

    width, height = 1320, 520
    open(output_path, "w") do io
        _svg_header(io, width, height)
        panel_width = 360
        starts = [40, 470, 900]
        for (panel_idx, (title, col)) in enumerate(metrics)
            labels_local = String[]
            values_local = Vector{Vector{Float64}}()
            for label in labels
                vals = Float64[]
                for v in df[df.solver_config .== label, col]
                    if !ismissing(v)
                        push!(vals, Float64(v))
                    end
                end
                isempty(vals) && continue
                push!(labels_local, label)
                push!(values_local, vals)
            end
            _boxplot_panel(
                io,
                starts[panel_idx],
                80,
                starts[panel_idx] + panel_width,
                430,
                labels_local,
                values_local,
                styles,
                title,
                "Count",
            )
        end
        _svg_footer(io)
    end
end

function _runtime_by_users(df::DataFrame, output_path::String)
    labels = sort(unique(df.solver_config))
    styles = _plot_style_map(labels)
    alphas = sort(unique(skipmissing(df.alpha)))
    width, height = 1240, 460

    open(output_path, "w") do io
        _svg_header(io, width, height)
        panel_w = 320
        starts = [60 + (i - 1) * 390 for i in 1:length(alphas)]
        for (panel_idx, alpha) in enumerate(alphas)
            sub = df[df.alpha .== alpha, :]
            x0, y0, x1, y1 = starts[panel_idx], 80, starts[panel_idx] + panel_w, 380
            ymax = maximum(skipmissing(sub.runtime))
            ymax = max(ymax, maximum(skipmissing(coalesce.(sub.runtime_limit, sub.runtime)))) * 1.1
            _draw_axes(io, x0, y0, x1, y1; xlabel="Number of Users", ylabel="Runtime [s]", title="Mean Runtime, α=$(round(alpha; digits=2))")
            _draw_y_ticks(io, x0, y0, y1, 0.0, ymax; nticks=5, fmt=x -> string(round(Int, x)))
            user_counts = sort(unique(sub.nusers))
            for (j, nusers) in enumerate(user_counts)
                x = x0 + (x1 - x0) * (j - 1) / max(1, length(user_counts) - 1)
                println(io, """<line x1="$x" y1="$y1" x2="$x" y2="$(y1+5)" stroke="#222" stroke-width="1"/>""")
                println(io, """<text x="$x" y="$(y1+20)" text-anchor="middle" class="tick">$nusers</text>""")
            end
            for label in labels
                style = styles[label]
                points = Tuple{Float64,Float64}[]
                for (j, nusers) in enumerate(user_counts)
                    vals = [
                        Float64(r) for r in sub[(sub.solver_config .== label) .& (sub.nusers .== nusers), :runtime]
                        if !ismissing(r)
                    ]
                    isempty(vals) && continue
                    x = x0 + (x1 - x0) * (j - 1) / max(1, length(user_counts) - 1)
                    y = y1 - (y1 - y0) * mean(vals) / max(1e-9, ymax)
                    push!(points, (x, y))
                end
                if length(points) >= 2
                    path = join(["$(i == 1 ? "M" : "L") $(p[1]) $(p[2])" for (i, p) in enumerate(points)], " ")
                    println(io, """<path d="$path" fill="none" stroke="$(style.color)" stroke-width="2.5" stroke-dasharray="$(style.dash)"/>""")
                end
                for (x, y) in points
                    println(io, _marker_svg(style.marker, x, y, style.color))
                end
            end
        end
        _draw_legend(io, styles, labels, 1060, 95)
        _svg_footer(io)
    end
end

function generate_hndp_benchmark_plots(csv_path::AbstractString; output_dir::AbstractString=joinpath(dirname(String(csv_path)), "plots"))
    mkpath(output_dir)
    df_raw = CSV.read(String(csv_path), DataFrame)
    df = _normalize_benchmark_df(df_raw)
    _performance_profile(df, joinpath(output_dir, "performance_profile.svg"))
    _avg_bnodes_plot(df, joinpath(output_dir, "average_bnodes.svg"))
    _variable_boxplots(df, joinpath(output_dir, "variable_boxplots.svg"))
    _runtime_by_users(df, joinpath(output_dir, "runtime_by_users.svg"))
    return output_dir
end

if abspath(PROGRAM_FILE) == @__FILE__
    if isempty(ARGS)
        error("Usage: julia --project=. examples/HNDP/hndp_benchmark_plots.jl <benchmark_csv> [output_dir]")
    end
    output_dir = length(ARGS) >= 2 ? ARGS[2] : joinpath(dirname(ARGS[1]), "plots")
    out = generate_hndp_benchmark_plots(ARGS[1]; output_dir=output_dir)
    println("Wrote plots to " * out)
end
