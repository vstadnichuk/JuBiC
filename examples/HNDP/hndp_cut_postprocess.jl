using CSV
using DataFrames
using Statistics

const CUT_LINE_REGEX = r"ScalarConstraint\{AffExpr, MathOptInterface\.LessThan\{Float64\}\}\((.*), MathOptInterface\.LessThan\{Float64\}\(([-+]?\d+(?:\.\d+)?)\)\)"
const BLC_USER_REGEX = r"optL2_U(\d+)"
const BLCLAG_USER_REGEX = r"subObj\[U(\d+)\]"
const ARC_COEFF_REGEX = r"([-+]?\d+(?:\.\d+)?) x\[\((\d+),\s*(\d+)\)\]"
const CASE_REGEX = r"(layered_U(\d+)_A([0-9.]+)_S(\d+))"

const POST_COLORS = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#8c564b", "#17becf"]
const POST_DASHES = ["", "8,4", "2,4", "10,3,2,3", "4,2,1,2", "12,4"]

_svg_escape(x) = replace(string(x), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;")

function _solver_config_from_path(path::AbstractString)
    occursin("blclag_blc_ws1", path) && return "BlCLag+BlC ws=on"
    occursin("blclag_blc_ws0", path) && return "BlCLag+BlC ws=off"
    occursin("blc", path) && return "BlC"
    return "Unknown"
end

function _parse_case_metadata(path::AbstractString)
    m = match(CASE_REGEX, path)
    isnothing(m) && error("Could not infer case metadata from path: $path")
    return (
        case_id = m.captures[1],
        nusers = parse(Int, m.captures[2]),
        alpha = parse(Float64, m.captures[3]),
        seed = parse(Int, m.captures[4]),
    )
end

function _collect_mastercut_files(benchmark_root::AbstractString)
    files = String[]
    for (root, _, names) in walkdir(joinpath(String(benchmark_root), "runs"))
        for name in names
            if name == "mastercuts_collection.txt"
                push!(files, joinpath(root, name))
            end
        end
    end
    return sort(files)
end

function _parse_mastercuts(files::Vector{String})
    rows = Dict{String,Any}[]
    for filepath in files
        solver_config = _solver_config_from_path(filepath)
        solver_config == "Unknown" && continue
        meta = _parse_case_metadata(filepath)
        cut_idx = 0
        for line in eachline(filepath)
            cut_idx += 1
            m = match(CUT_LINE_REGEX, line)
            isnothing(m) && continue
            expr = m.captures[1]
            rhs = parse(Float64, m.captures[2])
            user_match = solver_config == "BlC" ? match(BLC_USER_REGEX, expr) : match(BLCLAG_USER_REGEX, expr)
            isnothing(user_match) && continue
            user_id = parse(Int, user_match.captures[1])
            for arc_match in eachmatch(ARC_COEFF_REGEX, expr)
                coeff = parse(Float64, arc_match.captures[1])
                coeff <= 0 && continue
                push!(rows, Dict(
                    "case_id" => meta.case_id,
                    "nusers" => meta.nusers,
                    "alpha" => meta.alpha,
                    "seed" => meta.seed,
                    "solver_config" => solver_config,
                    "user_id" => user_id,
                    "cut_index" => cut_idx,
                    "arc_from" => parse(Int, arc_match.captures[2]),
                    "arc_to" => parse(Int, arc_match.captures[3]),
                    "arc" => "($(arc_match.captures[2]), $(arc_match.captures[3]))",
                    "coeff" => coeff,
                    "rhs" => rhs,
                    "line" => line,
                ))
            end
        end
    end
    return DataFrame(rows)
end

function _summary_stats(values::Vector{Float64})
    vals = sort(values)
    return (
        n = length(vals),
        mean = mean(vals),
        median = median(vals),
        min = minimum(vals),
        q25 = quantile(vals, 0.25),
        q75 = quantile(vals, 0.75),
        max = maximum(vals),
    )
end

function _group_summary(df::DataFrame, group_cols::Vector{Symbol})
    out = Dict{String,Any}[]
    for sub in groupby(df, group_cols)
        vals = collect(Float64.(sub.coeff))
        stats = _summary_stats(vals)
        row = Dict{String,Any}()
        for col in group_cols
            row[string(col)] = sub[1, col]
        end
        row["n_coefficients"] = stats.n
        row["mean_coeff"] = stats.mean
        row["median_coeff"] = stats.median
        row["min_coeff"] = stats.min
        row["q25_coeff"] = stats.q25
        row["q75_coeff"] = stats.q75
        row["max_coeff"] = stats.max
        push!(out, row)
    end
    return DataFrame(out)
end

function _plot_style_map(labels::Vector{String})
    unique_labels = sort(unique(labels))
    return Dict(
        label => (
            color = POST_COLORS[mod1(i, length(POST_COLORS))],
            dash = POST_DASHES[mod1(i, length(POST_DASHES))]
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

function _draw_y_ticks(io, x0, x1, y0, y1, ymin, ymax; nticks=5, fmt=x -> round(x; digits=1))
    ymax == ymin && (ymax = ymin + 1)
    for i in 0:nticks
        v = ymin + (ymax - ymin) * i / nticks
        y = y1 - (y1 - y0) * i / nticks
        println(io, """<line x1="$x0" y1="$y" x2="$(x0-5)" y2="$y" stroke="#222" stroke-width="1"/>""")
        println(io, """<line x1="$x0" y1="$y" x2="$x1" y2="$y" stroke="#e5e5e5" stroke-width="1"/>""")
        println(io, """<text x="$(x0-8)" y="$(y+4)" text-anchor="end" class="tick">$(_svg_escape(fmt(v)))</text>""")
    end
end

function _draw_legend(io, styles::AbstractDict{String,<:NamedTuple}, labels::Vector{String}, x, y)
    for (idx, label) in enumerate(labels)
        style = styles[label]
        yy = y + 22 * (idx - 1)
        println(io, """<line x1="$x" y1="$yy" x2="$(x+28)" y2="$yy" stroke="$(style.color)" stroke-width="2.5" stroke-dasharray="$(style.dash)"/>""")
        println(io, """<text x="$(x+38)" y="$(yy+4)" class="legend">$(_svg_escape(label))</text>""")
    end
end

function _boxplot_stats(vals::Vector{Float64})
    s = sort(vals)
    return minimum(s), quantile(s, 0.25), median(s), quantile(s, 0.75), maximum(s)
end

function _bigm_boxplot(df::DataFrame, output_path::AbstractString)
    labels = sort(unique(String.(df.solver_config)))
    styles = _plot_style_map(labels)
    values = [collect(Float64.(df[df.solver_config .== label, :coeff])) for label in labels]
    all_vals = reduce(vcat, values)
    ymin, ymax = minimum(all_vals), maximum(all_vals)
    width, height = 980, 560
    x0, y0, x1, y1 = 90, 80, 860, 460
    boxw = 0.55 * (x1 - x0) / max(1, length(labels))

    open(output_path, "w") do io
        _svg_header(io, width, height)
        _draw_axes(io, x0, y0, x1, y1; xlabel="Solver / Parameter Combination", ylabel="Big-M Coefficient", title="Distribution of Big-M Coefficients")
        _draw_y_ticks(io, x0, x1, y0, y1, ymin, ymax; nticks=5, fmt=x -> string(round(Int, x)))
        for (i, label) in enumerate(labels)
            qmin, q1, q2, q3, qmax = _boxplot_stats(values[i])
            style = styles[label]
            cx = x0 + (i - 0.5) * (x1 - x0) / length(labels)
            ymap(v) = y1 - (y1 - y0) * (v - ymin) / max(1e-9, (ymax - ymin))
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
        _svg_footer(io)
    end
end

function _bigm_alpha_means(summary_by_solver_alpha::DataFrame, output_path::AbstractString)
    labels = sort(unique(String.(summary_by_solver_alpha.solver_config)))
    alphas = sort(unique(Float64.(summary_by_solver_alpha.alpha)))
    styles = _plot_style_map(labels)
    width, height = 1180, 560
    x0, y0, x1, y1 = 90, 80, 860, 460
    ymax = maximum(Float64.(summary_by_solver_alpha.mean_coeff)) * 1.1

    open(output_path, "w") do io
        _svg_header(io, width, height)
        _draw_axes(io, x0, y0, x1, y1; xlabel="Alpha", ylabel="Mean Big-M Coefficient", title="Mean Big-M Coefficient by Alpha")
        _draw_y_ticks(io, x0, x1, y0, y1, 0.0, ymax; nticks=5, fmt=x -> string(round(Int, x)))
        for (j, alpha) in enumerate(alphas)
            x = x0 + (x1 - x0) * (j - 1) / max(1, length(alphas) - 1)
            println(io, """<line x1="$x" y1="$y1" x2="$x" y2="$(y1+5)" stroke="#222" stroke-width="1"/>""")
            println(io, """<text x="$x" y="$(y1+20)" text-anchor="middle" class="tick">$alpha</text>""")
        end
        for label in labels
            style = styles[label]
            points = Tuple{Float64,Float64}[]
            for (j, alpha) in enumerate(alphas)
                sub = summary_by_solver_alpha[(summary_by_solver_alpha.solver_config .== label) .& (summary_by_solver_alpha.alpha .== alpha), :]
                isempty(sub) && continue
                x = x0 + (x1 - x0) * (j - 1) / max(1, length(alphas) - 1)
                y = y1 - (y1 - y0) * Float64(sub.mean_coeff[1]) / max(1e-9, ymax)
                push!(points, (x, y))
            end
            if length(points) >= 2
                path = join(["$(i == 1 ? "M" : "L") $(p[1]) $(p[2])" for (i, p) in enumerate(points)], " ")
                println(io, """<path d="$path" fill="none" stroke="$(style.color)" stroke-width="2.5" stroke-dasharray="$(style.dash)"/>""")
            end
            for (x, y) in points
                println(io, """<circle cx="$x" cy="$y" r="4" fill="$(style.color)" stroke="#111" stroke-width="1"/>""")
            end
        end
        _draw_legend(io, styles, labels, 900, 120)
        _svg_footer(io)
    end
end

function _top_arc_barplot(df::DataFrame, output_path::AbstractString; top_n::Int=20)
    top = first(sort(df, [:score, :high_occurrences, :mean_coeff], rev=true), min(top_n, nrow(df)))
    width, height = 1240, 620
    x0, y0, x1, y1 = 220, 80, 1140, 560
    ymax = maximum(Float64.(top.score)) * 1.1

    open(output_path, "w") do io
        _svg_header(io, width, height)
        _draw_axes(io, x0, y0, x1, y1; xlabel="Importance Score (sum coeff)", ylabel="Arc", title="Most Important Arcs for U50, BlCLag+BlC ws=on")
        _draw_y_ticks(io, x0, x1, y0, y1, 0.0, ymax; nticks=5, fmt=x -> string(round(Int, x)))
        barh = 0.72 * (y1 - y0) / max(1, nrow(top))
        for (i, row) in enumerate(eachrow(top))
            cy = y0 + (i - 0.5) * (y1 - y0) / nrow(top)
            w = (x1 - x0) * Float64(row.score) / max(1e-9, ymax)
            println(io, """<rect x="$x0" y="$(cy - barh/2)" width="$w" height="$barh" fill="#1f77b4" fill-opacity="0.65" stroke="#111" stroke-width="1.2"/>""")
            println(io, """<text x="$(x0-10)" y="$(cy+4)" text-anchor="end" class="tick">$(_svg_escape(string(row.arc)))</text>""")
            println(io, """<text x="$(x0 + w + 8)" y="$(cy+4)" class="tick">score=$(round(Int, row.score)), users=$(Int(row.users_count))</text>""")
        end
        _svg_footer(io)
    end
end

function _performance_profile_audit(results_csv::AbstractString, output_dir::AbstractString)
    df = CSV.read(String(results_csv), DataFrame)
    labels = sort(unique(String.(df.solver_config)))
    case_ids = sort(unique(String.(df.case_id)))
    rows = Dict{String,Any}[]
    best_counts = Dict(label => 0 for label in labels)
    for case_id in case_ids
        sub = df[df.case_id .== case_id, :]
        effective = Dict{String,Float64}()
        for row in eachrow(sub)
            runtime = Float64(row.runtime)
            limit = Float64(row.runtime_limit)
            effective[String(row.solver_config)] = lowercase(String(row.status)) == "optimal" ? runtime : 2.0 * limit
        end
        best = minimum(values(effective))
        for label in labels
            ratio = effective[label] / best
            push!(rows, Dict("case_id" => case_id, "solver_config" => label, "effective_runtime" => effective[label], "ratio" => ratio))
            abs(ratio - 1.0) <= 1e-9 && (best_counts[label] += 1)
        end
    end
    audit_df = DataFrame(rows)
    CSV.write(joinpath(output_dir, "performance_profile_audit.csv"), audit_df)
    summary_rows = Dict{String,Any}[]
    for label in labels
        vals = sort(Float64.(audit_df[audit_df.solver_config .== label, :ratio]))
        push!(summary_rows, Dict(
            "solver_config" => label,
            "best_case_count" => best_counts[label],
            "median_ratio" => median(vals),
            "max_ratio" => maximum(vals),
            "ratio_le_2_count" => count(<=(2.0), vals),
        ))
    end
    summary_df = DataFrame(summary_rows)
    CSV.write(joinpath(output_dir, "performance_profile_audit_summary.csv"), summary_df)
    return summary_df
end

function _important_arcs_u50_ws1(raw_df::DataFrame)
    sub = raw_df[(raw_df.nusers .== 50) .& (raw_df.solver_config .== "BlCLag+BlC ws=on"), :]
    isempty(sub) && return DataFrame(), DataFrame()

    user_thresholds = Dict{Int,Float64}()
    for g in groupby(sub, :user_id)
        user_thresholds[Int(g.user_id[1])] = quantile(Float64.(g.coeff), 0.9)
    end

    per_user_rows = Dict{String,Any}[]
    for g in groupby(sub, [:user_id, :arc])
        vals = Float64.(g.coeff)
        user_id = Int(g.user_id[1])
        high_thr = user_thresholds[user_id]
        high_occ = count(>=(high_thr), vals)
        push!(per_user_rows, Dict(
            "user_id" => user_id,
            "arc" => String(g.arc[1]),
            "count" => length(vals),
            "high_occurrences" => high_occ,
            "mean_coeff" => mean(vals),
            "max_coeff" => maximum(vals),
            "sum_coeff" => sum(vals),
            "alphas_present" => join(sort(unique(Float64.(g.alpha))), ", "),
            "score" => sum(vals),
        ))
    end
    per_user_df = DataFrame(per_user_rows)

    overall_rows = Dict{String,Any}[]
    for g in groupby(per_user_df, :arc)
        push!(overall_rows, Dict(
            "arc" => String(g.arc[1]),
            "users_count" => nrow(g),
            "high_occurrences" => sum(Int.(g.high_occurrences)),
            "mean_coeff" => mean(Float64.(g.mean_coeff)),
            "max_coeff" => maximum(Float64.(g.max_coeff)),
            "score" => sum(Float64.(g.score)),
        ))
    end
    overall_df = DataFrame(overall_rows)
    return per_user_df, overall_df
end

function run_hndp_cut_postprocess(benchmark_root::AbstractString; output_dir::AbstractString=joinpath(String(benchmark_root), "postprocess"))
    mkpath(output_dir)
    cut_files = _collect_mastercut_files(benchmark_root)
    isempty(cut_files) && error("No mastercuts_collection.txt files found under $(benchmark_root)")

    raw_df = _parse_mastercuts(cut_files)
    isempty(raw_df) && error("No positive cut coefficients parsed from the master cut collections under $(benchmark_root)")
    CSV.write(joinpath(output_dir, "bigm_coefficients_raw.csv"), raw_df)

    summary_by_solver = _group_summary(raw_df, [:solver_config])
    summary_by_solver_alpha = _group_summary(raw_df, [:solver_config, :alpha])
    CSV.write(joinpath(output_dir, "bigm_summary_by_solver.csv"), summary_by_solver)
    CSV.write(joinpath(output_dir, "bigm_summary_by_solver_alpha.csv"), summary_by_solver_alpha)

    _bigm_boxplot(raw_df, joinpath(output_dir, "bigm_boxplot.svg"))
    _bigm_alpha_means(summary_by_solver_alpha, joinpath(output_dir, "bigm_mean_by_alpha.svg"))

    per_user_df, overall_df = _important_arcs_u50_ws1(raw_df)
    if !isempty(per_user_df)
        CSV.write(joinpath(output_dir, "u50_ws1_important_arcs_per_user.csv"), sort(per_user_df, [:user_id, :score, :high_occurrences], rev=[false, true, true]))
        CSV.write(joinpath(output_dir, "u50_ws1_important_arcs_overall.csv"), sort(overall_df, [:score, :high_occurrences], rev=true))
        _top_arc_barplot(overall_df, joinpath(output_dir, "u50_ws1_important_arcs_top20.svg"))
    end

    results_csv = joinpath(String(benchmark_root), "results", "benchmark_results.csv")
    audit_summary = isfile(results_csv) ? _performance_profile_audit(results_csv, output_dir) : DataFrame()

    return (
        raw_df = raw_df,
        summary_by_solver = summary_by_solver,
        summary_by_solver_alpha = summary_by_solver_alpha,
        important_arcs_per_user = per_user_df,
        important_arcs_overall = overall_df,
        performance_profile_audit = audit_summary,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && error("Usage: julia --project=. examples/HNDP/hndp_cut_postprocess.jl <benchmark_root> [output_dir]")
    benchmark_root = ARGS[1]
    output_dir = length(ARGS) >= 2 ? ARGS[2] : joinpath(benchmark_root, "postprocess")
    run_hndp_cut_postprocess(benchmark_root; output_dir=output_dir)
    println("Wrote HNDP cut postprocessing results to " * output_dir)
end
