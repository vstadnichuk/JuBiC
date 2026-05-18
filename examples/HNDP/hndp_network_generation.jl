using JSON
using Random
using Graphs

include("hndp_instances.jl")

"""
    HNDPGeneratedNetwork

Container for a generated HNDP network instance together with the metadata
describing how it was produced. The `instance` field stores the network and user
data only; solver-specific JuBiC models are intentionally built later.
"""
struct HNDPGeneratedNetwork
    name::String
    instance::HNDPwC
    metadata::Dict{String,Any}
end

"""
    load_hndp_network_generation_config(json_path)

Load the JSON configuration used by the new HNDP network generator.
"""
function load_hndp_network_generation_config(json_path::String)
    return JSON.parsefile(json_path)
end

"""
    generate_hndp_networks(json_path)
    generate_hndp_networks(config)

Generate a batch of HNDP network instances from a JSON path or parsed
configuration dictionary.

The current implementation supports two logical instance types:
- `"constrained_shortest_path"`: a single-layer HNDP network
- `"competition"`: a layered HNDP network where only the competitor layer is
  interdicted or constructed

For each entry in `instances`, the generator builds the Cartesian product of the
declared parameter values and the selected parameter seeds. Different topology
families are represented by different entries in the JSON file, which keeps the
pipeline easy to extend later with additional random topology generators.
"""
function generate_hndp_networks(json_path::String)
    return generate_hndp_networks(load_hndp_network_generation_config(json_path))
end

function generate_hndp_networks(config::Dict{String,Any})
    specs = get(config, "instances", nothing)
    specs === nothing && throw(ArgumentError("The network generation config requires an 'instances' field."))

    default_parameter_seeds = _int_vector(get(config, "parameter_seeds", [0]), "parameter_seeds")
    generated = HNDPGeneratedNetwork[]

    for raw_spec in specs
        spec = Dict{String,Any}(raw_spec)
        append!(generated, _expand_hndp_generation_spec!(spec, default_parameter_seeds))
    end

    if haskey(config, "instance_output")
        write_hndp_generated_networks(generated, Dict{String,Any}(config["instance_output"]))
    end

    return generated
end

"""
    visit_hndp_networks(json_path, visitor)
    visit_hndp_networks(config, visitor)

Stream generated HNDP networks to `visitor` one by one instead of returning all
instances at once. This is useful for large experiment pipelines where keeping
all generated instances in memory would be wasteful.
"""
function visit_hndp_networks(json_path::String, visitor::Function)
    return visit_hndp_networks(load_hndp_network_generation_config(json_path), visitor)
end

function visit_hndp_networks(config::Dict{String,Any}, visitor::Function)
    specs = get(config, "instances", nothing)
    specs === nothing && throw(ArgumentError("The network generation config requires an 'instances' field."))

    default_parameter_seeds = _int_vector(get(config, "parameter_seeds", [0]), "parameter_seeds")
    output_cfg = haskey(config, "instance_output") ? Dict{String,Any}(config["instance_output"]) : nothing

    for raw_spec in specs
        spec = Dict{String,Any}(raw_spec)
        _visit_hndp_generation_spec!(spec, default_parameter_seeds, output_cfg, visitor)
    end
    return nothing
end

function _get_optional_float_vector(cfg::Dict{String,Any}, key::String)
    haskey(cfg, key) || return nothing
    return _float_vector(cfg[key], key)
end

function _visit_hndp_generation_spec!(
    spec::Dict{String,Any},
    default_parameter_seeds::Vector{Int},
    output_cfg,
    visitor::Function,
)
    for generated in _expand_hndp_generation_spec!(spec, default_parameter_seeds)
        if !isnothing(output_cfg)
            write_hndp_generated_networks([generated], output_cfg)
        end
        visitor(generated)
    end
    return nothing
end

function _expand_hndp_generation_spec!(
    spec::Dict{String,Any},
    default_parameter_seeds::Vector{Int},
)
    instance_type = _required_string(spec, "instance_type")
    base_name = get(spec, "name", instance_type)
    parameter_seeds = _int_vector(get(spec, "parameter_seeds", default_parameter_seeds), "parameter_seeds")
    topology_ids = _string_vector(get(spec, "topologies", [get(spec, "topology_family", _default_topology_family(instance_type))]), "topologies")
    user_counts = _int_vector(get(spec, "nusers", [1]), "nusers")
    length_modes = _bool_vector(get(spec, "length_constrained", [true]), "length_constrained")

    generated_here = HNDPGeneratedNetwork[]

    if instance_type == "constrained_shortest_path"
        alpha_values = _get_optional_float_vector(spec, "alpha")
        length_slack_values = isnothing(alpha_values) ? _float_vector(get(spec, "length_slack", [1.0]), "length_slack") : alpha_values
        two_stage_modes = _bool_vector(get(spec, "two_stage", [false]), "two_stage")
        user_parameter_modes = _string_vector(get(spec, "user_parameter_mode", ["shared"]), "user_parameter_mode")

        for topology_id in topology_ids
            for nusers in user_counts
                for parameter_seed in parameter_seeds
                    for constrained in length_modes
                        for length_slack in length_slack_values
                            for two_stage in two_stage_modes
                                for user_parameter_mode in user_parameter_modes
                                    !constrained && length_slack != first(length_slack_values) && continue
                                    generated = _build_csp_generated_network(
                                        base_name,
                                        topology_id,
                                        nusers,
                                        parameter_seed,
                                        constrained,
                                        length_slack,
                                        two_stage,
                                        user_parameter_mode,
                                        spec,
                                    )
                                    push!(generated_here, generated)
                                end
                            end
                        end
                    end
                end
            end
        end
    elseif instance_type == "competition"
        alpha_values = _get_optional_float_vector(spec, "alpha")
        beta_values = _get_optional_float_vector(spec, "beta")
        length_slack_values = isnothing(alpha_values) ? _float_vector(get(spec, "length_slack", [1.0]), "length_slack") : alpha_values
        competitor_cost_factor_values = isnothing(beta_values) ? _float_vector(get(spec, "competitor_cost_factor", [0.8]), "competitor_cost_factor") : beta_values
        user_parameter_modes = _string_vector(get(spec, "user_parameter_mode", ["shared"]), "user_parameter_mode")
        decision_arc_counts = _int_vector(get(spec, "decision_arc_count", [-1]), "decision_arc_count")

        for topology_id in topology_ids
            for nusers in user_counts
                for parameter_seed in parameter_seeds
                    for constrained in length_modes
                        for length_slack in length_slack_values
                            for competitor_cost_factor in competitor_cost_factor_values
                                for decision_arc_count in decision_arc_counts
                                    for user_parameter_mode in user_parameter_modes
                                        !constrained && length_slack != first(length_slack_values) && continue
                                        generated = _build_competition_generated_network(
                                            base_name,
                                            topology_id,
                                            nusers,
                                            parameter_seed,
                                            constrained,
                                            length_slack,
                                            competitor_cost_factor,
                                            user_parameter_mode,
                                            spec;
                                            decision_arc_count=decision_arc_count,
                                        )
                                        push!(generated_here, generated)
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    else
        throw(ArgumentError("Unsupported HNDP instance_type '$instance_type'."))
    end

    return generated_here
end

function _build_csp_generated_network(
    base_name::String,
    topology_id::String,
    nusers::Int,
    parameter_seed::Int,
    constrained::Bool,
    length_slack::Float64,
    two_stage::Bool,
    user_parameter_mode::String,
    spec::Dict{String,Any},
)
    _validate_user_parameter_mode(user_parameter_mode)

    graph = _load_named_base_graph(topology_id)
    max_cost = Int(get(spec, "max_cost", 100))
    max_risk = Int(get(spec, "max_risk", 100))
    max_weight = Int(get(spec, "max_weight", 100))
    construction_cost = Int(get(spec, "construction_cost", 0))

    users, minweights = _build_single_layer_users(
        graph,
        nusers,
        parameter_seed,
        constrained,
        length_slack,
        max_cost,
        max_risk,
        max_weight,
        two_stage,
        user_parameter_mode,
    )
    edgeA = [(src(e), dst(e)) for e in edges(graph)]
    edge_price = Dict((src(e), dst(e)) => rand(MersenneTwister(parameter_seed + 10_000), 0:construction_cost) for e in edges(graph))
    instance = HNDPwC(graph, users, edgeA, edge_price, minweights)

    metadata = Dict{String,Any}(
        "name" => _generated_name(base_name, topology_id, nusers, parameter_seed, constrained; length_slack=length_slack, two_stage=two_stage, user_parameter_mode=user_parameter_mode),
        "instance_type" => "constrained_shortest_path",
        "objective_style" => "weighted",
        "topology_family" => topology_id,
        "nusers" => nusers,
        "parameter_seed" => parameter_seed,
        "length_constrained" => constrained,
        "alpha" => constrained ? length_slack : nothing,
        "length_slack" => constrained ? length_slack : nothing,
        "two_stage" => two_stage,
        "user_parameter_mode" => user_parameter_mode,
        "construction_cost" => construction_cost,
        "max_cost" => max_cost,
        "max_risk" => max_risk,
        "max_weight" => max_weight,
        "nnodes" => nv(graph),
        "narcs" => ne(graph),
        "decision_arcs" => length(edgeA),
    )

    return HNDPGeneratedNetwork(metadata["name"], instance, metadata)
end

"""
    write_hndp_generated_networks(generated, output_cfg)

Write generated HNDP networks to disk for debugging and inspection. Each
instance is written into its own folder. The main export is a GEXF file that
embeds node and edge attributes, including user-specific arc data. A companion
`users.json` file stores the realized OD pairs and path-length bounds.
"""
function write_hndp_generated_networks(
    generated::Vector{HNDPGeneratedNetwork},
    output_cfg::Dict{String,Any},
)
    base_folder = String(get(output_cfg, "folder", "examples/HNDP/generated_instances"))
    graph_format = lowercase(String(get(output_cfg, "graph_format", "gexf")))
    graph_format == "gexf" || throw(ArgumentError("Currently only the 'gexf' graph_format is supported."))

    mkpath(base_folder)
    for generated_network in generated
        _write_generated_network(generated_network, base_folder, graph_format)
    end
end

function _write_generated_network(
    generated_network::HNDPGeneratedNetwork,
    base_folder::String,
    graph_format::String,
)
    instance_folder = joinpath(base_folder, generated_network.name)
    mkpath(instance_folder)

    if graph_format == "gexf"
        _write_hndp_gexf(joinpath(instance_folder, "instance.gexf"), generated_network)
    end
    _write_hndp_users_json(joinpath(instance_folder, "users.json"), generated_network)
end

function _build_competition_generated_network(
    base_name::String,
    topology_id::String,
    nusers::Int,
    parameter_seed::Int,
    constrained::Bool,
    length_slack::Float64,
    competitor_cost_factor::Float64,
    user_parameter_mode::String,
    spec::Dict{String,Any},
    ;
    decision_arc_count::Int=-1,
)
    _validate_user_parameter_mode(user_parameter_mode)

    max_cost = Int(get(spec, "max_cost", 100))
    max_risk = Int(get(spec, "max_risk", 100))
    max_weight = Int(get(spec, "max_weight", 100))
    construction_cost = Int(get(spec, "construction_cost", 0))

    if topology_id == "layered_sioux_falls"
        base_graph = _load_sioux_falls_graph()
        graph, base_nodes, competitor_nodes, decision_arcs = _build_layered_competition_graph(base_graph)

        users, minweights = _build_two_layer_users(
            graph,
            base_nodes,
            nusers,
            parameter_seed,
            constrained,
            length_slack,
            max_cost,
            max_risk,
            max_weight,
            base_graph,
            competitor_nodes,
            decision_arcs,
            competitor_cost_factor,
            user_parameter_mode,
        )
        edge_price = Dict((src(e), dst(e)) => rand(MersenneTwister(parameter_seed + 20_000), 0:construction_cost) for e in edges(graph))
        instance = HNDPwC(graph, users, decision_arcs, edge_price, minweights)

        metadata = Dict{String,Any}(
            "name" => _generated_name(base_name, topology_id, nusers, parameter_seed, constrained; length_slack=length_slack, competitor_cost_factor=competitor_cost_factor, user_parameter_mode=user_parameter_mode, decision_arc_count=length(decision_arcs)),
            "instance_type" => "competition",
            "topology_family" => topology_id,
            "nusers" => nusers,
            "parameter_seed" => parameter_seed,
            "length_constrained" => constrained,
            "alpha" => constrained ? length_slack : nothing,
            "length_slack" => constrained ? length_slack : nothing,
            "beta" => competitor_cost_factor,
            "competitor_cost_factor" => competitor_cost_factor,
            "user_parameter_mode" => user_parameter_mode,
            "construction_cost" => construction_cost,
            "max_cost" => max_cost,
            "max_risk" => max_risk,
            "max_weight" => max_weight,
            "nnodes" => nv(graph),
            "narcs" => ne(graph),
            "decision_arcs" => length(decision_arcs),
        )

        return HNDPGeneratedNetwork(metadata["name"], instance, metadata)
    end

    if startswith(topology_id, "layered_")
        base_topology_id = topology_id[(findfirst('_', topology_id) + 1):end]
        base_graph = _load_named_base_graph(base_topology_id)
        graph, base_nodes, competitor_nodes, decision_arcs = _build_two_layer_competition_graph(base_graph)

        users, minweights = _build_layered_users(
            graph,
            base_nodes,
            nusers,
            parameter_seed,
            constrained,
            length_slack,
            max_cost,
            max_risk,
            max_weight,
            base_graph,
            competitor_nodes,
            decision_arcs,
            competitor_cost_factor,
            user_parameter_mode,
        )
        edge_price = Dict((src(e), dst(e)) => rand(MersenneTwister(parameter_seed + 20_000), 0:construction_cost) for e in edges(graph))
        instance = HNDPwC(graph, users, decision_arcs, edge_price, minweights)

        metadata = Dict{String,Any}(
            "name" => _generated_name(base_name, topology_id, nusers, parameter_seed, constrained; length_slack=length_slack, competitor_cost_factor=competitor_cost_factor, user_parameter_mode=user_parameter_mode, decision_arc_count=length(decision_arcs)),
            "instance_type" => "competition",
            "topology_family" => topology_id,
            "nusers" => nusers,
            "parameter_seed" => parameter_seed,
            "length_constrained" => constrained,
            "alpha" => constrained ? length_slack : nothing,
            "length_slack" => constrained ? length_slack : nothing,
            "beta" => competitor_cost_factor,
            "competitor_cost_factor" => competitor_cost_factor,
            "user_parameter_mode" => user_parameter_mode,
            "construction_cost" => construction_cost,
            "max_cost" => max_cost,
            "max_risk" => max_risk,
            "max_weight" => max_weight,
            "nnodes" => nv(graph),
            "narcs" => ne(graph),
            "decision_arcs" => length(decision_arcs),
        )

        return HNDPGeneratedNetwork(metadata["name"], instance, metadata)
    end

    base_graph = _load_named_base_graph(topology_id)
    graph = base_graph
    decision_arcs = _sample_decision_arcs(graph, parameter_seed, decision_arc_count)

    users, minweights = _build_single_layer_competition_users(
        graph,
        nusers,
        parameter_seed,
        constrained,
        length_slack,
        max_cost,
        max_risk,
        max_weight,
        decision_arcs,
        competitor_cost_factor,
        user_parameter_mode,
    )
    edge_price = Dict((src(e), dst(e)) => rand(MersenneTwister(parameter_seed + 20_000), 0:construction_cost) for e in edges(graph))
    instance = HNDPwC(graph, users, decision_arcs, edge_price, minweights)

    metadata = Dict{String,Any}(
        "name" => _generated_name(base_name, topology_id, nusers, parameter_seed, constrained; length_slack=length_slack, competitor_cost_factor=competitor_cost_factor, user_parameter_mode=user_parameter_mode, decision_arc_count=length(decision_arcs)),
        "instance_type" => "competition",
        "topology_family" => topology_id,
        "nusers" => nusers,
        "parameter_seed" => parameter_seed,
        "length_constrained" => constrained,
        "alpha" => constrained ? length_slack : nothing,
        "length_slack" => constrained ? length_slack : nothing,
        "beta" => competitor_cost_factor,
        "competitor_cost_factor" => competitor_cost_factor,
        "user_parameter_mode" => user_parameter_mode,
        "construction_cost" => construction_cost,
        "max_cost" => max_cost,
        "max_risk" => max_risk,
        "max_weight" => max_weight,
        "nnodes" => nv(graph),
        "narcs" => ne(graph),
        "decision_arcs" => length(decision_arcs),
    )

    return HNDPGeneratedNetwork(metadata["name"], instance, metadata)
end

function _required_string(cfg::Dict{String,Any}, key::String)
    haskey(cfg, key) || throw(ArgumentError("Missing required config key '$key'."))
    value = cfg[key]
    value isa String || throw(ArgumentError("Config entry '$key' must be a string."))
    return value
end

function _int_vector(value, key::String)
    values = value isa AbstractVector ? value : [value]
    return [Int(v) for v in values]
end

function _float_vector(value, key::String)
    values = value isa AbstractVector ? value : [value]
    return [Float64(v) for v in values]
end

function _bool_vector(value, key::String)
    values = value isa AbstractVector ? value : [value]
    return [Bool(v) for v in values]
end

function _string_vector(value, key::String)
    values = value isa AbstractVector ? value : [value]
    return [String(v) for v in values]
end

function _default_topology_family(instance_type::String)
    if instance_type == "constrained_shortest_path"
        return "sioux_falls"
    elseif instance_type == "competition"
        return "sioux_falls"
    end
    throw(ArgumentError("Unsupported HNDP instance_type '$instance_type'."))
end

function _load_sioux_falls_graph()
    file_path = "examples/data/SF_DNDP_10_base.txt"
    edge_list = parse_Sioux_Falls_csv(file_path)
    return create_graph_from_edges_SiouxFalls(edge_list)
end

function _load_named_base_graph(topology_id::String)
    if topology_id == "sioux_falls"
        return _load_sioux_falls_graph()
    elseif topology_id == "anaheim"
        return _load_tntp_graph("examples/data/Anaheim_net.tntp")
    elseif topology_id == "berlin_mitte_center"
        return _load_tntp_graph("examples/data/berlin-mitte-center_net.tntp")
    elseif topology_id == "ema"
        return _load_tntp_graph("examples/data/EMA_net.tntp")
    end
    throw(ArgumentError("Unsupported HNDP topology '$topology_id'."))
end

function _load_tntp_graph(file_path::String)
    edge_list, nnodes = parse_tntp_edges(file_path)
    graph = DiGraph(nnodes)
    for (i, j) in edge_list
        add_edge!(graph, i, j)
    end
    return graph
end

function parse_tntp_edges(file_path::String)
    lines = readlines(file_path)
    node_line = findfirst(x -> occursin("<NUMBER OF NODES>", x), lines)
    isnothing(node_line) && throw(ArgumentError("Could not find <NUMBER OF NODES> in TNTP file $file_path"))
    nnodes = parse(Int, split(strip(lines[node_line]))[end])

    start_idx = findfirst(x -> occursin("<END OF METADATA>", x), lines)
    isnothing(start_idx) && throw(ArgumentError("Could not find <END OF METADATA> in TNTP file $file_path"))

    edges = Tuple{Int,Int}[]
    for line in lines[(start_idx + 1):end]
        stripped = strip(line)
        isempty(stripped) && continue
        startswith(stripped, "~") && continue
        tokens = split(stripped, r"\s+")
        length(tokens) < 2 && continue
        if occursin(r"^\d+$", tokens[1]) && occursin(r"^\d+$", tokens[2])
            push!(edges, (parse(Int, tokens[1]), parse(Int, tokens[2])))
        end
    end
    return edges, nnodes
end

function _validate_user_parameter_mode(user_parameter_mode::String)
    if user_parameter_mode != "shared" && user_parameter_mode != "per_user"
        throw(ArgumentError("Unsupported user_parameter_mode '$user_parameter_mode'. Supported values are 'shared' and 'per_user'."))
    end
end

function _random_arc_matrices(
    graph::DiGraph,
    seed::Int;
    max_cost::Int,
    max_risk::Int,
    max_weight::Int,
    negative_risk::Bool,
)
    rng = MersenneTwister(seed)
    n = nv(graph)
    rcost = _random_nonnegative_matrix(rng, n, max_cost)
    rrisk = _random_nonnegative_matrix(rng, n, max_risk)
    rweight = _random_nonnegative_matrix(rng, n, max_weight)

    negative_risk && (rrisk .*= -1)

    for (i, j) in Iterators.product(1:n, 1:n)
        if !has_edge(graph, i, j)
            rcost[i, j] = 0
            rrisk[i, j] = 0
            rweight[i, j] = 0
        end
    end

    return rcost, rrisk, rweight
end

function _random_nonnegative_matrix(rng, n::Int, max_value::Int)
    max_value >= 0 || throw(ArgumentError("Maximum matrix values must be nonnegative."))
    if max_value == 0
        return zeros(Int, n, n)
    end
    return rand(rng, 1:max_value, n, n)
end

function _build_single_layer_users(
    graph::DiGraph,
    nusers::Int,
    seed::Int,
    constrained::Bool,
    length_slack,
    max_cost::Int,
    max_risk::Int,
    max_weight::Int,
    two_stage::Bool,
    user_parameter_mode::String,
)
    rng = MersenneTwister(seed)
    users = User[]
    feasible_targets = _feasible_od_targets(graph)

    if user_parameter_mode == "shared"
        rcost, rrisk, rweight = _random_arc_matrices(
            graph,
            seed;
            max_cost=max_cost,
            max_risk=max_risk,
            max_weight=max_weight,
            negative_risk=false,
        )
        if two_stage
            rrisk = copy(rcost)
        end
        user_weight = constrained ? rweight : nothing
        minweights = constrained ? floyd_warshall_shortest_paths(graph, rweight) : nothing
        for user_id in 1:nusers
            origin, destination = _sample_feasible_od_pair(rng, feasible_targets)
            bound = _compute_weight_bound(minweights, origin, destination, length_slack, nv(graph) * max_weight)
            push!(users, User("U$(user_id)", origin, destination, rrisk, rcost, user_weight, bound))
        end
        _validate_users_reachable!(graph, users)
        return users, minweights
    end

    for user_id in 1:nusers
        user_seed = seed + 1000 * user_id
        rcost, rrisk, rweight = _random_arc_matrices(
            graph,
            user_seed;
            max_cost=max_cost,
            max_risk=max_risk,
            max_weight=max_weight,
            negative_risk=false,
        )
        if two_stage
            rrisk = copy(rcost)
        end
        user_weight = constrained ? rweight : nothing
        user_minweights = constrained ? floyd_warshall_shortest_paths(graph, rweight) : nothing
        origin, destination = _sample_feasible_od_pair(rng, feasible_targets)
        bound = _compute_weight_bound(user_minweights, origin, destination, length_slack, nv(graph) * max_weight)
        push!(users, User("U$(user_id)", origin, destination, rrisk, rcost, user_weight, bound))
    end
    _validate_users_reachable!(graph, users)
    return users, nothing
end

function _sample_decision_arcs(graph::DiGraph, seed::Int, decision_arc_count::Int)
    all_arcs = Tuple{Int,Int}[(src(e), dst(e)) for e in edges(graph)]
    if decision_arc_count < 0 || decision_arc_count >= length(all_arcs)
        return all_arcs
    elseif decision_arc_count == 0
        return Tuple{Int,Int}[]
    end
    rng = MersenneTwister(seed + 50_000)
    perm = randperm(rng, length(all_arcs))
    return all_arcs[perm[1:decision_arc_count]]
end

function _reachable_destinations(graph::DiGraph, origin::Int)
    visited = falses(nv(graph))
    queue = [origin]
    visited[origin] = true
    head = 1
    while head <= length(queue)
        node = queue[head]
        head += 1
        for neigh in outneighbors(graph, node)
            if !visited[neigh]
                visited[neigh] = true
                push!(queue, neigh)
            end
        end
    end
    return [node for node in 1:nv(graph) if node != origin && visited[node]]
end

function _feasible_od_targets(graph::DiGraph)
    feasible = Dict{Int,Vector{Int}}()
    for origin in 1:nv(graph)
        dests = _reachable_destinations(graph, origin)
        isempty(dests) || (feasible[origin] = dests)
    end
    isempty(feasible) && throw(ArgumentError("The passed directed graph contains no feasible origin-destination pair with a directed path."))
    return feasible
end

function _sample_feasible_od_pair(rng, feasible_targets::Dict{Int,Vector{Int}})
    origins = collect(keys(feasible_targets))
    origin = rand(rng, origins)
    destination = rand(rng, feasible_targets[origin])
    return origin, destination
end

function _validate_users_reachable!(graph::DiGraph, users::Vector{User})
    feasible_targets = _feasible_od_targets(graph)
    for user in users
        if !haskey(feasible_targets, user.origin) || !(user.destination in feasible_targets[user.origin])
            throw(ArgumentError("Generated invalid OD pair for $(user.uname): no directed path from $(user.origin) to $(user.destination)."))
        end
    end
end

function _build_single_layer_competition_users(
    graph::DiGraph,
    nusers::Int,
    seed::Int,
    constrained::Bool,
    length_slack,
    max_cost::Int,
    max_risk::Int,
    max_weight::Int,
    decision_arcs,
    competitor_cost_factor::Float64,
    user_parameter_mode::String,
)
    rng = MersenneTwister(seed)
    users = User[]
    decision_arc_set = Set(decision_arcs)
    feasible_targets = _feasible_od_targets(graph)

    if user_parameter_mode == "shared"
        rcost, rrisk, rweight = _random_arc_matrices(
            graph,
            seed;
            max_cost=max_cost,
            max_risk=max_risk,
            max_weight=max_weight,
            negative_risk=true,
        )
        _apply_single_layer_competition_arc_rules!(graph, rcost, rrisk, rweight, decision_arc_set, competitor_cost_factor)
        user_weight = constrained ? rweight : nothing
        minweights = constrained ? floyd_warshall_shortest_paths(graph, rweight) : nothing
        for user_id in 1:nusers
            origin, destination = _sample_feasible_od_pair(rng, feasible_targets)
            bound = _compute_weight_bound(minweights, origin, destination, length_slack, nv(graph) * max_weight)
            push!(users, User("U$(user_id)", origin, destination, rrisk, rcost, user_weight, bound))
        end
        _validate_users_reachable!(graph, users)
        return users, minweights
    end

    for user_id in 1:nusers
        user_seed = seed + 1000 * user_id
        rcost, rrisk, rweight = _random_arc_matrices(
            graph,
            user_seed;
            max_cost=max_cost,
            max_risk=max_risk,
            max_weight=max_weight,
            negative_risk=true,
        )
        _apply_single_layer_competition_arc_rules!(graph, rcost, rrisk, rweight, decision_arc_set, competitor_cost_factor)
        user_weight = constrained ? rweight : nothing
        user_minweights = constrained ? floyd_warshall_shortest_paths(graph, rweight) : nothing
        origin, destination = _sample_feasible_od_pair(rng, feasible_targets)
        bound = _compute_weight_bound(user_minweights, origin, destination, length_slack, nv(graph) * max_weight)
        push!(users, User("U$(user_id)", origin, destination, rrisk, rcost, user_weight, bound))
    end
    _validate_users_reachable!(graph, users)
    return users, nothing
end

function _apply_single_layer_competition_arc_rules!(
    graph::DiGraph,
    rcost,
    rrisk,
    rweight,
    decision_arc_set,
    competitor_cost_factor::Float64,
)
    for (i, j) in Iterators.product(1:nv(graph), 1:nv(graph))
        has_edge(graph, i, j) || continue
        if (i, j) in decision_arc_set
            rcost[i, j] = round(Int, competitor_cost_factor * rcost[i, j])
        else
            rrisk[i, j] = 0
        end
    end
end

function _build_layered_users(
    graph::DiGraph,
    base_nodes,
    nusers::Int,
    seed::Int,
    constrained::Bool,
    length_slack,
    max_cost::Int,
    max_risk::Int,
    max_weight::Int,
    base_graph::DiGraph,
    competitor_nodes,
    decision_arcs,
    competitor_cost_factor::Float64,
    user_parameter_mode::String,
)
    rng = MersenneTwister(seed)
    users = User[]
    max_possible_weight = 4 * length(base_nodes) * max_weight
    feasible_targets = _feasible_od_targets(base_graph)

    if user_parameter_mode == "shared"
        rcost, rrisk, rweight = _random_arc_matrices(
            graph,
            seed;
            max_cost=max_cost,
            max_risk=max_risk,
            max_weight=max_weight,
            negative_risk=true,
        )
        _apply_competition_arc_rules!(
            graph,
            base_graph,
            base_nodes,
            competitor_nodes,
            decision_arcs,
            rcost,
            rrisk,
            rweight,
            competitor_cost_factor,
        )
        user_weight = constrained ? rweight : nothing
        minweights = constrained ? floyd_warshall_shortest_paths(graph, rweight) : nothing
        for user_id in 1:nusers
            base_origin, base_destination = _sample_feasible_od_pair(rng, feasible_targets)
            origin = base_nodes[base_origin]
            destination = base_nodes[base_destination]
            bound = _compute_weight_bound(minweights, origin, destination, length_slack, max_possible_weight)
            push!(users, User("U$(user_id)", origin, destination, rrisk, rcost, user_weight, bound))
        end
        return users, minweights
    end

    for user_id in 1:nusers
        user_seed = seed + 1000 * user_id
        rcost, rrisk, rweight = _random_arc_matrices(
            graph,
            user_seed;
            max_cost=max_cost,
            max_risk=max_risk,
            max_weight=max_weight,
            negative_risk=true,
        )
        _apply_competition_arc_rules!(
            graph,
            base_graph,
            base_nodes,
            competitor_nodes,
            decision_arcs,
            rcost,
            rrisk,
            rweight,
            competitor_cost_factor,
        )
        user_weight = constrained ? rweight : nothing
        user_minweights = constrained ? floyd_warshall_shortest_paths(graph, rweight) : nothing
        base_origin, base_destination = _sample_feasible_od_pair(rng, feasible_targets)
        origin = base_nodes[base_origin]
        destination = base_nodes[base_destination]
        bound = _compute_weight_bound(user_minweights, origin, destination, length_slack, max_possible_weight)
        push!(users, User("U$(user_id)", origin, destination, rrisk, rcost, user_weight, bound))
    end
    return users, nothing
end

function _build_two_layer_users(
    graph::DiGraph,
    base_nodes,
    nusers::Int,
    seed::Int,
    constrained::Bool,
    length_slack,
    max_cost::Int,
    max_risk::Int,
    max_weight::Int,
    base_graph::DiGraph,
    competitor_nodes,
    decision_arcs,
    competitor_cost_factor::Float64,
    user_parameter_mode::String,
)
    rng = MersenneTwister(seed)
    users = User[]
    max_possible_weight = 3 * length(base_nodes) * max_weight
    feasible_targets = _feasible_od_targets(base_graph)

    if user_parameter_mode == "shared"
        rcost, rrisk, rweight = _random_arc_matrices(
            graph,
            seed;
            max_cost=max_cost,
            max_risk=max_risk,
            max_weight=max_weight,
            negative_risk=true,
        )
        _apply_two_layer_competition_arc_rules!(
            graph,
            base_graph,
            competitor_nodes,
            decision_arcs,
            rcost,
            rrisk,
            rweight,
            competitor_cost_factor,
        )
        user_weight = constrained ? rweight : nothing
        minweights = constrained ? floyd_warshall_shortest_paths(graph, rweight) : nothing
        for user_id in 1:nusers
            base_origin, base_destination = _sample_feasible_od_pair(rng, feasible_targets)
            origin = competitor_nodes[base_origin]
            destination = competitor_nodes[base_destination]
            bound = _compute_weight_bound(minweights, origin, destination, length_slack, max_possible_weight)
            push!(users, User("U$(user_id)", origin, destination, rrisk, rcost, user_weight, bound))
        end
        return users, minweights
    end

    for user_id in 1:nusers
        user_seed = seed + 1000 * user_id
        rcost, rrisk, rweight = _random_arc_matrices(
            graph,
            user_seed;
            max_cost=max_cost,
            max_risk=max_risk,
            max_weight=max_weight,
            negative_risk=true,
        )
        _apply_two_layer_competition_arc_rules!(
            graph,
            base_graph,
            competitor_nodes,
            decision_arcs,
            rcost,
            rrisk,
            rweight,
            competitor_cost_factor,
        )
        user_weight = constrained ? rweight : nothing
        user_minweights = constrained ? floyd_warshall_shortest_paths(graph, rweight) : nothing
        base_origin, base_destination = _sample_feasible_od_pair(rng, feasible_targets)
        origin = competitor_nodes[base_origin]
        destination = competitor_nodes[base_destination]
        bound = _compute_weight_bound(user_minweights, origin, destination, length_slack, max_possible_weight)
        push!(users, User("U$(user_id)", origin, destination, rrisk, rcost, user_weight, bound))
    end
    return users, nothing
end

function _compute_weight_bound(minweights, origin, destination, length_slack, fallback_max_weight)
    if minweights === nothing || length_slack === nothing
        return nothing
    end
    length_slack = Float64(length_slack)
    0.0 <= length_slack <= 1.0 || throw(ArgumentError("length_slack must be in [0, 1]."))
    shortest_weight = minweights.dists[origin, destination]
    return ceil(Int, shortest_weight + length_slack * (fallback_max_weight - shortest_weight))
end

function _build_layered_competition_graph(base_graph::DiGraph)
    n = nv(base_graph)
    graph = DiGraph(3 * n)
    base_nodes = collect(1:n)
    layer1_nodes = collect((n + 1):(2 * n))
    layer2_nodes = collect((2 * n + 1):(3 * n))

    for node in base_nodes
        add_edge!(graph, node, node + n)
        add_edge!(graph, node + n, node)
        add_edge!(graph, node, node + 2n)
        add_edge!(graph, node + 2n, node)
    end

    for layer_offset in (n, 2n)
        for edge in edges(base_graph)
            add_edge!(graph, src(edge) + layer_offset, dst(edge) + layer_offset)
        end
    end

    decision_arcs = Tuple{Int,Int}[(src(e), dst(e)) for e in edges(graph) if src(e) in layer2_nodes && dst(e) in layer2_nodes]
    return graph, base_nodes, layer2_nodes, decision_arcs
end

function _build_two_layer_competition_graph(base_graph::DiGraph)
    n = nv(base_graph)
    graph = DiGraph(2 * n)
    layer1_nodes = collect(1:n)
    layer2_nodes = collect((n + 1):(2 * n))

    for node in 1:n
        add_edge!(graph, node, node + n)
        add_edge!(graph, node + n, node)
    end

    for edge in edges(base_graph)
        add_edge!(graph, src(edge), dst(edge))
        add_edge!(graph, src(edge) + n, dst(edge) + n)
    end

    decision_arcs = Tuple{Int,Int}[(src(e), dst(e)) for e in edges(graph) if src(e) in layer2_nodes && dst(e) in layer2_nodes]
    return graph, layer1_nodes, layer1_nodes, decision_arcs
end

function _apply_competition_arc_rules!(
    graph::DiGraph,
    base_graph::DiGraph,
    base_nodes,
    competitor_nodes,
    decision_arcs,
    rcost,
    rrisk,
    rweight,
    competitor_cost_factor::Float64,
)
    nbase = nv(base_graph)
    layer1_nodes = (nbase + 1):(2 * nbase)
    layer2_nodes = (2 * nbase + 1):(3 * nbase)
    decision_arc_set = Set(decision_arcs)

    for (i, j) in Iterators.product(1:nv(graph), 1:nv(graph))
        has_edge(graph, i, j) || continue

        if !(i in competitor_nodes || j in competitor_nodes)
            rrisk[i, j] = 0
        end

        if i in layer2_nodes && j in base_nodes
            rcost[i, j] = 0
            rrisk[i, j] = 0
            rweight[i, j] = 0
        elseif (i in layer1_nodes || i in layer2_nodes) && j in base_nodes
            rcost[i, j] = 0
            rweight[i, j] = 0
        end

        if (i, j) in decision_arc_set
            rcost[i, j] = round(Int, competitor_cost_factor * rcost[i, j])
        end
    end
end

function _apply_two_layer_competition_arc_rules!(
    graph::DiGraph,
    base_graph::DiGraph,
    competitor_nodes,
    decision_arcs,
    rcost,
    rrisk,
    rweight,
    competitor_cost_factor::Float64,
)
    nbase = nv(base_graph)
    decision_nodes = (nbase + 1):(2 * nbase)
    decision_arc_set = Set(decision_arcs)

    for (i, j) in Iterators.product(1:nv(graph), 1:nv(graph))
        has_edge(graph, i, j) || continue

        if (i, j) in decision_arc_set
            rcost[i, j] = round(Int, competitor_cost_factor * rcost[i, j])
            continue
        end

        rrisk[i, j] = 0

        if (i in competitor_nodes && j in decision_nodes) || (i in decision_nodes && j in competitor_nodes)
            rcost[i, j] = 0
            rweight[i, j] = 0
        end
    end
end

function _generated_name(
    base_name::String,
    topology_id::String,
    nusers::Int,
    parameter_seed::Int,
    constrained::Bool;
    length_slack=nothing,
    competitor_cost_factor=nothing,
    two_stage=nothing,
    user_parameter_mode=nothing,
    decision_arc_count=nothing,
)
    parts = [
        base_name,
        topology_id,
        "U$(nusers)",
        "S$(parameter_seed)",
        constrained ? "len" : "nolen",
    ]
    length_slack !== nothing && push!(parts, "LS$(length_slack)")
    competitor_cost_factor !== nothing && push!(parts, "CCF$(competitor_cost_factor)")
    decision_arc_count !== nothing && push!(parts, "K$(decision_arc_count)")
    two_stage !== nothing && push!(parts, two_stage ? "coop" : "bilevel")
    user_parameter_mode !== nothing && push!(parts, user_parameter_mode == "shared" ? "shared" : "peruser")
    return join(parts, "_")
end

function _write_hndp_users_json(path::String, generated_network::HNDPGeneratedNetwork)
    data = Dict(
        "instance_name" => generated_network.name,
        "instance_type" => generated_network.metadata["instance_type"],
        "users" => [
            Dict(
                "name" => string(user.uname),
                "origin" => user.origin,
                "destination" => user.destination,
                "weight_limit" => user.weighlimit,
            ) for user in generated_network.instance.users
        ],
    )
    open(path, "w") do io
        JSON.print(io, data, 2)
    end
end

function _write_hndp_gexf(path::String, generated_network::HNDPGeneratedNetwork)
    inst = generated_network.instance
    graph = inst.mygraph
    users = inst.users
    decision_arc_set = Set(inst.edgeA)
    node_attrs, edge_attrs = _gexf_attribute_specs(inst, users)

    open(path, "w") do io
        write(io, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
        write(io, "<gexf xmlns=\"http://www.gexf.net/1.2draft\" version=\"1.2\">\n")
        write(io, "  <meta lastmodifieddate=\"2026-04-13\">\n")
        write(io, "    <creator>JuBiC HNDP Network Generator</creator>\n")
        write(io, "    <description>" * _xml_escape(generated_network.name) * "</description>\n")
        write(io, "  </meta>\n")
        write(io, "  <graph mode=\"static\" defaultedgetype=\"directed\">\n")
        _write_gexf_attributes_block(io, "node", node_attrs)
        _write_gexf_attributes_block(io, "edge", edge_attrs)
        write(io, "    <nodes>\n")
        for node in vertices(graph)
            label = string(node)
            write(io, "      <node id=\"" * _xml_escape(string(node)) * "\" label=\"" * _xml_escape(label) * "\">\n")
            node_values = _node_attribute_values(node, inst, users)
            _write_gexf_attvalues(io, node_values, 8)
            write(io, "      </node>\n")
        end
        write(io, "    </nodes>\n")
        write(io, "    <edges>\n")
        for (edge_index, edge) in enumerate(edges(graph))
            src_node = src(edge)
            dst_node = dst(edge)
            write(
                io,
                "      <edge id=\"" * string(edge_index) * "\" source=\"" * string(src_node) * "\" target=\"" * string(dst_node) * "\">\n",
            )
            edge_values = _edge_attribute_values((src_node, dst_node), inst, users, decision_arc_set)
            _write_gexf_attvalues(io, edge_values, 8)
            write(io, "      </edge>\n")
        end
        write(io, "    </edges>\n")
        write(io, "  </graph>\n")
        write(io, "</gexf>\n")
    end
end

function _gexf_attribute_specs(inst::HNDPwC, users)
    node_attrs = [
        ("layer", "string"),
    ]
    edge_attrs = [
        ("decision_arc", "boolean"),
        ("construction_cost", "double"),
    ]
    for user in users
        uname = string(user.uname)
        push!(node_attrs, ("origin_" * uname, "boolean"))
        push!(node_attrs, ("destination_" * uname, "boolean"))
        push!(edge_attrs, ("cost_" * uname, "double"))
        push!(edge_attrs, ("risk_" * uname, "double"))
        push!(edge_attrs, ("weight_" * uname, "double"))
    end
    return node_attrs, edge_attrs
end

function _write_gexf_attributes_block(io, class_name::String, attrs)
    write(io, "    <attributes class=\"" * class_name * "\">\n")
    for (index, (title, attr_type)) in enumerate(attrs)
        write(
            io,
            "      <attribute id=\"" * string(index - 1) * "\" title=\"" * _xml_escape(title) * "\" type=\"" * attr_type * "\"/>\n",
        )
    end
    write(io, "    </attributes>\n")
end

function _write_gexf_attvalues(io, values, indent_spaces::Int)
    indent = " "^indent_spaces
    write(io, indent * "<attvalues>\n")
    for (attr_id, value) in values
        write(
            io,
            indent * "  <attvalue for=\"" * string(attr_id) * "\" value=\"" * _xml_escape(string(value)) * "\"/>\n",
        )
    end
    write(io, indent * "</attvalues>\n")
end

function _node_attribute_values(node, inst::HNDPwC, users)
    values = Vector{Tuple{Int,Any}}()
    push!(values, (0, _node_layer_label(node, inst)))
    attr_id = 1
    for user in users
        push!(values, (attr_id, node == user.origin))
        push!(values, (attr_id + 1, node == user.destination))
        attr_id += 2
    end
    return values
end

function _edge_attribute_values(edge::Tuple{Int,Int}, inst::HNDPwC, users, decision_arc_set)
    values = Vector{Tuple{Int,Any}}()
    push!(values, (0, edge in decision_arc_set))
    push!(values, (1, get(inst.edge_price, edge, 0)))
    attr_id = 2
    for user in users
        push!(values, (attr_id, _edge_matrix_value(user.mcost, edge)))
        push!(values, (attr_id + 1, _edge_matrix_value(user.mrisk, edge)))
        push!(values, (attr_id + 2, _edge_matrix_value(user.mweight, edge)))
        attr_id += 3
    end
    return values
end

function _edge_matrix_value(matrix, edge::Tuple{Int,Int})
    matrix === nothing && return 0
    return matrix[edge[1], edge[2]]
end

function _node_layer_label(node::Int, inst::HNDPwC)
    graph = inst.mygraph
    n = nv(graph)
    if n % 3 == 0 && !isempty(inst.edgeA) && maximum(max(edge[1], edge[2]) for edge in inst.edgeA) > n ÷ 3
        base_n = n ÷ 3
        if node <= base_n
            return "base"
        elseif node <= 2 * base_n
            return "layer_1"
        else
            return "layer_2"
        end
    end
    return "single_layer"
end

function _xml_escape(value::String)
    escaped = replace(value, "&" => "&amp;")
    escaped = replace(escaped, "<" => "&lt;")
    escaped = replace(escaped, ">" => "&gt;")
    escaped = replace(escaped, "\"" => "&quot;")
    escaped = replace(escaped, "'" => "&apos;")
    return escaped
end
