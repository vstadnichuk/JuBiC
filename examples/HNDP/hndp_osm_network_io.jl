# Data containers and lossless interchange for OSM-derived HNDP networks.
using JSON, Graphs, SparseArrays
if !isdefined(@__MODULE__, :HNDPwC)
    include("hndp_instances.jl")
end
if !isdefined(@__MODULE__, :HNDPGeneratedNetwork)
    struct HNDPGeneratedNetwork
        name::String
        instance::HNDPwC
        metadata::Dict{String,Any}
    end
end

function write_hndp_generated_networks(
    generated::Vector{HNDPGeneratedNetwork},
    output_cfg::Dict{String,Any},
)
    base_folder = String(get(output_cfg, "folder", "examples/HNDP/generated_instances"))
    graph_format = lowercase(String(get(output_cfg, "graph_format", "json")))
    graph_format in ("json", "gexf", "both") ||
        throw(ArgumentError("graph_format must be 'json', 'gexf', or 'both'."))

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

    if graph_format in ("gexf", "both")
        _write_hndp_gexf(joinpath(instance_folder, "instance.gexf"), generated_network)
    end
    # `instance.json` is the canonical, lossless interchange format.  Keep the
    # older users.json sidecar for compatibility with existing inspection tools.
    _write_hndp_instance_json(joinpath(instance_folder, "instance.json"), generated_network)
    _write_hndp_users_json(joinpath(instance_folder, "users.json"), generated_network)
end

function read_hndp_instance(path::AbstractString)
    data = JSON.parsefile(String(path))
    Int(get(data, "format_version", 0)) in (1, 2) ||
        throw(ArgumentError("Unsupported HNDP instance format in $(path)."))

    nodes = data["nodes"]
    n = length(nodes)
    graph = DiGraph(n)
    for edge_data in data["edges"]
        add_edge!(graph, Int(edge_data["source"]), Int(edge_data["target"]))
    end

    edge_price = Dict{Tuple{Int,Int},Float64}()
    edgeA = Tuple{Int,Int}[]
    for edge_data in data["edges"]
        edge = (Int(edge_data["source"]), Int(edge_data["target"]))
        edge_price[edge] = Float64(get(edge_data, "construction_cost", 0.0))
        Bool(get(edge_data, "decision_arc", false)) && push!(edgeA, edge)
    end

    users = User[]
    for user_data in data["users"]
        cost = spzeros(Float64, n, n)
        risk = spzeros(Float64, n, n)
        weight = spzeros(Float64, n, n)
        for edge_data in user_data["arcs"]
            i, j = Int(edge_data["source"]), Int(edge_data["target"])
            cost[i, j] = Float64(edge_data["cost"])
            risk[i, j] = Float64(edge_data["risk"])
            weight[i, j] = Float64(edge_data["weight"])
        end
        limit = get(user_data, "weight_limit", nothing)
        limit = isnothing(limit) ? nothing : Float64(limit)
        push!(users, User(
            user_data["name"], Int(user_data["origin"]), Int(user_data["destination"]),
            risk, cost, get(user_data, "has_weight", true) ? weight : nothing, limit,
        ))
    end

    # The shortest-weight matrix is a cache, not part of the mathematical
    # instance. Recompute it when a resource limit is present.
    minweights = nothing
    if any(!isnothing(user.weighlimit) for user in users)
        minweights = floyd_warshall_shortest_paths(graph, users[1].mweight)
    end
    raw_groups = get(data, "decision_groups", nothing)
    groups = isnothing(raw_groups) ? nothing :
        [[(Int(arc[1]), Int(arc[2])) for arc in group] for group in raw_groups]
    budget = get(data, "decision_budget", nothing)
    instance = HNDPwC(graph, users, edgeA, edge_price, minweights, groups, budget)
    metadata = Dict{String,Any}(get(data, "metadata", Dict{String,Any}()))
    name = String(get(data, "instance_name", get(metadata, "name", basename(dirname(String(path))))))
    metadata["name"] = name
    metadata["instance_file"] = abspath(String(path))
    return HNDPGeneratedNetwork(name, instance, metadata)
end

function visit_hndp_saved_instances(paths::AbstractVector, visitor::Function)
    seen = Set{String}()
    for raw_path in paths
        path = abspath(String(raw_path))
        path in seen && continue
        push!(seen, path)
        visitor(read_hndp_instance(path))
    end
    return nothing
end

function visit_hndp_saved_instances(config::Dict{String,Any}, visitor::Function)
    paths = get(config, "saved_instances", nothing)
    paths === nothing && throw(ArgumentError("A saved-instance configuration requires 'saved_instances'."))
    if paths isa AbstractString
        paths = [paths]
    end
    return visit_hndp_saved_instances(paths, visitor)
end

function _write_hndp_instance_json(path::String, generated_network::HNDPGeneratedNetwork)
    inst = generated_network.instance
    graph = inst.mygraph
    edge_records = Any[]
    for edge in edges(graph)
        arc = (src(edge), dst(edge))
        push!(edge_records, Dict(
            "source" => arc[1],
            "target" => arc[2],
            "decision_arc" => arc in Set(inst.edgeA),
            "construction_cost" => Float64(get(inst.edge_price, arc, 0.0)),
        ))
    end

    user_records = Any[]
    for user in inst.users
        arcs = Any[]
        for edge in edges(graph)
            arc = (src(edge), dst(edge))
            push!(arcs, Dict(
                "source" => arc[1],
                "target" => arc[2],
                "cost" => Float64(_edge_matrix_value(user.mcost, arc)),
                "risk" => Float64(_edge_matrix_value(user.mrisk, arc)),
                "weight" => Float64(_edge_matrix_value(user.mweight, arc)),
            ))
        end
        push!(user_records, Dict(
            "name" => string(user.uname),
            "origin" => Int(user.origin),
            "destination" => Int(user.destination),
            "weight_limit" => user.weighlimit,
            "has_weight" => !isnothing(user.mweight),
            "arcs" => arcs,
        ))
    end

    data = Dict(
        "format" => "JuBiC-HNDP-instance",
        # Older readers must reject grouped-decision instances rather than
        # silently dropping their additional constraints.
        "format_version" => 2,
        "instance_name" => generated_network.name,
        "metadata" => generated_network.metadata,
        "nodes" => collect(1:nv(graph)),
        "edges" => edge_records,
        "users" => user_records,
        "decision_groups" => isnothing(inst.decision_groups) ? nothing :
            [[[arc[1], arc[2]] for arc in group] for group in inst.decision_groups],
        "decision_budget" => inst.decision_budget,
    )
    open(path, "w") do io
        JSON.print(io, data, 2)
    end
end

function _write_hndp_users_json(path::String, generated_network::HNDPGeneratedNetwork)
    data = Dict(
        "instance_name" => generated_network.name,
        "instance_type" => get(generated_network.metadata, "instance_type", "HNDPwC"),
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

function _xml_escape(value::String)
    escaped = replace(value, "&" => "&amp;")
    escaped = replace(escaped, "<" => "&lt;")
    escaped = replace(escaped, ">" => "&gt;")
    escaped = replace(escaped, "\"" => "&quot;")
    escaped = replace(escaped, "'" => "&apos;")
    return escaped
end

function _node_layer_label(node::Int, inst::HNDPwC)
    n = nv(inst.mygraph)
    if !isnothing(inst.decision_groups) && n % 4 == 0
        return ("physical/walk", "car", "transit", "bike")[cld(node, n ÷ 4)]
    end
    return "physical"
end
