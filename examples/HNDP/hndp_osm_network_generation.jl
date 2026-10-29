# OSM/GTFS-only network generation and saved-instance workflows for HNDP.
if !isdefined(@__MODULE__, :write_hndp_generated_networks)
    include("hndp_osm_network_io.jl")
end
if !isdefined(@__MODULE__, :generate_osm_gtfs)
    include("hndp_osm_gtfs_generation.jl")
end

load_hndp_osm_config(path::AbstractString) = JSON.parsefile(String(path))

function _osm_integer_values(value, name::String)
    values = value isa AbstractVector ? value : [value]
    isempty(values) && throw(ArgumentError("$name must not be empty."))
    all(x -> x isa Integer && !(x isa Bool), values) ||
        throw(ArgumentError("$name must contain integers."))
    return Int.(values)
end

function _visit_osm_spec!(raw_spec, default_seeds, output_cfg, visitor::Function)
    spec = Dict{String,Any}(raw_spec)
    source = get(spec, "instance_type", "osm_gtfs")
    source == "osm_gtfs" || throw(ArgumentError("Only instance_type=\"osm_gtfs\" is supported on this branch; use saved_instances for archived JSON networks."))
    for field in ("osm_path", "boundary_path", "gtfs_dir")
        get(spec, field, nothing) isa AbstractString ||
            throw(ArgumentError("An OSM instance requires the '$field' path."))
    end
    counts = _osm_integer_values(get(spec, "nusers", 100), "nusers")
    all(>(0), counts) || throw(ArgumentError("nusers must be positive."))
    seeds = _osm_integer_values(get(spec, "parameter_seeds", get(spec, "seed", default_seeds)), "parameter_seeds")
    base_name = String(get(spec, "name", get(spec, "city_name", "osm_gtfs")))
    ignored = Set(["name", "instance_type", "parameter_seeds", "seed", "nusers"])
    kwargs = Dict{Symbol,Any}(Symbol(k)=>v for (k,v) in spec if !(k in ignored))
    kwargs[:use_boundary_compression] = get(kwargs, :use_boundary_compression, true)
    for count in counts, seed in seeds
        generated = generate_osm_gtfs(; kwargs..., nusers=count, seed=seed)
        name = "$(base_name)_U$(count)_S$(seed)"
        generated.metadata["name"] = name
        generated.metadata["instance_type"] = "osm_gtfs"
        generated = HNDPGeneratedNetwork(name, generated.instance, generated.metadata)
        if !isnothing(output_cfg)
            write_hndp_generated_networks([generated], output_cfg)
        end
        visitor(generated)
    end
    return nothing
end

"""Visit cached OSM/GTFS instances one at a time, or read saved JSON instances.

Source paths are relative to the working directory; run from the repository
root as in the supplied examples. Other importer keyword arguments pass through
unchanged. Only user counts and seeds are expanded into a parameter grid.
"""
function visit_hndp_osm_networks(config::Dict{String,Any}, visitor::Function)
    if haskey(config, "saved_instances")
        haskey(config, "instances") && throw(ArgumentError("Use either instances or saved_instances, not both."))
        return visit_hndp_saved_instances(config, visitor)
    end
    specs = get(config, "instances", nothing)
    specs isa AbstractVector || throw(ArgumentError("An OSM config requires an 'instances' array or 'saved_instances'."))
    seeds = _osm_integer_values(get(config, "parameter_seeds", [0]), "parameter_seeds")
    output_cfg = haskey(config, "instance_output") ? Dict{String,Any}(config["instance_output"]) : nothing
    for spec in specs
        _visit_osm_spec!(spec, seeds, output_cfg, visitor)
    end
    return nothing
end

visit_hndp_osm_networks(path::AbstractString, visitor::Function) =
    visit_hndp_osm_networks(load_hndp_osm_config(path), visitor)

function generate_hndp_osm_networks(config::Dict{String,Any})
    generated = HNDPGeneratedNetwork[]
    visit_hndp_osm_networks(config, network -> push!(generated, network))
    return generated
end

generate_hndp_osm_networks(path::AbstractString) =
    generate_hndp_osm_networks(load_hndp_osm_config(path))
