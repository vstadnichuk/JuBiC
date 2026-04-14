using JSON
using JuMP
import MathOptInterface as MOI

const GBC_INSTANCE_IO_FORMAT = "JuBiC_GBC_JuMP_v1"

function output_GBC_solver_instance(
    instance::Instance,
    output_directory::AbstractString;
    master_filename::AbstractString="master.lp",
)
    instance.master isa Master || throw(ArgumentError("GBC solver instance export requires a JuBiC `Master`."))
    master = instance.master
    subs = instance.subproblems

    isnothing(master.partial_decomposition) || throw(ArgumentError("GBC solver instance export currently supports only masters without `partial_decomposition`."))
    isnothing(master.objL2) || throw(ArgumentError("GBC solver instance export currently supports only masters without `objL2` information."))
    all(sub -> sub isa SubSolverJuMP, subs) || throw(ArgumentError("GBC solver instance export currently supports only `SubSolverJuMP` subproblems."))

    mkpath(output_directory)
    master_path = joinpath(output_directory, String(master_filename))
    write_to_file(master.model, master_path)

    metadata = Dict{String,Any}(
        "format" => GBC_INSTANCE_IO_FORMAT,
        "master_file" => basename(master_path),
        "A" => [_serialize_gbc_key(a) for a in master.A],
        "sub_names" => collect(master.sub_names),
        "master_link_vars" => Dict(_serialize_gbc_key(a) => JuMP.name(master.link_vars[a]) for a in master.A),
        "subproblems" => Any[],
    )

    for sub_any in subs
        sub = sub_any::SubSolverJuMP
        sub_filename = "sub_$(sub.name).lp"
        write_to_file(sub.mip_model, joinpath(output_directory, sub_filename))
        push!(metadata["subproblems"], _serialize_gbc_subproblem(sub, sub_filename))
    end

    open(joinpath(output_directory, "metadata.json"), "w") do io
        write(io, JSON.json(metadata, 2))
    end
    return nothing
end

function read_GBC_solver_instance(
    input_directory::AbstractString,
    solver::SolverWrapper,
)
    metadata_path = joinpath(input_directory, "metadata.json")
    isfile(metadata_path) || throw(ArgumentError("No GBC metadata file found at $(metadata_path)."))
    metadata = JSON.parsefile(metadata_path)

    String(get(metadata, "format", "")) == GBC_INSTANCE_IO_FORMAT ||
        throw(ArgumentError("Unsupported GBC instance format '$(get(metadata, "format", nothing))'."))

    master_model = _read_jump_model_from_lp(
        joinpath(input_directory, String(metadata["master_file"])),
        solver,
    )
    set_silent(master_model)

    A = [String(a) for a in metadata["A"]]
    sub_names = [String(name) for name in metadata["sub_names"]]
    master_link_vars = Dict(
        a => _required_variable_by_name(master_model, String(metadata["master_link_vars"][a]))
        for a in A
    )
    master = Master(master_model, A, master_link_vars, sub_names)

    subs = SubSolver[]
    for raw_sub in metadata["subproblems"]
        sub_meta = Dict{String,Any}(raw_sub)
        push!(subs, _read_gbc_subproblem(input_directory, solver, sub_meta))
    end

    return Instance(master, subs)
end

function _serialize_gbc_subproblem(sub::SubSolverJuMP, sub_filename::AbstractString)
    return Dict{String,Any}(
        "name" => sub.name,
        "file" => String(sub_filename),
        "A" => [_serialize_gbc_key(a) for a in sub.A],
        "link_varsC" => Dict(_serialize_gbc_key(a) => JuMP.name(sub.link_varsC[a]) for a in sub.A),
        "y_vars" => Dict(_serialize_gbc_key(a) => JuMP.name(sub.y_vars[a]) for a in sub.A),
        "r_objterm" => _serialize_linear_expression(sub.r_objterm),
        "c_objterm" => _serialize_linear_expression(sub.c_objterm),
    )
end

function _read_gbc_subproblem(
    input_directory::AbstractString,
    solver::SolverWrapper,
    sub_meta::Dict{String,Any},
)
    sub_model = _read_jump_model_from_lp(joinpath(input_directory, String(sub_meta["file"])), solver)
    set_silent(sub_model)

    A = [String(a) for a in sub_meta["A"]]
    link_varsC = Dict(
        a => _required_variable_by_name(sub_model, String(sub_meta["link_varsC"][a]))
        for a in A
    )
    y_vars = Dict(
        a => _required_variable_by_name(sub_model, String(sub_meta["y_vars"][a]))
        for a in A
    )
    r_objterm = _deserialize_linear_expression(sub_model, Dict{String,Any}(sub_meta["r_objterm"]))
    c_objterm = _deserialize_linear_expression(sub_model, Dict{String,Any}(sub_meta["c_objterm"]))

    return SubSolverJuMP{String}(
        String(sub_meta["name"]),
        sub_model,
        A,
        link_varsC,
        y_vars,
        r_objterm,
        c_objterm,
        _default_extra_cuts(),
    )
end

function _serialize_linear_expression(expr)
    if expr isa Number
        return Dict{String,Any}("constant" => expr, "terms" => Any[])
    elseif expr isa JuMP.VariableRef
        return Dict{String,Any}(
            "constant" => 0.0,
            "terms" => Any[Dict("var" => JuMP.name(expr), "coef" => 1.0)],
        )
    elseif expr isa JuMP.GenericAffExpr
        return Dict{String,Any}(
            "constant" => expr.constant,
            "terms" => Any[
                Dict("var" => JuMP.name(var), "coef" => coeff) for (var, coeff) in expr.terms
            ],
        )
    end
    throw(ArgumentError("Only linear JuMP expressions can be exported in the current GBC solver instance format. Got $(typeof(expr))."))
end

function _deserialize_linear_expression(model::JuMP.Model, meta::Dict{String,Any})
    expr = AffExpr(Float64(get(meta, "constant", 0.0)))
    for term_any in get(meta, "terms", Any[])
        term = Dict{String,Any}(term_any)
        add_to_expression!(expr, Float64(term["coef"]), _required_variable_by_name(model, String(term["var"])))
    end
    return expr
end

function _serialize_gbc_key(a)
    return JSON.json(a)
end

function _required_variable_by_name(model::JuMP.Model, var_name::String)
    for candidate in _gbc_variable_name_candidates(var_name)
        var = variable_by_name(model, candidate)
        if !isnothing(var)
            return var
        end
    end
    throw(ArgumentError("Variable '$(var_name)' was not found when reconstructing a GBC solver instance."))
end

function _read_jump_model_from_lp(path::AbstractString, solver::SolverWrapper)
    isfile(path) || throw(ArgumentError("Expected LP file at $(path), but it does not exist."))
    lp_model = MOI.FileFormats.Model(format=MOI.FileFormats.FORMAT_LP)
    MOI.read_from_file(lp_model, String(path))
    jump_model = Model(() -> get_next_optimizer(solver))
    MOI.copy_to(jump_model, lp_model)
    return jump_model
end

function _gbc_variable_name_candidates(var_name::String)
    sanitized = replace(var_name, r"[^A-Za-z0-9_]" => "_")
    if sanitized == var_name
        return (var_name,)
    end
    return (var_name, sanitized)
end
