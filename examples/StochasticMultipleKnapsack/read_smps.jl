### Provides functionalities to work with the SMPS format from SIPLIB (https://www2.isye.gatech.edu/~sahmed/siplib/)
###  Code mainly generated with ChatGPT

using JuMP

# -------------------------
# Minimal MPS parser by reading msp file in JuMP
# -------------------------

struct Knapcksack_Master 
    my_model::JuMP.Model
    p_vars::Dict{String, JuMP.VariableRef} # name-> variable
    x_vars::Dict{String, JuMP.VariableRef} # name-> variable

    info_second_stage::Dict{String,NamedTuple{(:sense,:rhs,:q,:x), Tuple{String,Float64,Dict{String,Float64},Dict{String,Float64}}}}
end

function read_knapsack_mps_to_model(cor_path::AbstractString)
    ## read mps file
    mps_model = MOI.FileFormats.Model(format = MOI.FileFormats.FORMAT_MPS)
    MOI.read_from_file(mps_model, cor_path)

    ## transform to JuMP model
    jump = Model(()->Gurobi.Optimizer())
    MOI.copy_to(jump, mps_model)

    ## collect variable information 
    # Get all variables
    vars = all_variables(jump)
    # Group them by prefix
    p_vars = Dict{String, JuMP.VariableRef}(name(v) => v for v in vars if startswith(name(v), "p"))
    x_vars = Dict{String, JuMP.VariableRef}(name(v) => v for v in vars if startswith(name(v), "x"))
    q_vars = [v for v in vars if startswith(name(v), "q")]

    ## now collect information about the second stage constraints
    targets = ["c51", "c52", "c53", "c54", "c55"]

    # Collect info here
    saved = Dict{String,NamedTuple{(:sense,:rhs,:q,:x),
                                Tuple{String,Float64,Dict{String,Float64},Dict{String,Float64}}}}()

    for cname in targets
        cref = constraint_by_name(jump, cname)
        @assert cref !== nothing "No constraint named $cname in model"

        obj = JuMP.constraint_object(cref)
        f, s = obj.func, obj.set
        #@assert f isa MOI.ScalarAffineFunction{Float64} "Only linear scalar constraints are handled but $(typeof(f)) is not"

        # Gather coefficients of q* and x*
        q = Dict{String,Float64}()
        x = Dict{String,Float64}()
        for t in f.terms
            #println("t=$t, t0 = $(typeof(t[1])), t2=$(t[2])")
            vname = name(t[1])
            if startswith(vname, "q")
                q[vname] = get(q, vname, 0.0) + t[2]
            elseif startswith(vname, "x")
                x[vname] = get(x, vname, 0.0) + t[2]
            end
        end

        # Normalize sense and RHS: a'x + c ≤ ub  ->  a'x ≤ (ub - c)
        #println("s=$s")
        sense, rhs = if s isa MOI.LessThan{Float64}
            ("<=", s.upper - f.constant)
        elseif s isa MOI.GreaterThan{Float64}
            (">=", s.lower - f.constant)
        elseif s isa MOI.EqualTo{Float64}
            ("==", s.value - f.constant)
        else
            error("Unhandled set type $(typeof(s)) for $cname")
        end

        saved[cname] = (sense = sense, rhs = rhs, q = q, x = x)
    end

    # Remove the constraint from the model
    for cname in targets
        cref = constraint_by_name(jump, cname)
        delete(jump, cref)
    end

    # Remove the q variables 
    for qvar in q_vars
        delete(jump, qvar)
    end

    # return 
    return Knapcksack_Master(jump, p_vars, x_vars, saved)
end

# -------------------------
# Parse STO and a vector of scenarios
# -------------------------

function read_q_scenarios(path::AbstractString)::Vector{Dict{String,Float64}}
    scen_qs = Vector{Dict{String,Float64}}()
    current = nothing

    # regex for scenario header lines starting with "SC"
    sc_header = r"^\s*SC\b"
    # regex for lines like: "    q_00      obj                 17"
    qline = r"^\s*(\S+)\s+obj\s+([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s*$"

    open(path, "r") do io
        for line in eachline(io)
            if occursin(sc_header, line)
                # start a new scenario block
                if current !== nothing
                    push!(scen_qs, current)
                end
                current = Dict{String,Float64}()
            else
                m = match(qline, line)
                if m !== nothing
                    vname = m.captures[1]
                    if startswith(vname, "q")
                        current[vname] = parse(Float64, m.captures[2])
                    end
                end
            end
        end
    end

    # push the last scenario (if any)
    if current !== nothing
        push!(scen_qs, current)
    end

    return scen_qs
end

