using Test
using JuMP
using Gurobi
import MathOptInterface as MOI

include("hndp_model_generation.jl")

function multimodal_smoke(; solve::Bool=false)
    config = Dict{String,Any}(
        "parameter_seeds" => [17],
        "instances" => [
            Dict{String,Any}(
                "name" => "multimodal_smoke",
                "instance_type" => "multimodal_bike",
                "topologies" => ["sioux_falls"],
                "nusers" => [8],
                "station_limit" => 6,
                "station_budget" => 3,
            ),
            Dict{String,Any}(
                "name" => "multimodal_smoke",
                "instance_type" => "multimodal_expansion",
                "topologies" => ["sioux_falls"],
                "nusers" => [8],
                "station_limit" => 6,
                "expansion_station_count" => 3,
                "arc_budget" => 4,
            ),
        ],
    )
    generated = generate_hndp_networks(config)
    @test length(generated) == 2
    @test generated[1].metadata["stations"] == generated[2].metadata["stations"]
    base_stations = generated[1].metadata["stations"]
    for i in base_stations, j in base_stations
        i == j && continue
        @test generated[1].instance.users[1].mcost[i + 24, j + 24] == generated[2].instance.users[1].mcost[i + 24, j + 24]
    end
    for network in generated
        hndp = network.instance
        @test all(user.mweight === nothing for user in hndp.users)
        @test eltype(hndp.users[1].mcost) == Int
        @test eltype(hndp.users[1].mrisk) == Int
        @test all(value isa Int for value in values(hndp.edge_price))
        @test all(hndp.users[1].mcost[a...] > 0 for a in _hndp_all_arcs(hndp) if hndp.users[1].mcost[a...] != 0)
        @test length(hndp.decision_groups) == length(hndp.edgeA) ÷ (network.metadata["scenario"] == "bike" ? 2 : 1)
        if network.metadata["scenario"] == "bike"
            @test all(length(group) == 2 for group in hndp.decision_groups)
            @test network.metadata["availability_budget_count"] == 2 * network.metadata["decision_budget"]
        else
            @test all(length(group) == 1 for group in hndp.decision_groups)
            @test network.metadata["availability_budget_count"] == network.metadata["decision_budget"]
            @test all(arc[1] > 24 && arc[2] > 24 for arc in hndp.edgeA)
            b = length(network.metadata["stations"])
            a = length(network.metadata["additional_stations"])
            @test length(hndp.edgeA) == (b + a) * (b + a - 1) - b * (b - 1)
        end
        if solve
            model = build_hndp_sd_instance(
                hndp,
                GurobiSolver();
                big_m_mode=HNDP_BIGM_FIXED_NETWORK_PATH,
                availability_budget_count=network.metadata["availability_budget_count"],
            )
            set_time_limit_sec(model.master.mymip, 60.0)
            set_silent(model.master.mymip)
            optimize!(model.master.mymip)
            println(network.metadata["scenario"], ": status=", termination_status(model.master.mymip),
                ", objective=", has_values(model.master.mymip) ? objective_value(model.master.mymip) : "n/a",
                ", nodes=", network.metadata["nnodes"], ", arcs=", network.metadata["narcs"],
                ", decisions=", network.metadata["decision_arcs"])
            @test termination_status(model.master.mymip) in
                (MOI.OPTIMAL, MOI.TIME_LIMIT, MOI.FEASIBLE_POINT, MOI.INTERRUPTED)
        end
    end
    return generated
end

multimodal_smoke(solve=("--solve" in ARGS))
println("Multimodal Sioux Falls smoke test passed.")
