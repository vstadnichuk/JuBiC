using Test
if !isdefined(@__MODULE__, :_build_hndp_model_from_spec)
    include("../experiments/hndp_osm_experiment_runner.jl")
end
include("../experiments/osm_full_od_mip_benchmark.jl")
include("osm_smoke.jl")

function station_cost_fixture(budget)
    graph=DiGraph(8);cost=spzeros(Float64,8,8);risk=spzeros(Float64,8,8)
    groups=[[(1,7),(7,1)],[(2,8),(8,2)]]
    for (u,v,c) in [(1,2,10.0),(2,1,10.0),(7,8,1.0),(8,7,1.0),
                    (1,7,1.0),(7,1,0.1),(2,8,1.0),(8,2,0.1)]
        add_edge!(graph,u,v);cost[u,v]=c
    end
    risk[1,7]=-40.0;risk[2,8]=-40.0
    decisions,prices=hndp_sharing_station_decisions(groups;construction_cost=[3.0,5.0])
    user=User("u",1,2,risk,cost,nothing,nothing)
    HNDPwC(graph,[user],decisions,prices,nothing,groups,budget)
end

function station_decision_regression_tests(output_root)
    mkpath(output_root)
    @testset "Station construction helper" begin
        h=station_cost_fixture(2)
        @test length(h.edgeA)==4
        @test h.edge_price[(1,7)]==3 && h.edge_price[(2,8)]==5
        @test h.edge_price[(7,1)]==0 && h.edge_price[(8,2)]==0
        _,prices=hndp_sharing_station_decisions(h.decision_groups;construction_cost=9)
        @test sum(values(prices))==18
        @test_throws ArgumentError hndp_sharing_station_decisions(h.decision_groups;construction_cost=-1)
        @test_throws ArgumentError hndp_sharing_station_decisions(h.decision_groups;construction_cost=Inf)
        @test_throws ArgumentError hndp_sharing_station_decisions(h.decision_groups;construction_cost=[1])
        @test_throws ArgumentError hndp_sharing_station_decisions([[(1,7)]];construction_cost=1)
    end
    @testset "Coupling and entry-only cost through all four formulations" begin
        for kind in ("sd","blc","path_dominance","path_no_dominance"),budget in (0,1,2)
            original=station_cost_fixture(budget)
            directory=joinpath(output_root,"$(kind)_b$(budget)");mkpath(directory)
            path=joinpath(directory,"instance.json")
            _write_hndp_instance_json(path,HNDPGeneratedNetwork("station_fixture",original,Dict{String,Any}()))
            h=read_hndp_instance(path).instance
            @test Set(h.edgeA)==Set(original.edgeA)
            @test all(h.edge_price[a]==original.edge_price[a] for a in h.edgeA)
            @test h.decision_groups==original.decision_groups && h.decision_budget==budget
            spec=Dict{String,Any}("name"=>kind,"model_type"=>(kind in ("sd","blc") ? kind : "path"),
                "big_m_mode"=>"fixed_network_path_current_cost","subproblem_method"=>"mip",
                "parallelize"=>false,"use_decision_arc_dominance"=>kind=="path_dominance",
                "enumeration_time_limit"=>60.0,"fixed_path_bound_method"=>"labeling")
            cfg=Dict{String,Any}("name"=>"test","mip_solver"=>"Gurobi","runtime"=>60.0,
                "threads_master"=>1,"threads_sub_con"=>1,"parallel_separation"=>false,"integer_obj"=>false)
            instance,meta=_build_hndp_model_from_spec(h,Dict{String,Any}(),GurobiSolver(silent=true),spec,cfg)
            master=_hndp_master_model(instance);x=master[:x]
            construction=constraint_by_name(master,"construction_cost")
            @test normalized_coefficient(construction,x[(1,7)])==-3
            @test normalized_coefficient(construction,x[(2,8)])==-5
            @test normalized_coefficient(construction,x[(7,1)])==0
            @test normalized_coefficient(construction,x[(8,2)])==0
            resolved=_resolve_hndp_param_config(cfg,JuBiC._infer_solver_name(instance),directory,meta,false)
            stats=solve_instance!(instance,JuBiC._build_solver_params(JuBiC._infer_solver_name(instance),resolved))
            @test stats.data["Opt_status"]=="Optimal"
            @test stats.data["Opt"] ≈ (budget==2 ? -32.0 : 0.0) atol=1e-6
            @test value(x[(1,7)]) ≈ value(x[(7,1)]) atol=1e-6
            @test value(x[(2,8)]) ≈ value(x[(8,2)]) atol=1e-6
        end
        # A legacy entry-only decision set must not silently force stations open.
        good=station_cost_fixture(1)
        bad=HNDPwC(good.mygraph,good.users,first.(good.decision_groups),good.edge_price,nothing,good.decision_groups,1)
        @test_throws ArgumentError build_hndp_sd_instance(bad,GurobiSolver(silent=true))
        @test_throws ArgumentError build_hndp_blc_instance(bad,GurobiSolver(silent=true))
    end
    osm_smoke()
end
