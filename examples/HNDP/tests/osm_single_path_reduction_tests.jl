using Test
if !isdefined(@__MODULE__, :_build_hndp_model_from_spec)
    include("../experiments/hndp_osm_experiment_runner.jl")
end

function single_path_fixture(;budget=1,price=0.0,only_fixed=false)
    graph=DiGraph(4);cost=spzeros(Float64,4,4);risk=spzeros(Float64,4,4)
    for (u,v,c) in [(1,2,1.0),(2,4,10.0),(2,3,1.0),(3,2,1.0),(3,4,1.0)]
        add_edge!(graph,u,v);cost[u,v]=c
    end
    risk[1,2]=7.0;risk[2,3]=-5.0
    users=User[User("fixed",1,2,risk,cost,nothing,nothing)]
    only_fixed || push!(users,User("choice",2,4,risk,cost,nothing,nothing))
    groups=[[(2,3),(3,2)]]
    arcs,prices=hndp_sharing_station_decisions(groups;construction_cost=price)
    HNDPwC(graph,users,arcs,prices,nothing,groups,budget)
end

function single_path_reduction_tests(output_root)
    mkpath(output_root)
    @testset "Fixed singleton elimination and full reconstruction" begin
        for parallel in (false,true),dominance in (false,true),(budget,price) in [(0,0.0),(1,0.0),(1,3.0)]
            h=single_path_fixture(;budget,price)
            built=[];counts=[]
            for eliminate in (false,true)
                instance,_,c=build_hndp_path_instance(h,GurobiSolver();enumeration_time_limit=60.0,parallelize=parallel,use_decision_arc_dominance=dominance,eliminate_single_path_users=eliminate)
                m=_hndp_master_model(instance)
                @test_throws ArgumentError hndp_path_solution_weights(m)
                @test_throws ArgumentError mip_solution_values(m)
                set_silent(m);set_optimizer_attribute(m,"Threads",1);optimize!(m)
                @test termination_status(m)==MOI.OPTIMAL
                @test objective_value(m)≈(budget==0 ? 7.0 : 2.0+price) atol=1e-8
                weights=hndp_path_solution_weights(m)
                @test length(weights)==2 && weights["fixed"]==[1.0]
                @test all(isapprox(sum(v),1.0;atol=1e-8) for v in values(weights))
                exported_values=mip_solution_values(m)
                @test exported_values["λ_fixed[1]"]==1 && exported_values["pi_path_fixed"]==1 && exported_values["optL2_path_fixed"]==1
                @test hndp_path_reduction_metadata(m)["eliminated_fixed_path_users"]==Int(eliminate)*(budget==0 ? 2 : 1)
                @test isnothing(variable_by_name(m,"λ_fixed[1]"))==eliminate
                if eliminate
                    restored=hndp_reconstructed_fixed_user_solution(m)
                    @test restored["fixed"]["arcs"]==[(1,2)]
                    @test restored["fixed"]["leader_contribution"]==7
                end
                push!(built,m);push!(counts,c)
            end
            @test counts[1]==counts[2]
            @test num_variables(built[1])-num_variables(built[2])==3*(budget==0 ? 2 : 1)
            @test num_constraints(built[1];count_variable_in_set_constraints=false)-num_constraints(built[2];count_variable_in_set_constraints=false)==5*(budget==0 ? 2 : 1)
            reconstructed=mip_solution_values(built[2])
            point=Dict(v=>reconstructed[name(v)] for v in all_variables(built[1]))
            @test isempty(primal_feasibility_report(built[1],point;atol=1e-6))
            @test value(v->point[v],objective_function(built[1]))≈objective_value(built[2]) atol=1e-8
        end
    end
    @testset "All followers fixed, safe guards, and generic solution export" begin
        h=single_path_fixture(;only_fixed=true,price=3)
        instance,_,_=build_hndp_path_instance(h,GurobiSolver();enumeration_time_limit=60.0)
        m=_hndp_master_model(instance)
        @test hndp_path_reduction_metadata(m)["modeled_path_users"]==0
        folder=joinpath(output_root,"export_"*string(time_ns()));mkpath(folder)
        param=MIPparam(GurobiSolver(),true,folder,"lp",JuBiC.RunStats(),60.0,20261102,1)
        JuBiC.solve_with_MIP!(instance,param)
        @test objective_value(m)≈7.0 atol=1e-8
        exported=JSON.parsefile(joinpath(folder,"solution.json"))
        @test exported["λ_fixed[1]"]==1 && exported["optL2_path_fixed"]==1
        restored=JSON.parsefile(joinpath(folder,"reconstructed_solution.json"))
        @test restored["fixed"]["arcs"]==[[1,2]]
        @test restored["fixed"]["demand"]==1
        # A partial path set is not a proof of fixed behavior.
        dummy=Model();path=HNDPPath([(1,2)],1.0,7.0,0.0)
        @test !_try_eliminate_fixed_path_user!(dummy,only(h.users),[path],h.edgeA;complete=false)
        # A singleton requiring a decision arc must keep its availability constraints.
        for parallel in (false,true),budget in (0,1)
            graph=DiGraph(2);add_edge!(graph,1,2);c=spzeros(2,2);c[1,2]=1;r=spzeros(2,2);r[1,2]=-5
            user=User("dependent",1,2,r,c,nothing,nothing)
            dependent=HNDPwC(graph,[user],[(1,2)],Dict((1,2)=>0.0),nothing,[[(1,2)]],budget)
            inst,_,counts=build_hndp_path_instance(dependent,GurobiSolver();enumeration_time_limit=60.0,parallelize=parallel)
            mm=_hndp_master_model(inst);set_silent(mm);set_optimizer_attribute(mm,"Threads",1);set_optimizer_attribute(mm,"DualReductions",0);optimize!(mm)
            @test counts["dependent"]==(budget==0 ? 0 : 1)
            @test hndp_path_reduction_metadata(mm)["eliminated_fixed_path_users"]==0
            @test termination_status(mm)==(budget==0 ? MOI.INFEASIBLE : MOI.OPTIMAL)
        end
    end
end





