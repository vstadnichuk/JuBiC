using Test
if !isdefined(@__MODULE__, :hndp_path_cleanup_tests)
    include("osm_path_cleanup_tests.jl")
end

function hndp_budget_pruning_tests()
    @testset "Budget-group mask construction" begin
        a=(1,2);b=(2,3);c=(3,4)
        make(groups,budget,arcs=[a,b,c])=HNDPwC(DiGraph(4),User[],arcs,Dict(x=>0.0 for x in arcs),nothing,groups,budget)
        @test isnothing(_hndp_enumeration_budget_data(make(nothing,1),Set([a,b,c])))
        @test isnothing(_hndp_enumeration_budget_data(make([[a,b]],nothing),Set([a,b,c])))
        @test isnothing(_hndp_enumeration_budget_data(make([[a,b]],1),Set([a,b,c])))
        d=_hndp_enumeration_budget_data(make([[a,b],[c]],1),Set([a,b,c]))
        @test d.empty_mask isa UInt64
        @test d.arc_masks[a]==d.arc_masks[b]
        @test count_ones(d.arc_masks[a]|d.arc_masks[b])==1
        @test count_ones(d.arc_masks[a]|d.arc_masks[c])==2
        overlap=_hndp_enumeration_budget_data(make([[a,b],[b,c]],1),Set([a,b,c]))
        @test overlap.arc_masks[a]==overlap.arc_masks[b]==overlap.arc_masks[c]
        @test count_ones(overlap.arc_masks[a])==2
        @test_throws ArgumentError _hndp_enumeration_budget_data(make([[a]],-1),Set([a,b,c]))
        @test_throws ArgumentError _hndp_enumeration_budget_data(make([[(8,9)]],0),Set([a,b,c]))
        many=[(i,i+100) for i in 1:70]
        large=make([[arc] for arc in many],64,many)
        big=_hndp_enumeration_budget_data(large,Set(many))
        @test big.empty_mask isa BigInt
        @test count_ones(foldl(|,values(big.arc_masks);init=big.empty_mask))==70
        @test big.limit==64
        boundary=_hndp_enumeration_budget_data(make([[arc] for arc in many[1:64]],63,many[1:64]),Set(many))
        @test boundary.empty_mask isa UInt64
        @test count_ones(foldl(|,values(boundary.arc_masks);init=boundary.empty_mask))==64
    end
    @testset "Prefix pruning independent of dominance and model builder" begin
        graph=DiGraph(4);cost=spzeros(Float64,4,4);risk=spzeros(Float64,4,4)
        for (u,v,c) in [(1,2,1.0),(2,3,1.0),(3,4,1.0),(1,4,10.0)];add_edge!(graph,u,v);cost[u,v]=c;end
        a=(1,2);b=(3,4);risk[a...]=-1.0;u=User("budget",1,4,risk,cost,nothing,nothing)
        for groups in ([[a],[b]],[[a,b]]), budget in (0,1,2), dominance in (false,true)
            h=HNDPwC(graph,[u],[a,b],Dict(a=>0.0,b=>0.0),nothing,groups,budget)
            raw,_,tr=_enumerate_user_paths_v2(u,h,10.0,time()+60;use_decision_arc_dominance=dominance,check_decision_budget=false)
            clean,_,tc=_enumerate_user_paths_v2(u,h,10.0,time()+60;use_decision_arc_dominance=dominance)
            feasible=budget>=length(groups)
            @test !tr && !tc
            @test getfield.(raw,:cost)==[3.0,10.0]
            @test getfield.(clean,:cost)==(feasible ? [3.0,10.0] : [10.0])
            for parallel in (false,true)
                inst,_,counts=build_hndp_path_instance(h,GurobiSolver();enumeration_time_limit=60.0,parallelize=parallel,use_decision_arc_dominance=dominance)
                m=_hndp_master_model(inst);set_silent(m);set_optimizer_attribute(m,"Threads",1);optimize!(m)
                @test termination_status(m)==MOI.OPTIMAL
                @test isapprox(objective_value(m),feasible ? -1.0 : 0.0;atol=1e-9)
                @test counts["budget"]==(feasible ? 2 : 1)
            end
        end
        # Ungrouped decisions are not terms of the group's opening budget.
        h=HNDPwC(graph,[u],[a,b],Dict(a=>0.0,b=>0.0),nothing,[[b]],0)
        short=User("ungrouped",1,2,risk,cost,nothing,nothing)
        p,_,to=_enumerate_user_paths_v2(short,h,2.0,time()+60;use_decision_arc_dominance=false)
        @test !to && length(p)==1 && only(p).arcs==[a]
        # A legitimately empty feasible path set must yield an infeasible MIP.
        g=DiGraph(2);add_edge!(g,1,2);c=spzeros(2,2);c[1,2]=1;r=spzeros(2,2)
        single=User("infeasible",1,2,r,c,nothing,nothing)
        hh=HNDPwC(g,[single],[(1,2)],Dict((1,2)=>0.0),nothing,[[(1,2)]],0)
        for parallel in (false,true),dominance in (false,true)
            inst,_,counts=build_hndp_path_instance(hh,GurobiSolver();enumeration_time_limit=60.0,parallelize=parallel,use_decision_arc_dominance=dominance)
            mm=_hndp_master_model(inst);set_silent(mm);set_optimizer_attribute(mm,"Threads",1);set_optimizer_attribute(mm,"DualReductions",0);optimize!(mm)
            @test counts["infeasible"]==0
            @test termination_status(mm)==MOI.INFEASIBLE
        end
    end
end
