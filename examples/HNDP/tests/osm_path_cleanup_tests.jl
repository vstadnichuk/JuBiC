using Test
if !isdefined(@__MODULE__, :_cleanup_strictly_dominated_hndp_paths)
    include("../models/hndp_model_generation.jl")
end

function hndp_path_cleanup_tests()
    a=(1,2);b=(5,6);c=(7,8)
    make(arcs,cost;risk=0.0,weight=0.0)=HNDPPath(arcs,cost,risk,weight)
    @testset "Complete-path strict dominance cleanup" begin
        expensive=make([a,(2,3),(3,5),b],10.0;risk=-100.0)
        cheap=make([a,(2,4),(4,5),b],4.0)
        fixed=make([(1,6)],12.0)
        @test _cleanup_strictly_dominated_hndp_paths([expensive,cheap,fixed],[a,b])==[cheap,fixed]
        @test _cleanup_strictly_dominated_hndp_paths([fixed,expensive,cheap],[a,b])==[fixed,cheap]
        @test _cleanup_strictly_dominated_hndp_paths([expensive,fixed,cheap],[a,b])==[fixed,cheap]
        # Equal follower costs remain, even with different revenue/resources/modes.
        tie=make([a,(2,9),(9,5),b],4.0;risk=-1.0,weight=5.0)
        @test _cleanup_strictly_dominated_hndp_paths([cheap,tie],[a,b])==[cheap,tie]
        @test _cleanup_strictly_dominated_hndp_paths([cheap,cheap],[a,b])==[cheap,cheap]
        alternative=make([c],2.0)
        @test _cleanup_strictly_dominated_hndp_paths([cheap,alternative],[a,b,c])==[cheap,alternative]
        subset=make([a],3.0;weight=10.0)
        @test _cleanup_strictly_dominated_hndp_paths([cheap,subset],[a,b])==[subset]
        always=make([(1,6)],1.0)
        @test _cleanup_strictly_dominated_hndp_paths([cheap,always],[a,b])==[always]
        @test isempty(_cleanup_strictly_dominated_hndp_paths(HNDPPath[],[a,b]))
        @test _cleanup_strictly_dominated_hndp_paths([cheap],[a,b])==[cheap]
        # More than 64 decision arcs must not collide or truncate masks.
        many=[(i,i+100) for i in 1:100]
        p=make([many[1],many[100]],5.0);q=make([many[100]],2.0);r=make([many[1]],1.0)
        @test _cleanup_strictly_dominated_hndp_paths([p,q,r],many)==[q,r]
        # A dominated witness may be replaced transitively; survivor order is stable.
        chain=[make([a,b],9.0),make([a],6.0),make(Tuple{Int,Int}[],3.0)]
        @test _cleanup_strictly_dominated_hndp_paths(chain,[a,b])==chain[3:3]
    end
    @testset "DFS discovery order and integration" begin
        graph=DiGraph(6);cost=spzeros(Float64,6,6);risk=spzeros(Float64,6,6)
        for (u,v,value) in [(1,2,2.0),(2,3,4.0),(3,5,4.0),(2,4,1.0),(4,5,1.0),(5,6,0.0),(1,6,12.0)]
            add_edge!(graph,u,v);cost[u,v]=value
        end
        risk[1,2]=-1.0
        user=User("toy",1,6,risk,cost,nothing,nothing)
        h=HNDPwC(graph,[user],[a,b],Dict(a=>0.0,b=>0.0),nothing,[[a],[b]],2)
        raw,_,to=_enumerate_user_paths_v2(user,h,12.0,time()+60;use_decision_arc_dominance=false)
        clean,_,tc=_enumerate_user_paths_v2(user,h,12.0,time()+60;use_decision_arc_dominance=true)
        @test !to && !tc
        @test getfield.(raw,:cost)==[10.0,4.0,12.0]
        @test getfield.(clean,:cost)==[4.0,12.0]
        @test getfield.(_cleanup_strictly_dominated_hndp_paths(raw,[a,b]),:arcs)==getfield.(clean,:arcs)
        for parallel in (false,true),budget in (0,1,2),dominance in (false,true)
            hh=HNDPwC(graph,[user],[a,b],Dict(a=>0.0,b=>0.0),nothing,[[a],[b]],budget)
            observed=Ref{Any}(nothing)
            observer=(event,index,user,data)->(event==:finished && (observed[]=data.paths);nothing)
            instance,_,counts=build_hndp_path_instance(hh,GurobiSolver();enumeration_time_limit=60.0,parallelize=parallel,use_decision_arc_dominance=dominance,enumeration_observer=observer)
            m=_hndp_master_model(instance);set_silent(m);set_optimizer_attribute(m,"Threads",1);optimize!(m)
            @test termination_status(m)==MOI.OPTIMAL
            @test isapprox(objective_value(m), (budget==2 ? -1.0 : 0.0); atol=1e-9)
            @test counts["toy"]==(budget<2 ? 1 : dominance ? 2 : 3)
            parallel && @test length(observed[])==counts["toy"]
        end
        # An already-expired enumeration remains explicitly incomplete.
        partial,_,timedout=_enumerate_user_paths_v2(user,h,12.0,time()-1;use_decision_arc_dominance=true)
        @test timedout && isempty(partial)
    end
end
