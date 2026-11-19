using Test
if !isdefined(@__MODULE__, :station_cost_fixture)
    include("osm_station_decision_tests.jl")
end

"""Run with Julia --threads=4 to cover concurrent path and MIP-bound workers."""
function osm_parallel_path_regression_tests()
@testset "Zero-cost stations and parallel enumeration" begin
 @test Threads.nthreads()>=1
 for dominance in (true,false),boundmethod in (:labeling,:mip),budget in (1,2)
  base=station_cost_fixture(budget)
  users=[User("u$(i)",1,2,i*base.users[1].mrisk,base.users[1].mcost,nothing,nothing) for i in 1:8]
  h=HNDPwC(base.mygraph,users,base.edgeA,Dict(a=>0.0 for a in base.edgeA),nothing,base.decision_groups,budget)
  outcomes=[];counts=[]
  for parallel in (false,true)
   instance,elapsed,c=build_hndp_path_instance(h,GurobiSolver();enumeration_time_limit=60.0,parallelize=parallel,use_decision_arc_dominance=dominance,fixed_path_bound_method=boundmethod)
   @test 0<=elapsed<60 && length(c)==8
   m=_hndp_master_model(instance);set_silent(m);set_optimizer_attribute(m,"Threads",1);optimize!(m)
   @test termination_status(m)==MOI.OPTIMAL
   @test objective_value(m)≈(budget==2 ? -1440.0 : 0.0) atol=1e-6
   @test value(m[:construction_cost_var])≈0 atol=1e-8
   @test sum(value(m[:x][first(g)]) for g in h.decision_groups)<=budget+1e-6
   @test all(isapprox(value(m[:x][first(g)]),value(m[:x][last(g)]);atol=1e-6) for g in h.decision_groups)
   push!(outcomes,objective_value(m));push!(counts,c)
  end
  @test outcomes[1]≈outcomes[2] atol=1e-6
  @test counts[1]==counts[2]
 end
end
end
