using CSV, JSON, Dates, Statistics
# Experiment worker: started by run_experiments.py, it receives one case at a time through
# <study>/runtime/<experiment>/command.json and writes the results into the case directory.
# SLP_RUNTIME (set by the controller) is the runtime folder of the current experiment group.
include(joinpath(@__DIR__,"paths.jl"))
include(joinpath(SRC, "experiments", "hndp_osm_experiment_runner.jl"))
include(joinpath(SRC, "experiments", "osm_full_od_mip_benchmark.jl"))
include(joinpath(@__DIR__,"benchmark_factory.jl"))
include(joinpath(@__DIR__,"path_shapes.jl"))
include(joinpath(@__DIR__,"audit_station.jl"))
const ROOT=get(ENV,"SLP_RUNTIME",joinpath(B,"runtime","worker"));mkpath(ROOT);const OUT=joinpath(B,"inputs")
function save(path,x)
 temp=path*".tmp"
 open(temp,"w") do io;JSON.print(io,x,2);end
 # Atomic rename without Base.rename's copy fallback, retrying reader locks.
 for attempt in 1:200
  code=ccall(:jl_fs_rename,Int32,(Cstring,Cstring),temp,path)
  code==0 && return
  attempt==200 && Base.uv_error("replace checkpoint",code)
  sleep(.005)
 end
end
function mode_mask(path,k)
 mask=UInt8(0)
 for (u,v) in path.arcs
  l=cld(u,k);l==cld(v,k) && (mask|=UInt8(1<<(l-1)))
 end
 mask
end
function modes_from_solution(mip,users,masks,directory,meta)
 lookup=Dict(string(u.uname)=>(i,u) for (i,u) in enumerate(users))
 mass=zeros(16,2);flow=zeros(length(users));revenue=0.0
 for (uid,weights) in hndp_path_solution_weights(mip)
  i,u=lookup[uid];class=endswith(uid,"_1") ? 1 : 2
  for (pi,z) in enumerate(weights)
   mask=masks[i][pi]
   mass[Int(mask)+1,class]+=z*u.mrisk.demand;flow[i]+=z
  end
 end
 rows=NamedTuple[]
 for c in 1:2,mask in 0:15
  powered=mask&14
  mode= count_ones(powered)>=2 ? "multi_powered" : powered==2 ? "car" : powered==4 ? "transit" : powered==8 ? "bike_sharing" : "walking"
  denominator=sum(mass[:,c]);allmass=sum(mass)
  push!(rows,(car_available=c==1,mode_mask=mask,mode=mode,any_mode_multimodal=count_ones(mask)>=2,multi_powered=count_ones(powered)>=2,demand=mass[mask+1,c],share_of_original_total=mass[mask+1,c],share_of_modeled_total=mass[mask+1,c]/allmass,share_within_class=mass[mask+1,c]/denominator))
 end
 CSV.write(joinpath(directory,"solution_modes.csv"),rows)
 Dict("modeled_demand"=>sum(mass),"maximum_user_path_sum_error"=>maximum(abs.(flow.-1)),"negative_path_flow_mass"=>sum(min(z,0.0) for z in mass),"modal_denominator"=>"original total demand = 1; conditional modeled/class shares also provided","intrazonal_owner_mass"=>meta["intrazonal_owner_mass"],"intrazonal_nonowner_mass"=>meta["intrazonal_nonowner_mass"])
end
function run_case(job)
 city=job["city"];k=Int(job["size"]);dominance=Bool(job["dominance"]);dir=job["directory"];mkpath(dir)
 state=joinpath(dir,"state.json")
 save(state,Dict("phase"=>"loading_instance","time"=>string(now())))
 agg=get(job,"aggregation_root",joinpath(OUT,"aggregation"))
 demand=get(job,"demand_root",joinpath(OUT,"demand"))
 generated=build_osm_hour_full_od_benchmark(city,k;aggregation_root=agg,network_root=joinpath(OUT,"city_networks_walk_connected"),demand_root=demand,station_construction_cost=0.0,candidate_count=Int(get(job,"candidate_count",20)),station_budget=Int(get(job,"station_budget",10)))
 h=generated.instance;save(joinpath(dir,"instance.json"),generated.metadata)
 n=length(h.users);masks=Vector{Vector{UInt8}}(undef,n);starts=zeros(n);lock=ReentrantLock()
 shapes=Vector{Vector{HourPathShape}}(undef,n);decisions=Set(h.edgeA)
 catalogue=Dict{Tuple{Bool,String,String,String},Tuple{Float64,Float64}}()
 complete_by_class=zeros(Int,2);mass_by_class=zeros(2)
 shape_trace=open(joinpath(dir,"path_shape_users.csv"),"w")
 shape_header_written=Ref(false)
 function checkpoint_shapes()
  save(joinpath(dir,"catalogue_checkpoint.json"),Dict("completed_users_by_class"=>complete_by_class,"completed_demand_by_class"=>mass_by_class,"histograms"=>hour_hist_rows(catalogue),"scope"=>"completed users only; equal weights within each user's path alternatives"))
 end
 trace=open(joinpath(dir,"user_events.csv"),"w")
 println(trace,"event,index,user,car_available,unix_time,case_elapsed_seconds,bound_seconds,enumeration_seconds,worker_seconds,path_count,timed_out,expanded_states")
 begin_time=time();deadline=begin_time+3600.0
 function phase(name;kwargs...)
  record=Dict{String,Any}("phase"=>name,"benchmark_start_unix"=>begin_time,"time"=>string(now()),"elapsed_seconds"=>time()-begin_time)
  for (key,value) in kwargs;record[string(key)]=value;end
  save(state,record)
 end
 phase("enumerating")
 finished=Ref(0)
 function observer(event,index,user,data)
  stamp=time();uname=isnothing(user) ? "" : string(user.uname);car=endswith(uname,"_1")
  if event==:finished
   shapes[index]=[hour_path_shape(path,k,decisions) for path in data.paths]
   masks[index]=[s.movement for s in shapes[index]]
  end
  Base.lock(lock) do
   if event==:started
    starts[index]=stamp
    println(trace,"started,$index,$uname,$car,$stamp,$(stamp-begin_time),,,,,,")
   elseif event==:finished
    println(trace,"finished,$index,$uname,$car,$stamp,$(stamp-begin_time),$(data.bound_runtime),$(data.enumeration_runtime),$(data.runtime),$(length(data.paths)),$(data.timed_out),")
    row=hour_user_shape_row(index,user,shapes[index],!data.timed_out)
    shape_header_written[] || println(shape_trace,join(keys(row),","));shape_header_written[]=true
    println(shape_trace,join(values(row),","))
    if !data.timed_out
     cls=car ? 1 : 2;complete_by_class[cls]+=1;mass_by_class[cls]+=user.mrisk.demand
     hour_merge_shapes!(catalogue,shapes[index],car,user.mrisk.demand)
    end
    finished[]+=1
    finished[]%1000==0 && checkpoint_shapes()
    finished[]%100==0 && phase("enumerating";users_finished=finished[],users_total=n)
   elseif event==:progress
    println(trace,"progress,$index,$uname,$car,$stamp,$(stamp-begin_time),,,$(data.elapsed_seconds),$(data.paths_found),,$(data.expanded_states)")
   elseif event==:assembly_started
    phase("assembling";users_finished=finished[],users_total=n)
   elseif event==:assembled
    index%100==0 && phase("assembling";users_assembled=index,users_total=n)
   end
   finished[]%100==0 && (flush(trace);flush(shape_trace))
  end
 end
 result=Dict{String,Any}("formulation"=>"station_path","dominance_rule"=>"stations","city"=>city,"size"=>k,"dominance"=>dominance,"users"=>n,"threads"=>Threads.nthreads(),"time_budget_seconds"=>3600,"construction_cost"=>0,"candidate_count"=>Int(get(job,"candidate_count",20)),"station_budget"=>Int(get(job,"station_budget",10)))
 try
  solver=GurobiSolver(silent=false)
  instance,enumtime,counts=build_hndp_station_path_instance(h,solver;enumeration_time_limit=max(0.0,deadline-time()),use_station_dominance=dominance,eliminate_single_path_users=true,enumeration_observer=observer)
  result["enumeration_seconds"]=enumtime;result["build_seconds"]=time()-begin_time
  result["enumerated_users"]=length(counts);result["path_count"]=sum(values(counts));result["enumeration_complete"]=length(counts)==n && all(>(0),values(counts))
  if !result["enumeration_complete"] || time()>=deadline
   result["status"]=result["enumeration_complete"] ? "time_limit_during_assembly" : "time_limit_during_enumeration"
   result["model_seconds"]=time()-begin_time;save(joinpath(dir,"result.json"),result);return
  end
  mip=_hndp_master_model(instance)
  merge!(result,hndp_path_reduction_metadata(mip))
  result["variables"]=num_variables(mip);result["constraints"]=num_constraints(mip;count_variable_in_set_constraints=false)
  save(joinpath(dir,"build.json"),result)
  phase("solving")
  set_optimizer_attribute(mip,"LogFile",joinpath(dir,"gurobi.log"));set_optimizer_attribute(mip,"TimeLimit",max(.01,deadline-time()));set_optimizer_attribute(mip,"Threads",Threads.nthreads());set_optimizer_attribute(mip,"Seed",20261102);set_optimizer_attribute(mip,"SoftMemLimit",14.0)
  best_incumbent=Ref(Inf)
  function callback(cbdata,cb_where)
   if cb_where==Gurobi.GRB_CB_MIPSOL
    obj=Ref{Float64}(0);Gurobi.GRBcbget(cbdata,cb_where,Gurobi.GRB_CB_MIPSOL_OBJ,obj)
    if obj[]<best_incumbent[]-1e-9
     Gurobi.load_callback_variable_primal(cbdata,cb_where)
     save(joinpath(dir,"incumbent.json"),Dict("objective"=>obj[],"elapsed_seconds"=>time()-begin_time,"open_candidate_ids"=>[g[1][1] for g in h.decision_groups if callback_value(cbdata,mip[:x][g[1]])>.5]))
     best_incumbent[]=obj[]
    end
   end
   time()>=deadline && Gurobi.GRBterminate(cbdata.model)
  end
  MOI.set(mip,Gurobi.CallbackFunction(),callback)
  solvebegin=time();optimize!(mip)
  result["solve_wall_seconds"]=time()-solvebegin;result["model_seconds"]=time()-begin_time
  result["status"]=string(termination_status(mip));result["objective_bound"]=objective_bound(mip);result["nodes"]=node_count(mip)
  if has_values(mip)
   result["objective"]=objective_value(mip);result["gap"]=relative_gap(mip)
   result["open_candidate_ids"]=[g[1][1] for g in h.decision_groups if value(mip[:x][g[1]])>.5]
   result["opened_stations"]=length(result["open_candidate_ids"])
   save(joinpath(dir,"result.json"),result)
   phase("postprocessing")
   result["modal_statistics"]=modes_from_solution(mip,h.users,masks,dir,generated.metadata)
   solution_hist=Dict{Tuple{Bool,String,String,String},Tuple{Float64,Float64}}()
   weights=hndp_path_solution_weights(mip)
   for (i,u) in enumerate(h.users)
    hour_merge_shapes!(solution_hist,shapes[i],endswith(string(u.uname),"_1"),u.mrisk.demand;group="solution",weights=weights[string(u.uname)])
   end
   CSV.write(joinpath(dir,"solution_path_shapes.csv"),hour_hist_rows(solution_hist))
  end
  save(joinpath(dir,"result.json"),result)
 catch error
  result["status"]="error";result["error"]=sprint(showerror,error,catch_backtrace());result["model_seconds"]=time()-begin_time
  save(joinpath(dir,"result.json"),result);showerror(stdout,error,catch_backtrace());println();flush(stdout)
 finally
  checkpoint_shapes();CSV.write(joinpath(dir,"catalogue_path_shapes.csv"),hour_hist_rows(catalogue))
  close(shape_trace);close(trace);phase("done")
 end
end
println("RUNTIME_IMPORTED_WAITING_FOR_INPUTS");flush(stdout)
while !isfile(joinpath(OUT,"aggregation/aachen/COMPLETE.json"));sleep(2);end
# Validate the instrumented threaded implementation before accepting jobs.
@assert hour_shape_smoke_tests()
isdefined(@__MODULE__,:station_cost_fixture) || include(joinpath(SRC, "tests", "osm_station_decision_tests.jl"))
include(joinpath(SRC, "tests", "osm_station_path_tests.jl"))
# A failing top-level testset throws, so the worker does not become ready.
@testset "station path model" begin hndp_station_path_tests() end
let h=station_cost_fixture(2),events=Tuple{Symbol,Int}[],guard=ReentrantLock()
 observer=(event,index,user,data)->Base.lock(guard) do;push!(events,(event,index));end
 _,_,counts=build_hndp_station_path_instance(h,GurobiSolver();enumeration_time_limit=30.0,enumeration_observer=observer)
 @assert length(counts)==1 && count(e->e[1]==:finished,events)==1 && any(e->e[1]==:assembly_started,events)
end
let out=OUT
 args=(aggregation_root=joinpath(out,"aggregation"),network_root=joinpath(out,"city_networks_walk_connected"),demand_root=joinpath(out,"demand"))
 old=build_osm_full_od_benchmark("aachen",50;args...)
 new=build_osm_hour_full_od_benchmark("aachen",50;args...)
 @assert old.metadata["candidate_ids"]==new.metadata["candidate_ids"]
 for cls in 1:2
  @assert old.instance.users[cls].mcost.fixed_distances[1:50,1:50]==new.instance.users[cls].mcost.fixed_distances
 end
 save(joinpath(ROOT,"factory_validation.json"),Dict("passed"=>true,"method"=>"parallel physical-source Dijkstra equals Floyd-Warshall on both Aachen50 class profiles","threads"=>Threads.nthreads()))
end
let out=OUT
 args=(aggregation_root=joinpath(out,"aggregation"),network_root=joinpath(out,"city_networks_walk_connected"),demand_root=joinpath(out,"demand"))
 small=build_osm_hour_full_od_benchmark("aachen",50;args...,candidate_count=20,station_budget=5)
 large=build_osm_hour_full_od_benchmark("aachen",50;args...,candidate_count=40,station_budget=30)
 @assert issubset(Set(small.metadata["candidate_ids"]),Set(large.metadata["candidate_ids"]))
 @assert length(large.instance.decision_groups)==40
 @assert large.metadata["decision_budget"]==30 && small.metadata["decision_budget"]==5
 @assert large.metadata["nusers"]==small.metadata["nusers"]==4900
 save(joinpath(ROOT,"candidate_validation.json"),Dict("passed"=>true,"nested_candidates"=>true,"same_users"=>true,"budgets"=>[5,30]))
end
save(joinpath(ROOT,"worker_ready.json"),Dict("ready"=>true,"time"=>string(now())))
println("WORKER_READY");flush(stdout)
function serve_jobs()
 last_job_id=""
 while true
  commandfile=joinpath(ROOT,"command.json")
  if !isfile(commandfile);sleep(.2);continue;end
  job=try
   JSON.parsefile(commandfile)
  catch err
   if err isa SystemError || err isa Base.IOError
    sleep(.1);continue
   end
   rethrow()
  end
  job_id=String(job["id"])
  if job_id==last_job_id;sleep(.2);continue;end
  get(job,"stop",false) && break
  last_job_id=job_id
  println("CASE_START $job_id");flush(stdout)
  if get(job,"kind","")=="audit"
   run_audit(String(job["city"]),50;maxseconds=3600.0)
  else
   run_case(job)
  end
  GC.gc();GC.gc()
  save(joinpath(ROOT,"done.json"),Dict("id"=>job_id,"time"=>string(now())))
  println("CASE_DONE $job_id");flush(stdout)
 end
end
serve_jobs()
