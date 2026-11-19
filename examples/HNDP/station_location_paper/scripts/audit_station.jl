# Dominance audit for the station-based rerun.
# For every user of the 50-location instances, compare
#   raw:     all simple paths within the fixed-network cost bound (no dominance, no budget rule),
#   station: station dominance + budget rule + station-based complete-path cleanup (new default),
#   legacy:  decision-arc dominance + budget rule (previous implementation, for reference only),
# and certify that every removed budget-feasible raw path has a retained replacement that requires
# a subset of its stations and has no greater follower cost and no worse leader contribution.
stationsubset(a,b)=issubset(a,b)
function run_audit(city,k;maxseconds=3600.0)
 dir=joinpath(ROOT,"audit_$(city)_$(k)");mkpath(dir)
 save(joinpath(dir,"state.json"),Dict("phase"=>"loading_instance","city"=>city,"size"=>k,"time"=>string(now())))
 gen=build_osm_hour_full_od_benchmark(city,k;aggregation_root=joinpath(OUT,"aggregation"),network_root=joinpath(OUT,"city_networks_walk_connected"),demand_root=joinpath(OUT,"demand"),station_construction_cost=0.0)
 h=gen.instance;data=hndp_station_data(h);budget=Int(h.decision_budget)
 rows=Vector{Any}(undef,length(h.users));detailfile=open(joinpath(dir,"paths_and_witnesses.jsonl"),"w");done=Threads.Atomic{Int}(0);lk=ReentrantLock();started=time();deadline=started+maxseconds
 Threads.@threads :dynamic for i in eachindex(h.users)
  u=h.users[i];bound=u.mcost.fixed_distances[u.origin,u.destination]
  raw,rt,rto=_enumerate_user_paths_station_dominance(u,h,bound,deadline;use_station_dominance=false,check_decision_budget=false,station_data=data)
  dom,dt,dto=_enumerate_user_paths_station_dominance(u,h,bound,deadline;use_station_dominance=true,check_decision_budget=true,station_data=data)
  leg,lt,lto=_enumerate_user_paths_v2(u,h,bound,deadline;use_decision_arc_dominance=true,check_decision_budget=true)
  signature(p)=Tuple(p.arcs)
  rawsig=Set(signature(p) for p in raw);domsig=Set(signature(p) for p in dom)
  rs=[hndp_path_stations(p,data) for p in raw];ds=[hndp_path_stations(p,data) for p in dom]
  # Leader key: lower is better for the optimistic upper level (risk = demand-weighted gamma).
  key(p)=(p.cost,p.risk)
  removed=findall(p->!(signature(p) in domsig),raw)
  witnesses=Any[length(rs[q])>budget ? "budget_infeasible" : findfirst(j->stationsubset(ds[j],rs[q]) && key(dom[j])<=key(raw[q]),eachindex(dom)) for q in removed]
  failures=count(isnothing,witnesses)
  not_in_raw=count(p->!(signature(p) in rawsig),dom)
  infeasible=count(w->w=="budget_infeasible",witnesses)
  row=(user=string(u.uname),car_available=endswith(string(u.uname),"_1"),demand=u.mrisk.demand,raw_paths=length(raw),dominance_paths=length(dom),legacy_arc_dominance_paths=length(leg),raw_seconds=rt,dominance_seconds=dt,legacy_seconds=lt,raw_timeout=rto,dominance_timeout=dto,legacy_timeout=lto,dominance_not_in_raw=not_in_raw,removed_paths=length(removed),budget_infeasible_paths=infeasible,certificate_failures=failures)
  rows[i]=row
  pack(paths,ss)=[Dict("arcs"=>p.arcs,"cost"=>p.cost,"risk"=>p.risk,"stations"=>collect(ss[j])) for (j,p) in enumerate(paths)]
  detail=Dict("user"=>string(u.uname),"raw"=>pack(raw,rs),"dominance"=>pack(dom,ds),"removed_raw_indices"=>removed,"station_witness_dom_indices"=>witnesses,"legacy_arc_dominance_paths"=>length(leg))
  lock(lk) do;JSON.print(detailfile,detail);println(detailfile);end
  n=Threads.atomic_add!(done,1)+1
  if n%100==0
   lock(lk) do
    save(joinpath(dir,"state.json"),Dict("phase"=>"enumerating_and_certifying","benchmark_start_unix"=>started,"city"=>city,"size"=>k,"completed_users"=>done[],"total_users"=>length(h.users),"elapsed_seconds"=>time()-started,"time"=>string(now())))
   end
  end
 end
 CSV.write(joinpath(dir,"users.csv"),rows)
 close(detailfile)
 fields=[:raw_paths,:dominance_paths,:legacy_arc_dominance_paths,:raw_seconds,:dominance_seconds,:legacy_seconds,:removed_paths,:dominance_not_in_raw,:budget_infeasible_paths,:certificate_failures]
 summary=Dict{String,Any}(string(f)=>sum(getproperty(r,f) for r in rows) for f in fields)
 summary["city"]=city;summary["size"]=k;summary["users"]=length(rows);summary["dominance_rule"]="stations"
 summary["incomplete_users"]=count(r->r.raw_timeout || r.dominance_timeout,rows);summary["wall_seconds"]=time()-started;summary["zero_demand_users"]=count(r->r.demand==0,rows)
 summary["one_path_users_raw"]=count(r->r.raw_paths==1,rows);summary["one_path_users_dominance"]=count(r->r.dominance_paths==1,rows)
 summary["users_station_below_legacy"]=count(r->r.dominance_paths<r.legacy_arc_dominance_paths,rows)
 summary["instance"]=gen.metadata
 save(joinpath(dir,"summary.json"),summary)
 println("AUDIT_DONE ",JSON.json(summary));flush(stdout);GC.gc();summary
end
