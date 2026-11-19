# Station-based path model and path enumeration

This guide describes the formulation solved in all experiments of the paper, the
path-based reformulation (PBR) with bilevel feasible paths, and how one experiment
case is run. The code is
[hndp_station_path_model.jl](models/hndp_station_path_model.jl).

## Stations

`hndp_station_data(hndp)` turns the decision groups of an instance into stations:
each group (bike entry and return arc of one location) is one station with one
binary opening variable, indexed by its first arc. Only these stations count
towards the opening budget. `hndp_path_stations(path, data)` returns the set
`S(p)` of stations a path uses.

## Formulation

`build_hndp_station_path_instance(hndp, solver; enumeration_time_limit, ...)`
builds one MIP. Per user `u` with enumerated path set `P_u`:

- path variables `λ_p ≥ 0` with `Σ_p λ_p = 1` (continuous; an optimal vertex picks
  one path);
- station-wise linking `Σ_{p : s ∈ S(p)} λ_p ≤ x_s` for every station used by a path;
- the follower cost `θ_u = Σ_p c_p λ_p` and, for every path `q`, the optimality
  constraint `θ_u ≤ c_q + M_q Σ_{s ∈ S(q)} (1 − x_s)` with
  `M_q = max_{p ∈ P_u} c_p − c_q` (Lemma 5 of the paper).

The budget constraint `Σ_s x_s ≤ B` limits the open stations. The leader minimizes
the negative demand-weighted entry revenue of the chosen paths plus construction
costs (zero in the experiments). Fixed arcs are omitted from teh model.

Users whose only path uses no station are removed from the MIP and their constant
contribution is added to the objective (`eliminate_single_path_users=true`).
`hndp_path_solution_weights(model)` restores their path weights for reporting.

## Enumeration of bilevel feasible paths

`_enumerate_user_paths_station_dominance` enumerates, per user, the simple paths
from origin to destination by depth-first search (Algorithm 1 of the paper):

1. **Cost bound.** The cost of the cheapest path without any station is an upper
   bound on the follower's cost. A partial path is
   discarded once its cost, plus the cheapest completion to the destination, exceeds
   this bound. Completion costs come from one reverse Dijkstra per destination and
   user class, shared by all users with that destination.
2. **Budget rule.** A partial path is discarded once it needs more stations than the
   budget allows.
3. **Station dominance.** At every node, a partial path is discarded if another
   partial path to that node needs a subset of its stations and is cheaper, or
   equally expensive with a leader contribution that is not worse.
4. **Cleanup.** After the search, every complete path for which a strictly cheaper
   complete path with a subset of its stations exists is removed. 

Users are enumerated in parallel on all Julia threads; the MIP is then assembled
serially. If the time limit is reached during enumeration, the builder returns an
empty model, and the case is reported as a timeout.

`use_station_dominance=false` keeps only the cost bound and the budget rule. The
dominance audit additionally compares against a raw enumeration with neither rule;
see below.

## One experiment case

`run_experiments.py` starts one Julia worker (`station_location_paper/scripts/worker.jl`),
which first runs the station-model tests and refuses to start if one fails. For
each case it:

1. builds the instance with `build_osm_hour_full_od_benchmark`
   (see [DEMAND_AND_AGGREGATION.md](DEMAND_AND_AGGREGATION.md));
2. builds the model with `build_hndp_station_path_instance`, recording per-user
   bound and enumeration times and path counts (`user_events.csv`);
3. solves it with Gurobi on 8 threads (seed 20261102) within the remaining time;
4. writes the result, the open stations and the modal split of the solution.

The one-hour limit covers enumeration, assembly and optimization. Results are in
`study/raw/<experiment>_cases/<case>/` (`result.json`, `gurobi.log`,
`solution_modes.csv`, path statistics).

## Dominance audit

`audit_station.jl` (run as part of `run_experiments.py aggregation`) enumerates every
user of the 50-location instances three times:

- **raw:** all simple paths within the cost bound, without budget rule or dominance;
- **station:** the enumeration above;
- **legacy:** dominance on individual decision arcs instead of stations, for reference.

It checks that every path removed by station dominance has a retained witness that is
cheaper and needs a subset of its stations. Table 3 of the paper and the dominance times
in Table EC.1 of the e-companion come from this audit.

## Tests

[osm_station_path_tests.jl](tests/osm_station_path_tests.jl) covers station
dominance, the budget rule, the cleanup, the model size and the LP bound against
toy instances. Run them with `station_location_paper/scripts/run_tests.jl`.
