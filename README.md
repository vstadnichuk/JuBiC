# JuBiC: Reproduction Branch for "Bilevel Feasible Paths for the Station Location Problem"

This branch contains the code and the experiment scripts of the manuscript

> V. Stadnichuk, *Bilevel Feasible Paths for the Station Location Problem*.

The manuscript is currently under review. The draft can be provided on request
(contact below). This README describes how to reproduce the computational results
of the paper. The multimodal network instance generator used in the paper is described in
[`examples/HNDP/INSTANCE_GENERATOR.pdf`](examples/HNDP/INSTANCE_GENERATOR.pdf). For JuBiC in general, see the `main` branch and the documentation in
[`docs/src/index.md`](./docs/src/index.md).

---

## Where the paper's methods are implemented

| Paper | Code |
|---|---|
| Path-based reformulation (PBR) with station-wise linking and the big-*M* values of Lemma 5 | [`examples/HNDP/models/hndp_station_path_model.jl`](examples/HNDP/models/hndp_station_path_model.jl), `build_hndp_station_path_instance` |
| Enumeration of bilevel feasible paths (Algorithm 1: cost bound, budget rule, station dominance, cleanup) | same file, `_enumerate_user_paths_station_dominance` and `_cleanup_station_dominated_hndp_paths` |
| Tests (station dominance, budget pruning, cleanup, model size, LP strength) | [`examples/HNDP/tests/osm_station_path_tests.jl`](examples/HNDP/tests/osm_station_path_tests.jl) |
| OSM/GTFS instance generator (Section 5 and [`INSTANCE_GENERATOR.pdf`](examples/HNDP/INSTANCE_GENERATOR.pdf)) | [`examples/HNDP/OSM_NETWORKS.md`](examples/HNDP/OSM_NETWORKS.md) and [`examples/HNDP/DEMAND_AND_AGGREGATION.md`](examples/HNDP/DEMAND_AND_AGGREGATION.md) |
| Experiments, analysis and figures (Section 6) | [`examples/HNDP/station_location_paper/`](examples/HNDP/station_location_paper/) |

---

## Requirements

- Julia 1.11.7 and a Gurobi license; the paper used Gurobi 13.0.3. Instantiate the
  project environment of this repository once:

  ```powershell
  julia --project=. -e "using Pkg; Pkg.instantiate()"
  ```

- Python 3 with the packages in
  [`station_location_paper/requirements.txt`](examples/HNDP/station_location_paper/requirements.txt)
  (pyosmium for the OSM extracts, NumPy for the transit lengths, ReportLab for the figures).
- The paper's runs used 8 threads, about 32 GiB of RAM and Windows 11. The experiment
  controller enforces a memory guard of 16 GiB and a time limit of one hour per instance.
  Running times depend on the hardware, so expect them to differ from the paper.
- Disk space: about 1 GB for the downloaded sources, about 0.3 GB for the generated instances, about 30 GB for the all-pairs
  reference caches that are built during instance preparation, and up to 15 GB for the
  raw experiment logs.

All scripts read their parameters from
[`station_location_paper/settings.json`](examples/HNDP/station_location_paper/settings.json)
and write into a study folder, by default `examples/HNDP/station_location_paper/study/`
(ignored by git). Set `SLP_STUDY_DIR` to use another location, and `JULIA` if `julia` is
not on the `PATH`. All commands below are run from `examples/HNDP/station_location_paper/`.

---

## Step 0: Check the implementation

```powershell
julia --project=../../.. --threads=8 scripts/run_tests.jl
```

The experiment worker runs the same tests on start-up and does not start if one fails.

## Step 1: Generate the city networks

```powershell
python scripts/download_sources.py
python scripts/prepare_sources.py
julia --project=../../.. --threads=8 scripts/generate_networks.jl
python ../network_generation/osm_transit_arc_lengths.py --raw-root study/source_data --network-root study/inputs/city_networks_walk_connected --output-root study/inputs/geometry/transit_lengths aachen cologne heidelberg karlsruhe nuremberg dresden
```

1. `download_sources.py` downloads the BBBike road extracts, the GTFS feeds and the municipal
   boundaries listed in [`osm_city_sources.json`](examples/HNDP/sources/osm_city_sources.json) to
   `study/source_data/sources/<city>/`.
2. `prepare_sources.py` cuts the roads and the transit trips to the study rectangle of each city
   (municipal bounding box plus 5 km).
3. `generate_networks.jl` builds the layered networks with the JuBiC generator
   (Section 5 and `INSTANCE_GENERATOR.pdf`) and recovers the physical road lengths.
4. `osm_transit_arc_lengths.py` recovers the physical lengths of the transit arcs from the
   GTFS shapes.

The paper uses the OSM data of 24 October 2026 and the GTFS feeds downloaded on
28 October 2026. Note that the registry records the checksums of these files. The providers only
serve their current data, so a new download yields networks that reflect the current OSM
and timetable data. You will find the records which files differ
from the paper's snapshot in `study/source_data/download_report.json`.

## Step 2: Prepare the instances

```powershell
julia --project=../../.. --threads=8 scripts/full_network_city_profiles.jl
julia --project=../../.. --threads=8 scripts/prepare_instances.jl
```

`full_network_city_profiles.jl` computes the full-network demand and the modal split with
all stations open (Table 2 of the paper, Table 2 of `INSTANCE_GENERATOR.pdf`). `prepare_instances.jl` then selects the
retained locations (50 to 400) and the nested potential-station sets, and aggregates the
demand (Section 5). Both scripts build their all-pairs reference caches in `study/work/` on
the first run and reuse them afterwards.

## Step 3: Run the experiments

```powershell
python scripts/run_experiments.py aggregation
python scripts/run_experiments.py scaling
```

| Command | Paper | Cases | Case time in the paper |
|---|---|---|---|
| `aggregation` | Sections 6.1 and 6.2 (Table 3, Table EC.1, Figures 1 and 2) | 6 dominance audits at 50 locations; 36 instances with 50 to 400 locations, 20 potential stations, budget 10 | about 2 h |
| `scaling` | Section 6.3 (Figures 3 and 4) | 300 locations: 20/80/140 potential stations with 25/50/75 % budgets, and 20 to 300 potential stations with a 25 % budget | about 17 h |

The controller starts one Julia worker and passes it one case at a time. The worker records
per-user enumeration times, Gurobi logs and the solution in
`study/raw/<experiment>_cases/<case>/`. Each case has a total time limit of one hour for
enumeration, model construction and optimization. Larger instances of a city are skipped
after a failure, as described in Section 6. The controller can be stopped at any time;
completed cases are reused when it is started again.

## Step 4: Evaluate the designs on the full networks

```powershell
julia --project=../../.. --threads=8 scripts/replay_designs_full_network.jl
python scripts/resolve_full_network_ties.py
```

The first script evaluates the optimal designs at 50 and 300 locations on the full
networks (modal split "F" in Figure 2). The second script resolves modal ties of the
full-network profiles of Step 2 (Table 2 of the paper, Table 2 of `INSTANCE_GENERATOR.pdf`).

## Step 5: Tables and figures

```powershell
python scripts/analyze_aggregation.py
python scripts/analyze_scaling.py
python scripts/network_topology.py
python scripts/figures_aggregation.py
python scripts/figures_scaling.py
python scripts/figure_network_maps.py
```

| Paper | Source in `study/` |
|---|---|
| Table 1 (network topology) | `data/network_topology.csv` (scope `full`) |
| Table 2 and Table 2 of `INSTANCE_GENERATOR.pdf` (users, modal split with all stations open) | `data/full_network_city_profiles_resolved.csv` |
| Table 3 and Table EC.1 of the e-companion (paths and enumeration times with and without dominance) | `data/aggregation/dominance50.csv` |
| Figure 1 (paths per user vs. aggregation) | `figures/aggregation_paths_per_user.pdf` |
| Figure 2 (modal split, aggregated vs. full network) | `figures/aggregation_modal.pdf` |
| Figure 3 (runtime grid) | `figures/candidate_runtime.pdf` |
| Figure 4 (station scaling) | `figures/station_scaling.pdf`, data in `data/scaling/station_scaling.csv` |
| Section 6.3 statistics (enumeration share, per-user times) | `data/scaling/solver_cases.csv`, `data/scaling/user_enumeration_summary.csv` |
| Figure 1 of `INSTANCE_GENERATOR.pdf` (network maps) | `figures/network_maps.pdf` |

---

## Contact

Vladimir Stadnichuk, `vladimir.stadnichuk@uni-kassel.de`

## Citation

Until the paper is published, please cite JuBiC as

> Stadnichuk and Koster (2024), *Solving Multi-Follower Mixed-Integer Bilevel
> Problems with Binary Linking Variables*, Optimization Online.
> [Link](https://optimization-online.org/?p=28877)

## License

JuBiC is released as an open-source project. See the [LICENSE](./LICENSE) file for details.
