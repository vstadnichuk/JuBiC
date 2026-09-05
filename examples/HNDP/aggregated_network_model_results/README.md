# Aggregated-network model comparison

This folder preserves the final aggregated-network benchmark used during the
HNDP experiments on Cologne, Heidelberg, and Karlsruhe.

The main result is [`network_balanced_hardness_sweep_fixed.csv`](./network_balanced_hardness_sweep_fixed.csv).
It compares six models on spatially balanced networks with 50--300 retained
stations, 100--1,000 users, a 50% station budget, and a 600-second limit:

- strong duality with fixed-network and `n-1` big-M values;
- Benders-like cuts with A* and the same two big-M variants;
- path enumeration with and without decision-arc dominance.

The run used the Cologne, Heidelberg, and Karlsruhe OSM/GTFS-derived networks
and seed `20260829`. The required source data are intentionally kept outside
the repository under the ignored `data/osm/` directory.

To reproduce the model sweep from the repository root:

```powershell
julia --project=. examples/HNDP/aggregated_network_model_results/reproduce_aggregated_model_runs.jl
```

The reproduction script writes a new CSV to `tmp_compare/` and does not modify
the archived result in this folder. It uses all available Julia threads and
requires a working Gurobi installation for the MIP-based models.

The helper files beside the script are kept together so the reproduction
entry point remains self-contained relative to this folder.
