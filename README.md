# JuBiC: Paper Reproduction Branch

This branch is dedicated to the paper [*Solving Multi-Follower Mixed-Integer
Bilevel Problems with Binary Linking Variables*](https://optimization-online.org/?p=28877).

---

## Getting Started

1. Install Julia and instantiate the project environment:

```powershell
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

2. Run the unit tests:

```powershell
julia --project=. test/runtests.jl
```

Note that you need an active Gurobi license.

---

## Reproducing Results

The [`benchmark`](./benchmark) directory contains the saved test instances and reproduction scripts.
The network-design benchmarks use alpha = 0.2 for length-constrained cases
and also include unrestricted cases. 

- [`benchmark/positive_risk`](./benchmark/positive_risk) contains the 60 saved
  Sioux Falls instances with nonnegative risk and the scripts
  [`run_positive_risk.jl`](./benchmark/positive_risk/run_positive_risk.jl) for the parallelization tests and
  [`run_positive_risk_additional.jl`](./benchmark/positive_risk/run_positive_risk_additional.jl) 
  for accelaration techniques tests.
- [`benchmark/negative_risk`](./benchmark/negative_risk) contains the 60 saved
  Sioux Falls instances with negative risk on decision arcs
  and [`run_negative_risk.jl`](./benchmark/negative_risk/run_negative_risk.jl) for reproducing accelaration techniques tests.
- [`benchmark/bobilib`](./benchmark/bobilib) contains the BOBILib experiments reproduction script.

You can run a network-design reproduction from the repository root with Julia using
the saved instances in the corresponding `instances` subfolder:

```powershell
julia --project=. benchmark/positive_risk/run_positive_risk.jl
julia --project=. benchmark/negative_risk/run_negative_risk.jl
```

The positive-risk main script accepts `all` (the default) or one of its named
parallelization modes. The supplementary positive-risk script runs the
additional GBC variants. All network-design scripts resume existing output
where possible and read the committed instances rather than generating new
random networks.

Some instances can exhibit numerical sensitivity in ConnectorLP cut
generation. When this occurs, the GBC solver supports manually setting a
problem-specific upper bound on the ConnectorLP variable `s` through the
`connector_s_bound` parameter.

For additional HNDP background, see the [benchmark documentation](./docs/src/examples/hndp/benchmarks.md)
and the [HNDP solver overview](./examples/HNDP/HNDP_SOLVER_OVERVIEW.md).

The BOBILib instances are not included in the repository. Download the
benchmark-interdiction set separately and pass its root directory explicitly:

```powershell
julia --project=. benchmark/bobilib/run_bobilib_nd_ws_pd_sub15.jl `
  --interdiction-root "C:\path\to\benchmark-interdiction"
```

[`run_bobilib_nd_ws_pd_sub15.jl`](./benchmark/bobilib/run_bobilib_nd_ws_pd_sub15.jl)
runs both MiBS and the `ND_WS_PD` GBC configuration sequentially. ND_WS_PD
uses one master thread, one worker, and 15 Gurobi threads in that worker. The
script resumes completed solver/instance pairs and writes both result types to
the combined CSV in `benchmark/bobilib`. Use `--source-csv` only when you want
to restrict the run to the instance names listed in an existing semicolon-
separated result CSV (e.g. if you want to continue a previous a run); otherwise all compatible instances below the supplied root are used.

---

## Disclaimer

JuBiC is still under active development. Results should be interpreted with
care, and computational results should be independently verified whenever
possible.

If you encounter difficulties or have questions, contact:
Vladimir Stadnichuk, `vladimir.stadnichuk@uni-kassel.de`.

---

## Citation

If you use this branch, please cite:

> Stadnichuk and Koster (2024), *Solving Multi-Follower Mixed-Integer
> Bilevel Problems with Binary Linking Variables*, Optimization Online.
> [Link](https://optimization-online.org/?p=28877)

---

## License

JuBiC is released as an open-source project. See the [LICENSE](./LICENSE) file
for details.
