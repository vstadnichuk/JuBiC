# JuBiC: Paper Reproduction Branch

This branch is dedicated to the paper [*Solving Multi-Follower Mixed-Integer
Bilevel Problems with Binary Linking Variables*](https://optimization-online.org/?p=28877).
It preserves the JuBiC implementation of the decomposition methods described
in the paper and provides the examples and tests needed to work with them.

If you want to work with JuBiC beyond this paper, check out the latest version
of the [`main`](https://github.com/vstadnichuk/JuBiC/tree/main) branch. The
general documentation starts at [`docs/src/index.md`](./docs/src/index.md).

To build the rendered HTML documentation locally, run:

```powershell
julia --project=. docs/make.jl
```

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

3. Read the documentation for the problem class and implemented solvers:

- [Problem formulation](./docs/src/index.md#problem-formulation)
- [Getting started](./docs/src/getting_started.md)
- [Generalized Benders cuts](./docs/src/solvers/gbc.md)
- [Benders-like cuts](./docs/src/solvers/blc.md)

The [`examples`](./examples) directory contains runnable demonstrations and
the HNDP application used to exercise the JuBiC solvers.

---

## Reproducing Results

The paper studies decomposition algorithms for multi-follower mixed-integer
bilevel problems with binary linking variables. The documented examples and
the HNDP benchmark tooling provide the reproducible starting point for these
experiments.

For HNDP experiments, see the [benchmark documentation](./docs/src/examples/hndp/benchmarks.md)
and the [HNDP solver overview](./examples/HNDP/HNDP_SOLVER_OVERVIEW.md).
Use the supplied scripts and manifests from the repository root with the
project environment enabled.

Exact runtimes and solver output can vary with the Julia version, CPU, number
of threads, Gurobi installation, and license settings.

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
