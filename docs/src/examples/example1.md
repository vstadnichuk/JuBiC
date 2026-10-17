# Example 1: A Small Bilevel Model

The [Getting Started](../getting_started.md) tutorial gives the mathematics,
model construction, and GBC/MiBS solve routes. The standalone
[example1.jl](example1.jl) builds that model and solves it with GBC using the
current `SubSolverJuMP` constructor.

From the repository root, with a Gurobi license:

```sh
julia --project=. docs/src/examples/example1.jl
```

The expected objective is 1, up to numerical tolerance. Output is created in a
unique repository-local scratch folder; no pre-existing logs are required.
