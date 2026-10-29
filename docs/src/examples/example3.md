# Example 3: A Routing Follower Oracle

JuBiC's `AStarSolver` interface lets an application solve routing followers
without constructing a new JuMP MIP for every oracle call. The maintained
[example3.jl](example3.jl) uses the HNDP implementation with two users on a
small single-layer Sioux Falls network with nonnegative routing and operator
costs.

```sh
julia --project=. --threads=2 docs/src/examples/example3.jl
```

The master still uses Gurobi, so a working Gurobi license is required. The
follower route is selected with `subproblem_method=HNDP_SUBPROBLEM_ASTAR`.
No additional graph package beyond JuBiC's dependencies is needed.

A* requires nonnegative transition costs for each objective it evaluates. GBC
also asks the oracle for a master-contribution lower bound, so negative
operator/profit coefficients in competition or multimodal instances can make
that route unsupported. Use a MIP follower oracle for those instances, or a
generator configuration with nonnegative operator costs.

For implementing your own states, resource transitions, costs, and dominance
rules, follow the [HNDP A* tutorial](hndp/astar.md) and the
[subsolver interface](../sub_solvers.md). The application implementation is in
`examples/HNDP/hndp_astar_wrapper.jl`.
