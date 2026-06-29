# `AStarSolver`

`AStarSolver` is JuBiC's generic wrapper for adding labeling or A*-style
subsolvers to JuBiC applications. It is primarily an interface layer: the
application provides the state representation, transition logic, costs,
dominance rule, and start/end conditions, while JuBiC can call the resulting
labeling routine through the standard subsolver interface.

This route is useful when a follower oracle is better implemented as a custom
labeling algorithm than as a JuMP model. The HNDP example uses this approach for
path-based follower solves.

For implementation details, see the [HNDP A* Subsolver Tutorial](../examples/hndp/astar.md).

## Related Pages

- [SubSolver Interface](../sub_solvers.md)
- [HNDP A* Subsolver Tutorial](../examples/hndp/astar.md)
