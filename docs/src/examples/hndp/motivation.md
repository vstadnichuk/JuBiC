# HNDP Motivation

Besides the basic implementations of the core algorithms, JuBiC provides
concrete examples showing how these algorithms can be applied to application
models. The example studied here is motivated by network design settings where
design and routing decisions are made by different actors.

The first-level actor is an operator who decides which arcs of a network are
available. The second-level actors are users who react to this design by
choosing routes. This naturally leads to a bilevel model: the operator controls
network availability, while each user solves a routing problem on the resulting
network.

## Mathematical Model

Let ``G = (V, E)`` be a directed graph and let ``A \subseteq E`` be the set of
operator-controlled decision arcs. The operator chooses

```math
x_a \in \{0,1\} \qquad a \in A,
```

where ``x_a = 1`` means that arc ``a`` is available. Arcs in ``E \setminus A`` are
fixed and are always available.

For each user ``u``, let ``o_u`` be the origin, ``d_u`` the destination,
``c^u_a`` the user cost of arc ``a``, ``r^u_a`` the operator-side contribution
induced by arc ``a``, and optionally ``w^u_a`` an arc weight. If the instance is
length constrained, user ``u`` also receives a path-weight bound ``W_u``.

In the arc-flow view, user ``u`` solves

```math
\begin{aligned}
\min_{y^u}\quad & \sum_{a \in E} c^u_a y^u_a \\
\text{s.t.}\quad
& y^u \text{ is an } o_u\text{-}d_u \text{ path flow}, \\
& y^u_a \le x_a \qquad a \in A, \\
& \sum_{a \in E} w^u_a y^u_a \le W_u \qquad \text{if a weight bound is used.}
\end{aligned}
```

The first-level objective contains construction costs and the user-side
contribution:

```math
\min_x \sum_{a \in A} p_a x_a
      + \sum_u \sum_{a \in E} r^u_a y^u_a .
```

## HNDP Documentation Outline

The remaining HNDP pages describe how this abstract model is used in the example
pipeline:

- [Instances](instances.md) explains how concrete graph topologies, costs,
  risks, weights, budgets, and user data are generated.
- [Solver Models](solvers.md) describes the HNDP formulations implemented on
  top of JuBiC solvers.
- [A* Subsolver Tutorial](astar.md) explains the custom HNDP labeling subsolver
  and how similar subsolvers can be implemented.
- [Benchmark Pipeline](benchmarks.md) explains how JSON configurations are used
  to run reproducible HNDP experiment sweeps.
