# Extending JuBiC

JuBiC is implemented as a generic tool for multi-follower generalized
interdiction games. In particular, the decomposition methods are designed so
that users can implement their own subsolvers when the follower problems have
application-specific structure.

This is an important design goal of the package. At the same time, we want to
be transparent about the current state of the code and documentation. JuBiC grew
quickly and currently contains a comparatively broad set of methods that follow
different algorithmic philosophies. As a result, not every solver route fits
neatly into one uniform code pattern, and some implementation details may still
feel unintuitive.

We are a small team, and documenting and maintaining all parts of the code at
the level we would like is currently beyond our available time. We are working
on improving this, but progress is gradual.

## Recommended Starting Point

For implementing a custom solver or subsolver, the most practical starting point
is currently the HNDP example:

- [HNDP Motivation](examples/hndp/motivation.md)
- [HNDP Instances](examples/hndp/instances.md)
- [HNDP Solver Models](examples/hndp/solvers.md)
- [HNDP A* Subsolver Tutorial](examples/hndp/astar.md)
- [HNDP Benchmark Pipeline](examples/hndp/benchmarks.md)

The HNDP code shows how an application-specific problem family is generated,
converted into JuBiC wrappers, solved with different JuBiC methods, and tested
through the benchmark pipeline. In many cases, the most efficient way to build a
new application is to start from this example and adapt the relevant parts to
the new setting.

The core files worth inspecting are:

- `examples/HNDP/hndp_instances.jl`
- `examples/HNDP/hndp_network_generation.jl`
- `examples/HNDP/hndp_model_generation.jl`
- `examples/HNDP/hndp_astar_wrapper.jl`
- `examples/HNDP/hndp_experiment_runner.jl`

## Developer Perspective

From a developer perspective, implementing new components directly inside the
JuBiC interface is currently most likely worth the effort in two situations.

First, you want to use `GBC` or `BlCLag` for your application. These solvers
contain several acceleration techniques, solver-interface details, and numerical
stability safeguards. Reimplementing those reliably can take substantial time,
even when using an AI coding agent.

Second, you want to benchmark a new approach against the solver methods already
available in JuBiC. In that case, implementing your solver through the JuBiC
interface should be comparatively simple, and it allows you to reuse common
instance generation, parameter handling, and benchmark output.

This is not meant to discourage other applications. If you find another reason
to build on JuBiC, please do so and let us know. We would like to expand this
list as new use cases appear.

We also welcome applications or solvers that use JuBiC externally. If desired,
we can link such projects from this page.

## Contact

If you have concrete questions, bug reports, or suggestions for improving the
documentation or code structure, please contact us:

- `vladimir.stadnichuk@uni-kassel.de`

We welcome feedback, especially when it points to a concrete part of the code or
documentation that was hard to use.

## Disclaimer

JuBiC is still a research-code project developed partly in limited available
time. It may contain bugs. Although the repository contains a substantial unit
test suite, we strongly recommend verifying computational results whenever
possible.

If you detect a bug or a potential bug, please contact us. We will try to
investigate and fix it.
