# JuBiC

**JuBiC** is an open-source software package written in Julia for solving mixed-integer bilevel optimization problems with binary linking variables, with a strong focus on decomposition methods and custom approaches for bilevel network design problems.

The theoretical foundations of these algorithms are described in [Vladimir Stadnichuk's dissertation](https://publications.rwth-aachen.de/record/1039116/).

We would like to thank **Cedrick Wieber**, who contributed significantly to the implementation.

---

## Getting Started

1. Read the [HTML documentation](https://vstadnichuk.github.io/JuBiC/).

2. Check the [`test`](./test) folder for **basic examples** that demonstrate the core functionality.
   - Run `julia --project=. test/runtests.jl` from the repository root to verify the installation (requires Gurobi and MiBS).

3. Explore the [`examples`](./examples) folder for advanced applications:
   - **HNDP**: Bilevel network design problems.
     - Demonstrates GBC, BlC, BlCLag, MiBS, and compact MIP-based formulations.
     - Start with the [HNDP guide](./examples/HNDP/README.md); `run_hndp_manifest.jl` runs a configuration file or a three-file experiment directory.

---

## Implemented Solvers

1. **GBC solver**
   - Implements generalized Benders cuts in the first-level variable space.
   - Supports JuMP, MiBS, and custom subsolvers.

2. **BlC solver**
   - Implements Benders-like cuts for bilevel optimization.
   - Requires user-provided big-M values for the Benders-like cuts.

3. **BlCLag solver**
   - Uses the BlC structure but generates cut coefficients from the problem structure.
   - Requires bilevel-capable subsolvers.

4. **MiBS solver**
   - Direct wrapper for the external bilevel solver **MiBS** through BilevelJuMP.

5. **MIP solver**
   - Wrapper for solving compact JuMP MIP models through JuBiC's common interface.

---

## Paper-Specific Branches

The following branches contain code associated with three papers:

- [*Solving Multi-Follower Mixed-Integer Bilevel Problems with Binary Linking Variables*](https://optimization-online.org/2024/12/solving-multi-follower-mixed-integer-bilevel-problems-with-binary-linking-variables/) (2024):
  [`solving-multi-follower-mixed-integer-bilevel-problems-with-binary-linking-variables`](https://github.com/vstadnichuk/JuBiC/tree/solving-multi-follower-mixed-integer-bilevel-problems-with-binary-linking-variables).
- [*Automated Benders-like Cut Generation and its Application to the Bilevel Network Design Problem*](https://www.euro-online.org/inoc2026/conference-program/schedule/extended_abstracts/stadnichuk.pdf) (INOC 2026 extended abstract):
  [`automated-benders-like-cut-generation-and-its-application-to-the-bilevel-network-design-problem`](https://github.com/vstadnichuk/JuBiC/tree/automated-benders-like-cut-generation-and-its-application-to-the-bilevel-network-design-problem).
- [*A Catalog of Formulations for the Multi-Follower Discrete Bilevel Network Design Problem*](https://optimization-online.org/2026/06/a-catalog-of-formulations-for-the-multi-follower-discrete-bilevel-network-design-problem/) (2026):
  [`catalog-formulations-multi-follower-discrete-bilevel-network-design`](https://github.com/vstadnichuk/JuBiC/tree/catalog-formulations-multi-follower-discrete-bilevel-network-design).

---

## Disclaimer

JuBiC is still under active development. Therefore, results should be interpreted with care. While extensive experiments have been conducted to validate the solvers, we **strongly recommend cross-verifying results**.

If you encounter difficulties or have any questions, feel free to contact:
📧 **vladimir.stadnichuk@uni-kassel.de**

---

## Citation

If you use **JuBiC** in your research, please cite:

> Stadnichuk and Koster (2024), *Solving Multi-Follower Mixed-Integer Bilevel Problems with Binary Linking Variables*, Optimization Online. [Link](https://optimization-online.org/?p=28877)

---

## License

JuBiC is released as an **open-source project**. See the [LICENSE](./LICENSE) file for details.
