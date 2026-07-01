# JuBiC

**JuBiC** is an open-source software package written in Julia for testing and developing different kinds of **MIP-based decomposition algorithms**.  

The algorithms implemented in JuBiC are based on the dissertation of **Vladimir Stadnichuk**, with the main ideas summarized in our preprint:  
> Stadnichuk and Koster (2024), *Solving Multi-Follower Mixed-Integer Bilevel Problems with Binary Linking Variables*, Optimization Online. [Link](https://optimization-online.org/?p=28877)

We would like to thank **Cedrick Wieber** who significantly contributed to the implementation.  

While JuBiC was originally developed for **bilevel optimization problems**, it can also be naturally applied to **two-stage optimization problems**.

---

## Documentation

The current documentation starts at [`docs/src/index.md`](./docs/src/index.md).

To build the rendered HTML documentation locally, run:

```julia
julia --project=. docs/make.jl
```


---

## Paper-Specific Branches

Some papers have been written with JuBiC. For reproducibility, these papers may
have dedicated branches that preserve the state of the code used when the paper
was published.

Currently, the available paper-specific branch is:

- [`catalog-formulations-multi-follower-discrete-bilevel-network-design`](https://github.com/vstadnichuk/JuBiC/tree/catalog-formulations-multi-follower-discrete-bilevel-network-design):
  branch for
  [*A Catalog of Formulations for the Multi-Follower Discrete Bilevel Network Design Problem*](https://optimization-online.org/?p=35437).

---

### Implemented Solvers

1. **GBC-Solver**
   - Implements generalized Benders cuts in the first-level variable space.
   - Supports JuMP, MiBS, and custom subsolvers.

2. **BlC-Solver**
   - Implements Benders-like cuts for bilevel optimization.
   - Requires user-provided big-M values for the Benders-like cuts.

3. **BlCLag-Solver**
   - Uses the BlC structure but generates cut coefficients from the problem structure.
   - Requires bilevel-capable subsolvers.

4. **MiBS-Solver**
   - Direct wrapper for the external bilevel solver **MiBS** through BilevelJuMP.

5. **MIP-Solver**
   - Wrapper for solving compact JuMP MIP models through JuBiC's common interface.

---

## Getting Started

1. Check the [`test`](./test) folder for **basic examples** that demonstrate the core functionalities.  
   - Run `runtest.jl` to verify that the installation was successful.  

2. Explore the [`examples`](./examples) folder for advanced applications:  
   - **HNDP**: Bilevel network design problems.  
     - Demonstrates GBC, BlC, BlCLag, MiBS, and compact MIP-based formulations.  
     - The function *testrun!* in the file *hndp_tests.jl* demonstrates how a large number of tests can be automated within the JuBiC framework.  

---

## Disclaimer

JuBiC is still under active development.  

- Further examples are in progress.  
- Results should be interpreted with care. While extensive experiments have been conducted to validate the solvers, we **strongly recommend cross-verifying results**.  

If you encounter difficulties or have any questions, feel free to contact:  
📧 **Vladimir.stadnichuk@om.rwth-aachen.de**  

---

## Citation

If you use **JuBiC** in your research, please cite:  

> Stadnichuk and Koster (2024), *Solving Multi-Follower Mixed-Integer Bilevel Problems with Binary Linking Variables*, Optimization Online. [Link](https://optimization-online.org/?p=28877)  

---

## License

JuBiC is released as an **open-source project**. See the [LICENSE](./LICENSE) file for details.  
