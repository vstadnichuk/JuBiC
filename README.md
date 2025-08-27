# JuBiC

**JuBiC** is an open-source software package written in Julia for testing and developing different kinds of **MIP-based decomposition algorithms**.  

The algorithms implemented in JuBiC are based on the dissertation of **Vladimir Stadnichuk**, with the main ideas summarized in our preprint:  
> Stadnichuk and Koster (2024), *Solving Multi-Follower Mixed-Integer Bilevel Problems with Binary Linking Variables*, Optimization Online. [Link](https://optimization-online.org/?p=28877)

We would like to thank **Cedrick Wieber** who significantly contributed to the implementation.  

While JuBiC was originally developed for **bilevel optimization problems**, it can also be naturally applied to **two-stage optimization problems**.

---

## Key Features

- Flexible interface for implementing your own decomposition methods or using the already available ones.  
- Solvers are structured around a **master solver** (for the leader/master problem) and **subsolvers** (for each follower/subproblem).  
- Subproblems can be solved with different solvers, offering high flexibility.  
- Currently includes **four solvers**:

### Implemented Solvers

1. **GBC-Solver**  
   - Implements the hierarchical decomposition from our preprint.  
   - Provides a wrapper for an MIP-based subproblem solver.  
   - Alternatively, you may implement your own subsolver for the Pricing Problem.  

2. **BLC-Solver**  
   - A basic implementation of a Benders-like decomposition for bilevel optimization.  
   - Requires the user to provide a function that generates the *big-M* values for Benders-like cuts.  

3. **MIBS-Solver**  
   - Wrapper for the bilevel solver **MIBS**.  

4. **MIP-Solver**  
   - Wrapper for solving standard **MIPs**.  

---

## Getting Started

1. Check the [`test`](./test) folder for **basic examples** that demonstrate the core functionalities.  
   - Run `runtest.jl` to verify that the installation was successful.  

2. Explore the [`examples`](./examples) folder for advanced applications:  
   - **HNDP**: Two-stage and bilevel network design problems.  
     - Demonstrates GBC-Solver and BLC-Solver. Also shows how MIP-Solver can be applied to solve the compact reformulations.  
     - The function *testrun!* in the file *hndp_tests.jl* demonstrates how a large number of tests can be automated within the JuBiC framework.  
   - **StochasticMultipleKnapsack**: Application of our solvers to the **two-stage stochastic multiple knapsack problem** from SIPLIB.  

---

## Disclaimer

JuBiC is still under active development.  

- A detailed documentation of the solvers and further examples are in progress.  
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
