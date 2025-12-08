# JuBiC Documentation

## Introduction

**JuBiC** is an open-source software package developed in Julia specifically designed for solving bilevel optimization problems and two-stage optimization problems that involve binary linking variables. These problems feature a hierarchical structure consisting of a first-level (leader) problem and one or more second-level (follower) problems. Such problems arise in various fields, including economics, engineering, and logistics, where decisions made at the first level influence the feasible region of the second-level problems.


### Problem Formulation

Let $\mathcal{A}$ denote a set of shared resources. For each resource $a \in \mathcal{A}$, the first-level problem uses binary linking variables $x_a \in \{0, 1\}$, and we denote the vector of all $x$-variables as $\mathbf{x}$. Let $K$ be the total number of second-level problems. For each second-level problem $k \in \{1, \dots, K\}$ and resource $a \in \mathcal{A}$, let $C_a^k > 0$ denote the resource capacity available to the problem, and let $y_a^k \geq 0$ represent the variable in the second-level problem that is linked to the corresponding first-level variable $x_a$. Similarly, $\mathbf{y}^k$ denotes the vector of all $y$-variables of the $k$-th second-level problem. Furthermore, let $\mathbf{h}$ and $\mathbf{z}^k$ denote the (bounded) feasible regions for purely first-level and purely second-level variables. The functions $r(\cdot), r_1(\cdot), \dots, r_K(\cdot)$ and $c_1(\cdot), \dots, c_K(\cdot)$ represent the (linear) objective functions for the first-level and all second-level problems, respectively.

The first-level problem is defined as:

```math
\begin{align*}
    \min \quad & r(\mathbf{h}, \mathbf{x}) + \sum_{k=1}^K r_k(\mathbf{y}^k, \mathbf{z}^k) & \\
    s.t. \quad & \mathbf{h} \in \mathcal{H}(\mathbf{x}) & \\
               & x_a \in \{0, 1\}   \quad &\forall a \in \mathcal{A} \\
               & (\mathbf{y}^k, \mathbf{z}^k) \in \mathcal{S}_k(\mathbf{x}) &\forall k \in \{1, \dots, K\}
\end{align*}
```

Here, $\mathcal{S}_k(\mathbf{x})$ represents the set of optimal solutions to the $k$-th second-level problem defined as: 

```math
\begin{align*}
    \min \quad & c_k(\mathbf{y}^k, \mathbf{z}^k) & \\
    s.t. \quad & y_a^k \leq C_a^k x_a \quad &\forall a \in \mathcal{A} \\
               & y_a^k \geq 0   \quad &\forall a \in \mathcal{A} \\
               & \mathbf{z}^k \in \mathcal{Z}_k(\mathbf{y}^k) &
\end{align*}
```

## Installation

You can install JuBiC via Julia's package manager as follows:

```julia-repl
julia> import Pkg

julia> Pkg.add("JuBiC")
```

Alternatively, the latest development version can be installed directly from GitHub:

```julia-repl
julia> import Pkg

julia> Pkg.add(url="https://github.com/vstadnichuk/JuBiC")
```

## Verify Installation

To ensure that everything is set up correctly after the installation, you can run the following command to execute the test provided within JuBiC:

```julia-repl
julia> Pkg.test("JuBiC")
```

## License

JuBiC is released as an **open-source project**. For detailed licensing information, please refer to the [LICENSE](https://github.com/vstadnichuk/JuBiC/blob/main/LICENSE).

## Citing `JuBiC`

If you use **JuBiC** in your research, please cite:

> Stadnichuk and Koster (2024), *Solving Multi-Follower Mixed-Integer Bilevel Problems with Binary Linking Variables*, Optimization Online. [Link](https://optimization-online.org/?p=28877)