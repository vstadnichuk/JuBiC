using JuBiC, Documenter

makedocs(
    sitename = "JuBiC",
    remotes = nothing,
    clean = true,
    format = Documenter.HTML(
        collapselevel = 1,
        repolink = "https://github.com/vstadnichuk/JuBiC",
    ),
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Using JuBiC" => "model_objects.md",
        "Solver Methods" => [
            "Generalized Benders Cuts (GBC)" => joinpath("solvers", "gbc.md"),
            "Benders-like Cuts (BlC)" => joinpath("solvers", "blc.md"),
            "Benders-like Cuts with Lagrangian Coefficients (BlCLag)" => joinpath("solvers", "blc_lag.md"),
            "Compact MIP Wrapper" => joinpath("solvers", "mip.md"),
            "Direct MiBS Wrapper" => joinpath("solvers", "mibs.md"),
        ],
        "SubSolvers" => [
            "Interface" => "sub_solvers.md",
            "SubSolverJuMP" => joinpath("subsolvers", "subsolver_jump.md"),
            "SubSolverBlCJuMP" => joinpath("subsolvers", "subsolver_blc_jump.md"),
            "SubSolverMiBS" => joinpath("subsolvers", "subsolver_mibs.md"),
            "AStarSolver" => joinpath("subsolvers", "astar.md"),
        ],
        "Numerics and Status Codes" => "numerics_and_status.md",
        "Examples" => [
            "HNDP" => [
                "Motivation and Instances" => joinpath("examples", "hndp", "motivation.md"),
                "Solver Models" => joinpath("examples", "hndp", "solvers.md"),
                "A* Subsolver Tutorial" => joinpath("examples", "hndp", "astar.md"),
                "Benchmark Pipeline" => joinpath("examples", "hndp", "benchmarks.md"),
            ],
        ],
        "Core Solver API" => "solver_api.md",
    ]
)
