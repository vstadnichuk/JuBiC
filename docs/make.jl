using JuBiC, Documenter

makedocs(
    sitename = "JuBiC",
    remotes = nothing,
    clean = true,
    format = Documenter.HTML(
        collapselevel = 1,
    ),
    pages = [
        "Home" => "index.md",
        "Examples" => joinpath.("examples", ["example1.md", "example2.md", "example3.md"]),
        "Solver Methods" => joinpath.("solvers", ["gbc.md", "blc.md", "blc_lag.md", "mip.md", "mibs.md"]),
        "Implementation of SubSolvers" => "sub_solvers.md"
    ]
)
