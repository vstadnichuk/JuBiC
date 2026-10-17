using JuBiC, Documenter

# The detailed HNDP guides have one maintained source beside the example code.
# Copy them into the site at build time and translate repository links for HTML.
const REPO_URL = "https://github.com/vstadnichuk/JuBiC"
const SOURCE_REF = get(ENV, "GITHUB_SHA", "v0.6")
const HNDP_SOURCE = normpath(joinpath(@__DIR__, "..", "examples", "HNDP"))
const HNDP_GENERATED = joinpath(@__DIR__, "src", "generated", "hndp")
mkpath(HNDP_GENERATED)
for (source, target) in [
    ("NETWORK_GENERATION.md", "generators.md"),
    ("HNDP_SOLVER_OVERVIEW.md", "formulations.md"),
]
    text = read(joinpath(HNDP_SOURCE, source), String)
    text = replace(text,
        "](NETWORK_GENERATION.md)" => "](generators.md)",
        "](HNDP_SOLVER_OVERVIEW.md)" => "](formulations.md)",
        "](hndp_network_generation.jl)" => "]($(REPO_URL)/blob/$(SOURCE_REF)/examples/HNDP/hndp_network_generation.jl)",
        "](hndp_model_generation.jl)" => "]($(REPO_URL)/blob/$(SOURCE_REF)/examples/HNDP/hndp_model_generation.jl)",
        "](../../src/solvers/solver_parameters.jl)" => "]($(REPO_URL)/blob/$(SOURCE_REF)/src/solvers/solver_parameters.jl)",
    )
    metadata = "```@meta\nEditURL = \"$(REPO_URL)/blob/$(SOURCE_REF)/examples/HNDP/$(source)\"\n```\n\n"
    write(joinpath(HNDP_GENERATED, target), metadata * text)
end
for name in ("multimodal_example.json", "network_generation_example.json")
    cp(joinpath(HNDP_SOURCE, name), joinpath(HNDP_GENERATED, name); force=true)
end

makedocs(
    sitename = "JuBiC",
    root = @__DIR__,
    build = get(ENV, "JUBIC_DOCS_BUILD", "build"),
    modules = [JuBiC],
    checkdocs = :none,
    doctest = false,
    warnonly = false,
    remotes = nothing,
    clean = true,
    format = Documenter.HTML(
        collapselevel = 1,
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://vstadnichuk.github.io/JuBiC/",
        edit_link = nothing,
        repolink = "https://github.com/vstadnichuk/JuBiC",
    ),
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Overview" => "model_objects.md",
        "Implemented Solvers" => [
            "Generalized Benders Cuts (GBC)" => joinpath("solvers", "gbc.md"),
            "Benders-like Cuts (BlC)" => joinpath("solvers", "blc.md"),
            "Benders-like Cuts with Lagrangian Coefficients (BlCLag)" => joinpath("solvers", "blc_lag.md"),
            "Compact MIP Wrapper" => joinpath("solvers", "mip.md"),
            "Direct MiBS Wrapper" => joinpath("solvers", "mibs.md"),
        ],
        "Implemented Follower Subsolvers" => [
            "Abstract Subsolver Interface" => "sub_solvers.md",
            "SubSolverJuMP" => joinpath("subsolvers", "subsolver_jump.md"),
            "SubSolverBlCJuMP" => joinpath("subsolvers", "subsolver_blc_jump.md"),
            "SubSolverMiBS" => joinpath("subsolvers", "subsolver_mibs.md"),
            "AStarSolver" => joinpath("subsolvers", "astar.md"),
        ],
        "Custom (Sub)Solver" => "extending_jubic.md",
        "Bilevel Network Design" => [
            "Motivation" => joinpath("examples", "hndp", "motivation.md"),
            "Instances" => joinpath("examples", "hndp", "instances.md"),
            "Generator Reference and Examples" => "generated/hndp/generators.md",
            "Solver Models" => joinpath("examples", "hndp", "solvers.md"),
            "Formulation and Parameter Reference" => "generated/hndp/formulations.md",
            "A* Subsolver Tutorial" => joinpath("examples", "hndp", "astar.md"),
            "Benchmark Pipeline" => joinpath("examples", "hndp", "benchmarks.md"),
        ],
        "Numerics and Status Codes" => "numerics_and_status.md",
        "Examples" => [
            "Small Bilevel Model" => "examples/example1.md",
            "Reading MPS/AUX Input" => "examples/example2.md",
            "Routing Follower Oracle" => "examples/example3.md",
        ],
        "Package Structure" => "package_structure.md",
        "Core Solver API" => "solver_api.md",
        "Building and Publishing the Docs" => "documentation.md",
    ]
)
