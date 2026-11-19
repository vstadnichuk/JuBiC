using JuBiC, Documenter

# The detailed HNDP guides have one maintained source beside the example code.
# Copy them into the site at build time and translate repository links for HTML.
const REPO_URL = "https://github.com/vstadnichuk/JuBiC"
const SOURCE_REF = get(ENV, "GITHUB_SHA", "bilevel-feasible-paths-for-the-station-location-problem")
const HNDP_SOURCE = normpath(joinpath(@__DIR__, "..", "examples", "HNDP"))
const HNDP_GENERATED = joinpath(@__DIR__, "src", "generated", "hndp")
rm(HNDP_GENERATED; force=true, recursive=true)  # drop copies of renamed or removed guides
mkpath(HNDP_GENERATED)
for (source, target) in [
    ("README.md", "overview.md"),
    ("OSM_NETWORKS.md", "osm_networks.md"),
    ("DEMAND_AND_AGGREGATION.md", "demand_and_aggregation.md"),
    ("STATION_PATH_MODEL.md", "station_path_model.md"),
]
    text = read(joinpath(HNDP_SOURCE, source), String)
    text = replace(text,
        # Guides link to each other, also to a heading (`](OSM_NETWORKS.md#sources)`).
        "](README.md" => "](overview.md",
        "](OSM_NETWORKS.md" => "](osm_networks.md",
        "](DEMAND_AND_AGGREGATION.md" => "](demand_and_aggregation.md",
        "](STATION_PATH_MODEL.md" => "](station_path_model.md",
        # Code and data files live in subfolders of examples/HNDP; link them on GitHub.
        "](models/" => "]($(REPO_URL)/blob/$(SOURCE_REF)/examples/HNDP/models/",
        "](network_generation/" => "]($(REPO_URL)/blob/$(SOURCE_REF)/examples/HNDP/network_generation/",
        "](demand/" => "]($(REPO_URL)/blob/$(SOURCE_REF)/examples/HNDP/demand/",
        "](experiments/" => "]($(REPO_URL)/blob/$(SOURCE_REF)/examples/HNDP/experiments/",
        "](sources/" => "]($(REPO_URL)/blob/$(SOURCE_REF)/examples/HNDP/sources/",
        "](tests/" => "]($(REPO_URL)/blob/$(SOURCE_REF)/examples/HNDP/tests/",
        "](../../src/solvers/solver_parameters.jl)" => "]($(REPO_URL)/blob/$(SOURCE_REF)/src/solvers/solver_parameters.jl)",
        "](station_location_paper/" => "]($(REPO_URL)/blob/$(SOURCE_REF)/examples/HNDP/station_location_paper/",
        "](../../README.md)" => "]($(REPO_URL)/blob/$(SOURCE_REF)/README.md)",
    )
    metadata = "```@meta\nEditURL = \"$(REPO_URL)/blob/$(SOURCE_REF)/examples/HNDP/$(source)\"\n```\n\n"
    write(joinpath(HNDP_GENERATED, target), metadata * text)
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
            "Station Location Paper" => "generated/hndp/overview.md",
            "City Networks from OSM and GTFS" => "generated/hndp/osm_networks.md",
            "Demand and Aggregation" => "generated/hndp/demand_and_aggregation.md",
            "Station-Based Path Model" => "generated/hndp/station_path_model.md",
            "Solver Models" => joinpath("examples", "hndp", "solvers.md"),
            "Formulation and Parameter Reference" => joinpath("examples", "hndp", "formulations.md"),
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
