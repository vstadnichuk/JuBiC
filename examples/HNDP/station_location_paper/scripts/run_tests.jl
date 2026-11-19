# Unit tests of the station-based path model and the enumeration of bilevel feasible paths.
#   julia --project=../../.. --threads=8 scripts/run_tests.jl
using Test
include(joinpath(@__DIR__, "paths.jl"))
include(joinpath(SRC, "tests", "osm_station_path_tests.jl"))
@testset "station path model" begin
    hndp_station_path_tests()
end
