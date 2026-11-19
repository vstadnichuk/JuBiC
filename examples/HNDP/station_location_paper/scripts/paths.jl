# Shared locations for the station location paper scripts.
#   CODE  this folder's parent (examples/HNDP/station_location_paper)
#   REPO  JuBiC repository root (override with JUBIC_REPO)
#   B     study folder holding inputs, raw results and derived data (override with SLP_STUDY_DIR)
#   SRC   JuBiC HNDP example sources used by the study scripts
const CODE = normpath(joinpath(@__DIR__, ".."))
const REPO = normpath(get(ENV, "JUBIC_REPO", joinpath(CODE, "..", "..", "..")))
const B = normpath(get(ENV, "SLP_STUDY_DIR", joinpath(CODE, "study")))
const SRC = joinpath(REPO, "examples", "HNDP")
const SETTINGS_FILE = joinpath(CODE, "settings.json")
const PROVENANCE = joinpath(SRC, "sources")   # osm_distance_demand_sources.json, osm_car_ownership_sources.json
mkpath(joinpath(B, "data"))
