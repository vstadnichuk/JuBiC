"""Shared locations for the station location paper scripts (see paths.jl).

CODE  examples/HNDP/station_location_paper
REPO  JuBiC repository root (override with JUBIC_REPO)
B     study folder with inputs, raw results and derived data (override with SLP_STUDY_DIR)
JULIA Julia executable (override with JULIA; default: julia on PATH)
"""
from pathlib import Path
import json, os, shutil

CODE = Path(__file__).resolve().parents[1]
REPO = Path(os.environ.get("JUBIC_REPO", CODE.parents[2])).resolve()
B = Path(os.environ.get("SLP_STUDY_DIR", CODE / "study")).resolve()
(B / "data").mkdir(parents=True, exist_ok=True)
JULIA = os.environ.get("JULIA") or shutil.which("julia") or "julia"
SETTINGS = json.loads((CODE / "settings.json").read_text())
CITIES = SETTINGS["cities"]
FIGURES = B / "figures"
