"""Download the OSM and GTFS sources of the six study cities.

    python scripts/download_sources.py [city ...]

For every city of settings.json, the files listed in examples/HNDP/sources/osm_city_sources.json are
downloaded to <study>/source_data/sources/<city>/:
  roads.osm.pbf         BBBike road extract (covers the municipal boundary plus a 5 km buffer)
  gtfs.zip              GTFS feed of the local transit authority (Dresden: also regional_gtfs.zip)
  boundary_full.json    municipal boundary relation from the OpenStreetMap API
  boundary.json         the boundary relation with way geometries (derived from boundary_full.json)

The registry records the SHA-256 checksums of the snapshot used in the paper (OSM data of
24 October 2026, GTFS feeds downloaded on 28 October 2026). The providers only serve their
current data, so new downloads usually differ from that snapshot; the script reports this for
every file but keeps the download. Existing files are not downloaded again.
"""
import hashlib, json, sys, urllib.request
from paths import B, CITIES, REPO

REGISTRY = json.loads((REPO / "examples/HNDP/sources/osm_city_sources.json").read_text(encoding="utf-8"))
SOURCES = B / "source_data" / "sources"
HEADERS = {"User-Agent": "JuBiC-station-location-reproduction"}


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def download(url, path):
    if path.exists():
        print("cached     ", path.relative_to(B), flush=True)
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".part")
    print("downloading", url, flush=True)
    with urllib.request.urlopen(urllib.request.Request(url, headers=HEADERS), timeout=600) as response, \
            partial.open("wb") as f:
        while block := response.read(1 << 20):
            f.write(block)
    partial.replace(path)


def boundary_geometry(record, folder):
    """Attach the way geometries to the boundary relation (input format of the generator)."""
    doc = json.loads((folder / "boundary_full.json").read_text(encoding="utf-8"))
    nodes = {x["id"]: x for x in doc["elements"] if x["type"] == "node"}
    ways = {x["id"]: x for x in doc["elements"] if x["type"] == "way"}
    relation = next(x for x in doc["elements"] if x["type"] == "relation" and x["id"] == record["relation_id"])
    for member in relation["members"]:
        if member["type"] == "way":
            member["geometry"] = [{"lat": nodes[n]["lat"], "lon": nodes[n]["lon"]} for n in ways[member["ref"]]["nodes"]]
    (folder / "boundary.json").write_text(json.dumps({"elements": [relation]}), encoding="utf-8")


def main(cities):
    records = {r["city"]: r for r in REGISTRY["cities"]}
    report = {}
    for city in cities:
        record = records[city]; folder = SOURCES / city
        files = [(f["file"].split("/")[-1], f["url"], f.get("sha256")) for f in record["downloaded_files"]]
        files.append(("boundary_full.json", record["boundary_url"], None))
        for name, url, expected in files:
            path = folder / name
            download(url, path)
            digest = sha256(path)
            same = None if expected is None else digest == expected
            report[f"{city}/{name}"] = dict(url=url, sha256=digest, same_as_paper_snapshot=same)
            if same is False:
                print(f"  note: {city}/{name} differs from the paper's snapshot (provider data changed)", flush=True)
        boundary_geometry(record, folder)
    (B / "source_data" / "download_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main(sys.argv[1:] or CITIES)
