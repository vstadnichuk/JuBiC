"""Convert the downloaded sources into the input format of the network generator.

    python scripts/prepare_sources.py [city ...]          (requires pyosmium, see requirements.txt)

For every city, reads <study>/source_data/sources/<city>/{roads.osm.pbf,gtfs.zip[,regional_gtfs.zip]}
and writes to the same folder
  roads_buffered.json   all OSM highways (and traffic signals) intersecting the city rectangle
                        (municipal bounding box plus 5 km, see osm_city_sources.json), in
                        Overpass JSON format with way geometries
  gtfs/stops.txt, gtfs/stop_times.txt
                        all trips that serve at least one stop inside the rectangle, including
                        their stops outside it; IDs are prefixed with the feed index
"""
import csv, io, json, pathlib, sys, zipfile
import osmium
from paths import B, CITIES, REPO

REGISTRY = {r["city"]: r for r in json.loads((REPO / "examples/HNDP/sources/osm_city_sources.json").read_text(encoding="utf-8"))["cities"]}
SOURCES = B / "source_data" / "sources"


class Roads(osmium.SimpleHandler):
    def __init__(self, bbox):
        super().__init__(); self.bbox = bbox; self.elements = []

    def inside(self, lat, lon):
        s, w, n, e = self.bbox
        return s <= lat <= n and w <= lon <= e

    def node(self, node):
        if node.tags.get("highway") == "traffic_signals" and self.inside(node.location.lat, node.location.lon):
            self.elements.append({"type": "node", "id": node.id, "lat": node.location.lat, "lon": node.location.lon,
                                  "tags": {"highway": "traffic_signals"}})

    def way(self, way):
        if "highway" not in way.tags:
            return
        geom = [{"lat": p.lat, "lon": p.lon} for p in way.nodes]
        if len(geom) < 2 or not any(self.inside(p["lat"], p["lon"]) for p in geom):
            return
        self.elements.append({"type": "way", "id": way.id, "tags": dict(way.tags), "geometry": geom})


def gtfs_subset(city, path, bbox):
    target = path / "gtfs"; target.mkdir(exist_ok=True)
    # Keep every row of trips touching the study rectangle, including outside stops,
    # to preserve boundary-spanning services and their observed timing.
    archives = [path / "gtfs.zip"] + ([path / "regional_gtfs.zip"] if (path / "regional_gtfs.zip").exists() else [])
    stops, times = [], []
    for prefix, archive in enumerate(archives):
        with zipfile.ZipFile(archive) as z:
            names = {pathlib.PurePosixPath(n).name: n for n in z.namelist()}
            def rows(name):
                return csv.DictReader(io.TextIOWrapper(z.open(names[name]), encoding="utf-8-sig", newline=""))
            selected = set(); allstops = {}
            for row in rows("stops.txt"):
                sid = row["stop_id"]; allstops[sid] = row
                try:
                    lat, lon = float(row["stop_lat"]), float(row["stop_lon"])
                    if bbox[0] <= lat <= bbox[2] and bbox[1] <= lon <= bbox[3]:
                        selected.add(sid)
                except (ValueError, KeyError):
                    pass
            trips = {row["trip_id"] for row in rows("stop_times.txt") if row["stop_id"] in selected}
            used = set(); count = 0
            for row in rows("stop_times.txt"):
                if row["trip_id"] not in trips:
                    continue
                used.add(row["stop_id"]); row["trip_id"] = f"{prefix}:" + row["trip_id"]
                row["stop_id"] = f"{prefix}:" + row["stop_id"]; times.append(row); count += 1
            for sid in sorted(used):   # deterministic row order
                if sid in allstops:
                    row = allstops[sid]; row["stop_id"] = f"{prefix}:" + sid; stops.append(row)
            print(city, "GTFS", archive.name, "trips", len(trips), "rows", count, flush=True)
    for name, data, fields in [("stops.txt", stops, ["stop_id", "stop_name", "stop_lat", "stop_lon"]),
                               ("stop_times.txt", times, ["trip_id", "arrival_time", "departure_time", "stop_id", "stop_sequence"])]:
        with (target / name).open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fields, extrasaction="ignore"); writer.writeheader(); writer.writerows(data)


def main(cities):
    for city in cities:
        record = REGISTRY[city]; path = SOURCES / city; dest = path / "roads_buffered.json"
        if not record["bbbike_covers_buffer"]:
            raise RuntimeError("The road extract does not cover the study rectangle: " + city)
        if not dest.exists():
            print("Preparing roads", city, flush=True)
            handler = Roads(record["bbox"])
            handler.apply_file(str(path / "roads.osm.pbf"), locations=True, idx="flex_mem")
            dest.write_text(json.dumps({"generator": "JuBiC PBF-to-Overpass geometry adapter", "buffer_km": 5,
                                        "bbox": record["bbox"], "elements": handler.elements},
                                       separators=(",", ":")), encoding="utf-8")
            print(city, "OSM elements", len(handler.elements), flush=True)
        if not (path / "gtfs/stop_times.txt").exists():
            gtfs_subset(city, path, record["bbox"])


if __name__ == "__main__":
    main(sys.argv[1:] or CITIES)
