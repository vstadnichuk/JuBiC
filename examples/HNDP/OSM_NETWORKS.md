# City networks from OSM and GTFS

This guide describes how the six full city networks of the paper are built
(Section 5 and Appendix F of the paper).

## Sources

| City | BBBike road extract | GTFS supplier | Car availability |
|---|---|---|---:|
| Aachen | Aachen | AVV | 64% |
| Cologne | Koeln | VRS | 71% |
| Heidelberg | Mannheim | VRN | 63.2% |
| Karlsruhe | Karlsruhe | KVV | 71% |
| Nuremberg | Nuernberg | VGN | 70% |
| Dresden | Dresden | GTFS.de local/regional | 63% |

The registry [osm_city_sources.json](sources/osm_city_sources.json) lists, per city,
the download URLs, the municipal boundary relation, the study rectangle (municipal
bounding box plus 5 km) and the checksums of the paper's snapshot (OSM of
24 October 2026, GTFS feeds downloaded on 28 October 2026).

- `download_sources.py` downloads the extracts, feeds and boundaries into
  `study/source_data/sources/<city>/` and reports which files differ from the snapshot.
- `prepare_sources.py` reads the PBF extracts with pyosmium and writes the road
  ways, the boundary geometry and the GTFS trips within the study rectangle.

## Network construction

`generate_networks.jl` calls `generate_osm_gtfs`
([hndp_osm_gtfs_generation.jl](network_generation/hndp_osm_gtfs_generation.jl))
with its default parameters for every city:

1. **Roads.** Each OSM way becomes one arc between its two end points, with the
   length of its geometry; intermediate geometry points are not nodes. Ways with
   the same end points give one arc per direction.
2. **Boundary.** Walking and bike arcs are kept only inside the municipal boundary.
   Car routes outside it are replaced by shortcut arcs from boundary nodes back into
   the city, so cars may still leave the city on their way. Only the largest
   connected car component is kept.
3. **Layers.** Each retained physical location has one node in each of the four
   layers walking, car, transit and bike. Road geometry gives the walking, car and
   bike arcs; GTFS stop sequences give the transit arcs. Transfer arcs connect the
   layers at the same location.
4. **Stations.** Candidate sharing stations are existing road nodes matched to a
   GTFS stop within 250 m. Each station is a decision group of one bike entry arc
   and one bike return arc that are opened together.
5. **Walking cleanup.** Only the largest strongly connected component of the
   walking layer is kept, with all four copies of each removed location
   ([hndp_osm_connectivity.jl](network_generation/hndp_osm_connectivity.jl)). Walking
   is therefore a fallback route for every OD pair.

| City | Locations before cleanup | Retained | Removed | Candidate stations |
|---|---:|---:|---:|---:|
| Aachen | 5,572 | 5,074 | 498 | 709 |
| Cologne | 18,992 | 16,801 | 2,191 | 820 |
| Heidelberg | 4,297 | 4,006 | 291 | 420 |
| Karlsruhe | 10,248 | 9,300 | 948 | 744 |
| Nuremberg | 8,202 | 7,661 | 541 | 1,021 |
| Dresden | 13,944 | 13,427 | 517 | 1,296 |

## Arc costs

All costs are integer euro cents. Time is valued at EUR 13 per person-hour
(Richter 2022), so one minute costs 21.67 cents.

| Arc type | Cost |
|---|---|
| Walking, bike or transit movement | Travel time, rounded to cents |
| Car movement | Rounded time cost plus 22 cents/km, rounded again (ties to even) |
| Car entry / exit | 11 / 103 cents (30 s / 90 s; the exit includes 71 cents parking) |
| Transit entry / exit | 105 cents plus waiting time / zero |
| Bike entry / return | 83 / 11 cents (60 s plus the 61-cent fee / 30 s) |

- **Walking:** 5 km/h.
- **Car:** speed is 0.5 times the OSM `maxspeed`, clamped to 10–110 km/h and
  rounded. Ways with a missing or unsupported value use 20 km/h. The 22 cents/km
  are the cambio COMFORT base distance rate. The 71-cent parking charge is derived
  from the resident parking permit formula.
- **Bike:** 20 km/h on cycleways; 18 on tracks, living streets and tertiary roads;
  17 on residential, primary and secondary roads; 16 on unclassified roads; 15 on
  paths and service roads; 12 on footways and pedestrian zones. Ways tagged for
  cycling get at least 19 km/h. The operator earns the 61-cent fee per entry.
- **Transit:** movement uses the median GTFS stop-to-stop time. The entry fee is
  EUR 63 per month over 60 rides. Waiting is half the pooled headway over a 16-hour
  service day, capped at 10 minutes. 

The car costs are written separately as `speed_overrides/<city>.csv` and are applied
on top of the saved network by every later step.

## Physical arc lengths

Costs alone do not give trip lengths, which the demand model needs for its distance
bands. The lengths are recovered separately:

- **Roads:** `recover_osm_road_lengths`
  ([osm_road_arc_lengths.jl](network_generation/osm_road_arc_lengths.jl)), called by
  `generate_networks.jl`, recomputes the length of every road arc, including the
  car shortcuts, from the original OSM way geometry.
- **Transit:** [osm_transit_arc_lengths.py](network_generation/osm_transit_arc_lengths.py)
  measures each transit arc along the GTFS shapes. It snaps stops to the node set
  before the walking cleanup, so `generate_networks.jl` also stores those
  coordinates. Feeds without shapes (Karlsruhe, Nuremberg, Dresden) use the
  straight distance between consecutive stops, tagged `gtfs_stop_chord_fallback`.
- **Transfers** between layers at the same location have length zero.

## Outputs

`study/inputs/city_networks_walk_connected/`:

- `graphs/<city>/edges.csv`, `metadata.json` and `physical_node_mapping.csv`;
- `speed_overrides/<city>.csv`, the car-layer arc costs.

`study/inputs/geometry/road_lengths/<city>.csv` and
`study/inputs/geometry/transit_lengths/<city>.csv` hold the physical lengths.

From the snapshot files, these outputs are identical to the paper's inputs. A new
download reflects the current OSM and timetable data instead.
