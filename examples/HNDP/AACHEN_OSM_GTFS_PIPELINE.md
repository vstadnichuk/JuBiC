# Reproducible OSM/GTFS city-instance pipeline

This note describes the current procedure for constructing a bounded,
multimodal JuBiC instance. It is intended as an internal recipe that can be
adapted to another city.

## 1. Obtain the source data

Create an ignored cache directory, for example `data/osm/`. Download:

1. An OSM road extract covering the study area **plus a geographic buffer**.
   The Aachen extract uses the Overpass query
   `[out:json][timeout:240];way["highway"](50.68,5.90,50.88,6.30);out geom;`.
2. The OSM administrative-boundary relation for the study area. Aachen uses
   relation `62564`, saved as `aachen_boundary.json`.
3. The city's GTFS feed. Aachen uses the official AVV feed
   `https://opendata.avv.de/current_GTFS/AVV_GTFS_Masten_mit_SPNV.zip`,
   extracted to `data/osm/avv_gtfs/`.

Raw downloads are not committed. Keep the exact URLs, query, date, and feed
version alongside any published experiment results.

## 2. Build the bounded road network

`aachen_osm_gtfs_generation.jl` reads the Overpass JSON and:

- retains usable motor-vehicle roads for the car layer;
- retains car-compatible and bicycle-designated roads for the bike layer;
- contracts each OSM way's intermediate shape points into one directed arc
  between its endpoints, summing haversine segment lengths in kilometres;
- respects `oneway` and basic access tags for cars; bike-usable roads are
  bidirectional by default, including roads that are one-way for motor traffic,
  unless OSM explicitly sets `oneway:bicycle=yes` (or an equivalent true value);
- restricts the result to the largest undirected car-connected component.

The buffered extract is important: nodes inside the city are retained, while
roads outside the boundary are used only to connect the city boundary. For
each inside node adjacent to an outside road, the importer searches the
buffered car graph for paths through outside nodes to another inside node. It
then replaces each such external path by one directed car shortcut whose
distance is the path distance. All road arcs with both endpoints inside the
boundary remain explicitly represented. Bike arcs are kept only when both
endpoints are inside; external shortcuts are car-only.

This is the current bounded-network entry point:

```julia
include("examples/HNDP/aachen_osm_gtfs_generation.jl")
instance = generate_aachen_buffered_osm_gtfs(
    osm_path="data/osm/aachen_buffer_highways.json",
    boundary_path="data/osm/aachen_boundary.json",
    gtfs_dir="data/osm/avv_gtfs",
)
```

For another city, replace the buffered OSM file, boundary relation, and GTFS
directory. The boundary polygon must be represented by an OSM relation with
outer member geometries, as in the Aachen download.

### Fallback when a complete buffered download is unavailable

If the buffered OSM download cannot be obtained, an existing rectangular OSM
extract may be used as a temporary fallback. Supply that extract together
with the administrative boundary and set `use_boundary_compression=true`:

```julia
instance = generate_aachen_osm_gtfs(
    osm_path="data/osm/comparative_cities/heidelberg.json",
    boundary_path="data/osm/comparative_cities/heidelberg_boundary.json",
    gtfs_dir="data/osm/comparative_cities/heidelberg_gtfs",
    use_boundary_compression=true,
)
```

The importer still keeps only internal bike and walking arcs and replaces
available outside-boundary car paths by car-only shortcuts. However, the
shortcuts are limited to the roads present in the rectangular extract, so
detours outside that rectangle can be missing. This option is suitable for
smoke tests and exploratory comparisons, but a complete buffered extract
should be used for final results.

## 3. Assign car, bicycle, and walking speeds

Car speeds are integer operating-speed estimates in km/h. If an OSM
`maxspeed` value is present, it is parsed (including `mph`) and multiplied by
`0.60` to approximate operating speed rather than using the legal maximum.
Otherwise the following road-class defaults are used:

```text
motorway 90    trunk 65       primary 45       secondary 40
tertiary 35    unclassified 30  residential 25 living_street 15
service 20
```

The values are deliberately lower than legal limits for urban roads. This is
consistent with traffic-simulation practice: legal speed and average/free
operating speed are not the same, especially where signals and junctions are
frequent.

Car travel time additionally includes `car_signal_delay_seconds` on explicit
internal car arcs incident to an OSM-tagged traffic-signal node. The default is
15 seconds. The current cached way-only extract does not include separate
traffic-signal nodes, so the parameter has no effect unless an extract with
signal nodes is supplied. External boundary shortcuts do not receive this
increment. The value is configurable and is not a measured
Aachen-wide average. For calibration, the German HBS2015 criteria use average
signalized-intersection delays of up to 20 seconds for LOS A, 35 seconds for
LOS B, and 50 seconds for LOS C. Aachen reports typical signal cycles of 75
seconds in low-load periods and 90 seconds at higher loads. Sources:
[HBS2015 summary](https://www.researchgate.net/publication/329774871_German_Highway_Capacity_Manual_HBS2015),
[City of Aachen signal information](https://www.aachen.de/in-aachen-leben/mobilitaet-verkehr/mobilitaetskonzepte/autoverkehrskonzept/ampeln/).

## 4. Convert time and monetary charges to generalized cost

The generated `mcost` matrix is measured in integer euro-cents. Travel time is
converted using mode-specific values of time, then direct charges are added:

```text
generalized cost [cents] = time [seconds] * VoT [€/hour] * 100 / 3600
                           + direct monetary charge [cents]
```

The default values are car 15 €/h, transit 10 €/h, bike 12 €/h, and walking
9 €/h. They are configurable appraisal estimates; the literature supports
using different values by mode and recommends distinguishing door-to-door,
access, waiting, and in-vehicle time. Transit services have no direct fare in
this Aachen representation because the Deutschlandticket is assumed to cover
them. Walking has no direct fare, but its time still has a value in the
generalized-cost objective.

Bike monetary costs are parameterized explicitly by `bike_entry_cost_cents`
and `bike_cost_per_km_cents`. The default is 61 cents per entry and 0 cents per
kilometre, representing Velocity Aachen's €224 annual Velo365 subscription
divided by 365 days. The subscription includes the first 30 minutes of each
rental. The older subscription switch
`bike_daily_subscription_cost_cents` remains supported for compatibility; set
it to a value to override the explicit entry parameter. The per-use alternative
can be represented with a 50-cent entry charge and an appropriate positive
per-kilometre parameter.
Car exit costs include the existing
parking-search time plus the default local parking cost of 71 cents per day.
This is based on a representative 1.8 m × 4.5 m vehicle: Aachen's
resident-parking fee is vehicle area × 30 €/m²/year plus 15 € administration,
or about 258 €/year and 0.71 €/day. The permit does not guarantee a space.
Set `car_local_parking_cost_cents_per_day=nothing` to use the alternative
visitor tariffs instead: 300 cents inside the 2 km inner-city radius and 100
cents outside it. These approximate Aachen's public
parking tariffs: 3 €/hour in tariff zone I and a 1 € minimum charge for 40
minutes in tariff zone II. Sources: [Aachen resident parking fees](https://serviceportal.aachen.de/suche/-/vr-bis-detail/dienstleistung/7241931/show), [Aachen parking tariffs](https://www.aachen.de/in-aachen-leben/mobilitaet-verkehr/parken/parktarife-strassenrand/),
[European transport appraisal guidance](https://ec.europa.eu/regional_policy/sources/guides/vademecum_2127/vademecum_2127_en.pdf),
[German mode-specific travel-time valuation study](https://www.mdpi.com/2071-1050/11/4/962).

Bike speeds are assigned by OSM infrastructure type, in km/h:

```text
cycleway 20 (at least 19 when a cycleway tag is present)
track 18, path 15, footway/pedestrian 12
living_street 18, residential 17, service 15
primary/secondary 17, tertiary 18, unclassified 16
```

`highway=steps` is excluded from the bike layer. Motorways, motorway links,
trunks, and trunk links are excluded from cycling unless an explicit bicycle
access/designation tag overrides that restriction. Walking uses a separate
5 km/h layer and includes ordinary accessible roads and paths, but excludes
steps, motorways, motorway links, trunks, and trunk links. Walking arcs are
bidirectional because pedestrian travel is not governed by vehicle one-way
rules in this first approximation.

## 4. Convert GTFS to transit services

Stops are matched to the nearest retained OSM node within `0.25` km. For each
GTFS trip, consecutive matched stops are connected even when one or more
intermediate stops lie outside the study area. Thus an in-area stop `A` and an
in-area stop `C` can receive an `A -> C` service using the observed GTFS time
across an outside stop `X`.

For every directed stop pair, the service travel time is the integer median of
all positive observed arrival/departure durations. Services without a valid
positive GTFS duration are ignored. Transit services are fixed arcs on a
separate layer; the original road topology is not copied into that layer.

## 5. Build the physical and mode layers

The generated graph has four layers, each containing one copy of every
retained physical node:

1. physical/walk: the origin/destination layer and walking network;
2. car: original in-boundary roads plus compressed external shortcuts;
3. public transport: fixed GTFS stop-to-stop services;
4. bike: the in-boundary bike-compatible road topology;

For a physical node `i`, the layer indices are `i`, `i+n`, `i+2n`, and
`i+3n`, respectively. Walking arcs are placed directly on the physical layer,
so no artificial zero-time walking-transfer arcs are needed. Car, transit, and
bike transfers connect their mode copies to the physical layer. Bike transfers
are the leader-controlled station decisions.

Car and bike travel times are integer seconds based on the assigned per-arc
operating speed:

```text
car_time  = round(3600 * distance_km / car_operating_speed_kmh)
bike_time = round(3600 * distance_km / bike_operating_speed_kmh)
walk_time = round(3600 * distance_km / 5)
```

All are bounded below by one second. Bike station decisions enable both transfer directions
at a station and are limited by `bike_station_budget`. Bike entry receives
the negative entry fee, and each bike arc receives a negative integer travel
fee based on whole minutes:

```text
bike_arc_profit = fee_cents_per_minute * max(1, ceil(bike_time / 60))
```

The default tariff parameters are `bike_entry_fee_cents=50` and
`bike_travel_fee_cents_per_minute=5`. They are estimates and should be
replaced when a city-specific tariff source is available.

Current transfer-time defaults are:

```text
physical -> car       30 seconds      vehicle assumed immediately available; walk to/access car
car -> physical       90 seconds      parking search / destination access
physical -> bike       60 seconds      locate, unlock, and collect bike
bike -> physical       30 seconds      return, lock, and finish rental
physical -> transit    half expected headway + allocated ticket cost
transit -> physical    0 seconds       alighting is treated as free
```

The default car parking-search value is based on German GPS studies reporting
approximately 90 seconds overall and 1 minute 39 seconds in Frankfurt; a
larger value may be appropriate for central Aachen or peak periods. The
parameter is intentionally exposed because parking search is highly spatially
and temporally variable. See Hagen et al.'s Germany/Frankfurt measurements and
the Munich GPS study in the sources below.

There is no single robust German statistic isolating the time from a person's
front door to a privately owned car parked at home. The current 30-second car
entry value is therefore an explicit modeling assumption for immediate vehicle
availability, while the 90-second exit value represents parking/search and
destination access. These should be separated in later refinements for garage,
driveway, residential-street, and parking-garage settings.

For transit, departures at each matched stop are counted from GTFS. Using a
16-hour service window, the estimated hourly frequency is `departures / 16`,
and the entry time is the expected half-headway:

```text
transit_entry_time = round(1800 / hourly_frequency)
```

The result is capped at `transit_wait_cap_seconds=600` (10 minutes). This
represents a frequency-based approximation in which a user plans departure
to avoid very long waits when the timetable is known. It prevents very low
frequency or one-trip services from creating implausible multi-hour transfer
times.

This is the standard random-arrival approximation for frequency-based transit.
It is less suitable for passengers who plan around an exact timetable or for
very infrequent services, so `transit_service_window_hours` remains a
parameter.

In generalized-cost mode, the Deutschlandticket subscription is allocated
across assumed monthly rides. The defaults are
`transit_monthly_cost_cents=6300` (€63/month) and
`transit_rides_per_month=60` (30 days × 2 rides/day), giving
`transit_trip_cost_cents=105` (€1.05) charged once on physical-to-transit
entry. This is an average subscription allocation, not a single-ride fare,
and all three parameters are configurable. The current €63 price is listed
by [Deutsche Bahn](https://www.bahn.de/angebot/regio/deutschland-ticket).

For the current AVV Aachen feed and the default 16-hour window, after applying
the 10-minute cap, the calculated entry waits over transit-served nodes have a
minimum of 8 seconds, a median of 88 seconds, a mean of 153.8 seconds (2.6
minutes), and a maximum of 600 seconds (10 minutes). The corresponding
uncapped implied headways are twice the values before capping; the capped
maximum is a modeling limit rather than an observed headway.

No reliable general-purpose study was found that isolates only the seconds
needed to unlock and dock a shared bicycle. The current 60/30-second values
are therefore explicit modeling estimates. Studies do show that bike-share
access/egress includes walking plus locating and unlocking/locking, and that
over 90% of metro-bike-share transfers in one dataset occurred within ten
minutes and 300 metres. These values should be replaced if operator-level
observations become available.

## 6. Reproducibility checklist

Record the following with each generated instance:

- city boundary relation and OSM buffered-extract query;
- OSM and GTFS download dates/versions;
- all speed, station, budget, tariff, and user-count parameters;
- random seed;
- source files and importer revision;
- reported node, car-arc, bike-arc, shortcut, transit-service, and user counts.

Run the construction smoke test with:

```powershell
julia --project=. examples/HNDP/aachen_buffered_osm_gtfs_smoke.jl
```

## Sources

- [TU Berlin, OpenStreetMap for traffic simulation](https://svn.vsp.tu-berlin.de/repos/public-svn/publications/vspwp/2011/11-10/2011-06-20_openstreetmap_for_traffic_simulation_sotm-eu.pdf) - distinguishes legal speed limits from average/free operating speeds and reports approximately 25 km/h on Berlin's primary network.
- [OpenStreetMap routing guidance](https://wiki.openstreetmap.org/wiki/Router) - discusses the effects of road class, pavement, curvature, gradients, and intersection conditions on routing speed.
- [OSM maxspeed defaults and road classes](https://wiki.openstreetmap.org/wiki/OSM_tags_for_routing/Maxspeed) - reference for interpreting explicit and implicit speed-limit values.
- [Boufous, Hatfield & Grzebieta, cycling speed on shared paths](https://www.sciencedirect.com/science/article/pii/S000145751730338X) - reports an 18 km/h mean on Sydney shared paths and identifies pedestrian separation as a speed factor.
- [Cycling Speeds in Urban Traffic](https://findingspress.org/article/141204.pdf) - reports infrastructure-specific urban bicycle speeds, including approximately 18.2 km/h on off-street paths and 19.7 km/h on bicycle lanes in its observations.
- [Hagen et al., cruising for parking in Germany](https://www.sciencedirect.com/science/article/pii/S0965856424000934) - GPS-based parking-search measurements, including 1.5 minutes overall and 1:53 in large city centers.
- [Munich parking-search GPS study](https://mediatum.ub.tum.de/1836928) - reports a 90-second mean and 45-second median, with longer durations in central areas.
- [TUM study on residential parking access](https://onlinelibrary.wiley.com/doi/10.1155/2022/5164257) - identifies access convenience, pickup time, and walking distance after parking as distinct residential-car access attributes.
- [Waiting time and headway modelling review](https://www.tandfonline.com/doi/full/10.1080/01441647.2020.1806942) - reviews the half-headway assumption and its limitations.
- [Metro-bike-share transfer study](https://www.mdpi.com/2071-1050/10/5/1526) - defines access/egress transfer time and reports that over 90% of observed transfers were completed within ten minutes and 300 metres.
