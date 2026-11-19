using Test
if !isdefined(@__MODULE__, :generate_hndp_osm_networks)
    include("../network_generation/hndp_osm_network_generation.jl")
end

function osm_smoke()
    @testset "OSM numeric speed limits" begin
        for (tag,expected) in [("30",30.0),("50",50.0),(" 50 km/h ",50.0),
            ("30 mph",48.28032),("20 knots",37.04),("12.5",12.5),
            ("12,5",12.5),("50;30",30.0),("30 MPH",48.28032)]
            @test _tagged_speed_kmh(Dict("maxspeed"=>tag)) ≈ expected
        end
        @test _tagged_speed_kmh(Dict("maxspeed"=>30)) == 30.0
        for tag in ["", "none", "walk", "signals", "DE:urban", "0", "-30", "50;unknown", "30oops"]
            @test _tagged_speed_kmh(Dict("maxspeed"=>tag)) === nothing
        end
        @test _tagged_speed_kmh(Dict()) === nothing
        @test _tagged_speed_kmh(Dict("maxspeed"=>missing)) === nothing
        @test _car_operating_speed(Dict("maxspeed"=>"30"),"primary",35) == 15
        @test _car_operating_speed(Dict("maxspeed"=>"unknown"),"primary",35) == 20
        @test _car_movement_cost(60, 1.0, 13, true) == 44
        @test _car_movement_cost(1, 0.25, 13, true) == 6
        @test _car_movement_cost(60, 1.0, 13, false) == 60
    end
    @testset "OSM/GTFS generation and interchange" begin
        mktempdir() do dir
            coords = [(50.0,6.0),(50.0,6.002),(50.002,6.002),(50.002,6.0)]
            point(p) = Dict("lat"=>p[1], "lon"=>p[2])
            ways = [Dict("type"=>"way", "tags"=>Dict("highway"=>"residential"),
                "geometry"=>point.([coords[i],coords[mod1(i+1,4)]])) for i in 1:4]
            osm = joinpath(dir,"roads.json")
            boundary = joinpath(dir,"boundary.json")
            write(osm, JSON.json(Dict("elements"=>ways)))
            ring = [(49.99,5.99),(49.99,6.01),(50.01,6.01),(50.01,5.99),(49.99,5.99)]
            write(boundary, JSON.json(Dict("elements"=>[Dict("members"=>[
                Dict("role"=>"outer","geometry"=>point.(ring))])])) )
            write(joinpath(dir,"stops.txt"), "stop_id,stop_lat,stop_lon\n1,50.0,6.0\n2,50.002,6.002\n")
            write(joinpath(dir,"stop_times.txt"), "trip_id,arrival_time,departure_time,stop_id,stop_sequence\nt1,08:00:00,08:00:00,1,1\nt1,08:02:00,08:02:00,2,2\n")
            spec = Dict{String,Any}("name"=>"fixture", "city_name"=>"fixture",
                "city_center"=>[50.001,6.001], "osm_path"=>osm,
                "boundary_path"=>boundary, "gtfs_dir"=>dir,
                "nusers"=>[1,2], "bike_station_count"=>2,"bike_station_budget"=>1,
                "bike_station_construction_cost"=>7.0)
            config = Dict{String,Any}("instances"=>[spec], "parameter_seeds"=>[17,18],
                "instance_output"=>Dict{String,Any}("folder"=>joinpath(dir,"out"),"graph_format"=>"both"))
            generated = generate_hndp_osm_networks(config)
            # A car-connected trunk endpoint has no walking arc and must be
            # removed by the default generation pipeline, not just by analysis.
            car_only=Dict("type"=>"way","tags"=>Dict("highway"=>"trunk"),
                "geometry"=>point.([coords[3],(50.004,6.002)]))
            cleanup_osm=joinpath(dir,"roads_with_car_only_endpoint.json")
            write(cleanup_osm,JSON.json(Dict("elements"=>vcat(ways,[car_only]))))
            cleaned=generate_osm_gtfs(;osm_path=cleanup_osm,boundary_path=boundary,
                gtfs_dir=dir,nusers=2,city_name="fixture",city_center=(50.001,6.001))
            @test cleaned.metadata["physical_nodes_before_walking_cleanup"]==5
            @test cleaned.metadata["walking_cleanup_removed_nodes"]==1
            @test cleaned.metadata["walking_connectivity_policy"]=="largest_walking_scc"
            @test nv(cleaned.instance.mygraph)==16
            distance_generated=generate_osm_gtfs(;osm_path=osm,boundary_path=boundary,
                gtfs_dir=dir,nusers=20,city_name="fixture",city_center=(50.001,6.001),
                demand_profile=Dict("upper_bounds_km"=>[0.2],"mass_units"=>[1,1]),
                demand_detour_factor=1.0)
            @test distance_generated.metadata["distance_demand"]["band_user_counts"]==[10,10]
            dc=distance_generated.metadata["coords"]
            @test count(u->OSMDemand.distance_km(dc[u.origin],dc[u.destination])<=0.2,
                distance_generated.instance.users)==10
            @test length(generated)==4
            @test length(unique(g.name for g in generated))==4
            for g in generated
                @test g.metadata["car_distance_cost_cents_per_km"]==22
                @test g.metadata["car_tagged_speed_factor"]==0.50
                @test g.metadata["car_untagged_speed_kmh"]==20
                @test nv(g.instance.mygraph)==16
                @test g.metadata["city"]=="fixture"
                @test g.metadata["transit_services"]==1
                @test length(g.instance.decision_groups)==2
                @test g.instance.decision_budget==1
                @test length(g.instance.edgeA)==4
                @test g.metadata["station_construction_cost_convention"]=="entry_only"
                @test g.metadata["bike_station_construction_costs"]==[7.0,7.0]
                for group in g.instance.decision_groups
                    @test all(a in g.instance.edgeA for a in group)
                    @test g.instance.edge_price[first(group)]==7.0
                    @test g.instance.edge_price[last(group)]==0.0
                end
                @test all(u -> has_path(g.instance.mygraph,u.origin,u.destination),g.instance.users)
                path = joinpath(dir,"out",g.name,"instance.json")
                restored = read_hndp_instance(path)
                @test collect(edges(restored.instance.mygraph))==collect(edges(g.instance.mygraph))
                @test restored.instance.decision_groups==g.instance.decision_groups
                @test Set(restored.instance.edgeA)==Set(g.instance.edgeA)
                @test all(restored.instance.edge_price[a]==g.instance.edge_price[a] for a in g.instance.edgeA)
                @test restored.instance.decision_budget==g.instance.decision_budget
                @test restored.instance.users[1].mcost isa SparseMatrixCSC
                @test restored.instance.users[1].mcost==g.instance.users[1].mcost
                @test isfile(joinpath(dirname(path),"instance.gexf"))
                saved = generate_hndp_osm_networks(Dict{String,Any}("saved_instances"=>[path,path]))
                @test length(saved)==1
            end
            @test_throws ArgumentError generate_hndp_osm_networks(Dict{String,Any}("instances"=>[Dict("instance_type"=>"unsupported")]))
            @test_throws ArgumentError generate_hndp_osm_networks(Dict{String,Any}("instances"=>[Dict()]))
            # Release memory-mapped CSV buffers before Windows removes the fixture.
            GC.gc()
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    osm_smoke()
end
