using Test
if !isdefined(@__MODULE__, :generate_hndp_osm_networks)
    include("hndp_osm_network_generation.jl")
end

function osm_smoke()
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
                "nusers"=>[1,2], "bike_station_count"=>2,"bike_station_budget"=>1)
            config = Dict{String,Any}("instances"=>[spec], "parameter_seeds"=>[17,18],
                "instance_output"=>Dict{String,Any}("folder"=>joinpath(dir,"out"),"graph_format"=>"both"))
            generated = generate_hndp_osm_networks(config)
            @test length(generated)==4
            @test length(unique(g.name for g in generated))==4
            for g in generated
                @test nv(g.instance.mygraph)==16
                @test g.metadata["city"]=="fixture"
                @test g.metadata["transit_services"]==1
                @test length(g.instance.decision_groups)==2
                @test g.instance.decision_budget==1
                @test all(u -> has_path(g.instance.mygraph,u.origin,u.destination),g.instance.users)
                path = joinpath(dir,"out",g.name,"instance.json")
                restored = read_hndp_instance(path)
                @test collect(edges(restored.instance.mygraph))==collect(edges(g.instance.mygraph))
                @test restored.instance.decision_groups==g.instance.decision_groups
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
