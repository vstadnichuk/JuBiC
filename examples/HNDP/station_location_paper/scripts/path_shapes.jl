# Path-catalogue and incumbent diagnostics. W/C/T/B are movement layers.
struct HourPathShape
    movement::UInt8
    visited::UInt8
    transfers::Int32
    fixed_transfers::Int32
    decision_transfers::Int32
    sharing_entries::Int32
    sharing_exits::Int32
    arcs::Int32
    sequence::String
end
function hour_path_shape(path,k,decisions)
    moving=UInt8(0);visited=UInt8(0);transfers=0;fixed=0;controlled=0;entries=0;exits=0
    seq=Int[]
    for (u,v) in path.arcs
        a=cld(u,k);b=cld(v,k)
        visited|=UInt8(1<<(a-1));visited|=UInt8(1<<(b-1))
        isempty(seq) && push!(seq,a)
        if a==b
            moving|=UInt8(1<<(a-1))
        else
            transfers+=1
            (u,v) in decisions ? (controlled+=1) : (fixed+=1)
            push!(seq,b)
            a==1 && b==4 && (entries+=1)
            a==4 && b==1 && (exits+=1)
        end
    end
    names=("W","C","T","B")
    HourPathShape(moving,visited,transfers,fixed,controlled,entries,exits,length(path.arcs),join((names[i] for i in seq),">"))
end
const HOUR_SHAPE_DIMENSIONS=("transfers_total","transfers_fixed","transfers_decision","movement_modes","powered_movement_modes","visited_layers","sharing_entries","sharing_exits","arcs","mode_sequence")
hour_shape_values(s)= (s.transfers,s.fixed_transfers,s.decision_transfers,count_ones(s.movement),count_ones(s.movement&UInt8(14)),count_ones(s.visited),s.sharing_entries,s.sharing_exits,s.arcs,s.sequence)
function hour_merge_shapes!(hist,shapes,car,demand;group=nothing,weights=nothing)
    isempty(shapes) && return
    label=isnothing(group) ? (length(shapes)==1 ? "single_path_user" : "multiple_path_user") : group
    for (i,s) in enumerate(shapes)
        mass=isnothing(weights) ? demand/length(shapes) : demand*weights[i]
        countweight=isnothing(weights) ? 1.0 : weights[i]
        for (metric,value) in zip(HOUR_SHAPE_DIMENSIONS,hour_shape_values(s))
            key=(car,label,metric,string(value));old=get(hist,key,(0.0,0.0))
            hist[key]=(old[1]+countweight,old[2]+mass)
        end
    end
end
function hour_hist_rows(hist)
    [(car_available=key[1],user_group=key[2],metric=key[3],value=key[4],path_weight=value[1],demand_weight=value[2]) for (key,value) in sort!(collect(hist);by=x->x[1])]
end
function hour_user_shape_row(index,user,shapes,complete)
    n=length(shapes);v(f)=n==0 ? (0,0,0.0) : (minimum(f(s) for s in shapes),maximum(f(s) for s in shapes),sum(f(s) for s in shapes)/n)
    t=v(s->s.transfers);f=v(s->s.fixed_transfers);d=v(s->s.decision_transfers);m=v(s->count_ones(s.movement));p=v(s->count_ones(s.movement&UInt8(14)));e=v(s->s.sharing_entries)
    (index=index,user=string(user.uname),car_available=endswith(string(user.uname),"_1"),complete=complete,demand=user.mrisk.demand,paths=n,transfers_min=t[1],transfers_max=t[2],transfers_mean=t[3],fixed_transfers_min=f[1],fixed_transfers_max=f[2],fixed_transfers_mean=f[3],decision_transfers_min=d[1],decision_transfers_max=d[2],decision_transfers_mean=d[3],movement_modes_min=m[1],movement_modes_max=m[2],movement_modes_mean=m[3],powered_modes_min=p[1],powered_modes_max=p[2],powered_modes_mean=p[3],sharing_entries_min=e[1],sharing_entries_max=e[2],sharing_entries_mean=e[3],paths_with_multiple_powered_modes=count(s->count_ones(s.movement&UInt8(14))>=2,shapes))
end
function hour_shape_smoke_tests()
    decisions=Set([(1,31),(31,1),(2,32),(32,2)])
    p=HNDPPath([(1,11),(11,12),(12,2),(2,22),(22,23),(23,3),(3,4)],1.0,0.0,0.0)
    s=hour_path_shape(p,10,decisions)
    @assert s.transfers==4 && s.fixed_transfers==4 && s.decision_transfers==0
    @assert s.movement==UInt8(7) && s.visited==UInt8(7) && s.sequence=="W>C>W>T>W"
    p=HNDPPath([(1,31),(31,32),(32,2)],1.0,0.0,0.0)
    s=hour_path_shape(p,10,decisions)
    @assert s.transfers==2 && s.fixed_transfers==0 && s.decision_transfers==2
    @assert s.sharing_entries==1 && s.sharing_exits==1 && s.movement==UInt8(8) && s.visited==UInt8(9)
    hist=Dict{Tuple{Bool,String,String,String},Tuple{Float64,Float64}}()
    hour_merge_shapes!(hist,[s,s],false,.4)
    @assert sum(r.demand_weight for r in hour_hist_rows(hist) if r.metric=="transfers_total")==.4
    @assert sum(r.path_weight for r in hour_hist_rows(hist) if r.metric=="transfers_total")==2
    true
end
