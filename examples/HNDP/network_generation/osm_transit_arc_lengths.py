"""Recover physical lengths for existing GTFS arcs; explicitly label fallbacks.

Uses NumPy for vectorized shape projection. No network requests are made.
Raw node IDs are snapped before applying the cleanup map, preserving the
original saved topology rather than resnapping stops to a different graph.
"""
from pathlib import Path
from collections import defaultdict,Counter
import argparse,bisect,csv,io,json,math,zipfile
import numpy as np
R=6371.0088
def hav(a,b):
    x,y=map(math.radians,a);X,Y=map(math.radians,b)
    return 2*R*math.asin(math.sqrt(min(1,math.sin((X-x)/2)**2+math.cos(x)*math.cos(X)*math.sin((Y-y)/2)**2)))
def number(s):
    try:return float(s)
    except (TypeError,ValueError):return None
def seconds(s):
    try:
        h,m,t=map(int,s.split(':'));return h*3600+m*60+t
    except (ValueError,AttributeError):return -1
class Snapper:
    def __init__(self,coords):
        self.coords=coords;self.grid=defaultdict(list)
        self.dy=.25/110;self.dx=.25/(110*min(math.cos(math.radians(p[0])) for p in coords))
        for i,p in enumerate(coords,1):self.grid[self.key(p)].append(i)
    def key(self,p):return math.floor(p[0]/self.dy),math.floor(p[1]/self.dx)
    def snap(self,p):
        a,b=self.key(p);best=(float('inf'),0)
        for x in (-1,0,1):
            for y in (-1,0,1):
                for i in self.grid.get((a+x,b+y),[]):best=min(best,(hav(p,self.coords[i-1]),i))
        return best[1] if best[0]<=.25 else 0
class Shape:
    def __init__(self,rows):
        rows.sort(key=lambda x:x[0]);lat=np.array([r[1] for r in rows]);lon=np.array([r[2] for r in rows])
        x=np.radians(lat);y=np.radians(lon)
        h=np.sin(np.diff(x)/2)**2+np.cos(x[:-1])*np.cos(x[1:])*np.sin(np.diff(y)/2)**2
        self.cumulative=np.r_[0,np.cumsum(2*R*np.arcsin(np.sqrt(np.clip(h,0,1))))]
        self.latitude=lat;self.longitude=lon
        values=[r[3] for r in rows]
        self.distance=np.array(values,dtype=float) if all(v is not None for v in values) else None
        if self.distance is not None and (np.any(np.diff(self.distance)<0) or self.distance[-1]<=self.distance[0]):self.distance=None
    def position(self,point,reported,minimum):
        if self.distance is not None and reported is not None and self.distance[0]<=reported<=self.distance[-1]:
            pos=float(np.interp(reported,self.distance,self.cumulative))
            if pos+1e-7>=minimum:return pos
        scale=math.cos(math.radians(point[0]));ax=(self.longitude[:-1]-point[1])*scale*math.pi*R/180
        ay=(self.latitude[:-1]-point[0])*math.pi*R/180
        dx=np.diff(self.longitude)*scale*math.pi*R/180;dy=np.diff(self.latitude)*math.pi*R/180
        denominator=dx*dx+dy*dy
        t=np.clip(np.divide(-(ax*dx+ay*dy),denominator,out=np.zeros_like(dx),where=denominator>0),0,1)
        residual=(ax+t*dx)**2+(ay+t*dy)**2
        pos=self.cumulative[:-1]+t*np.diff(self.cumulative)
        residual[pos+1e-7<minimum]=np.inf
        i=int(np.argmin(residual))
        return float(pos[i]) if residual[i]<=.5**2 else None
def ziprows(z,name):
    matching=[n for n in z.namelist() if n.lower().endswith('/'+name) or n.lower()==name]
    if not matching:return iter(())
    return csv.DictReader(io.TextIOWrapper(z.open(matching[0]),encoding='utf-8-sig',newline=''))
def recover(city,raw_root,network_root,output_root):
    print('Transit geometry',city,flush=True)
    raw=raw_root/'sources'/city;original=json.loads((raw_root/'graphs'/city/'metadata.json').read_text())
    clean=network_root/'graphs'/city;meta=json.loads((clean/'metadata.json').read_text());n=len(meta['coords'])
    remap={int(r['original_node']):int(r['new_node']) for r in csv.DictReader((clean/'physical_node_mapping.csv').open())}
    expected={}
    for e in csv.DictReader((clean/'edges.csv').open()):
        u,v=int(e['source']),int(e['target'])
        if 2*n<u<=3*n and 2*n<v<=3*n:expected[u-2*n,v-2*n]=int(e['cost_cents'])
    snap=Snapper(original['coords']);stops={}
    for row in csv.DictReader((raw/'gtfs/stops.txt').open(encoding='utf-8-sig')):
        lat,lon=number(row.get('stop_lat')),number(row.get('stop_lon'))
        if lat is None or lon is None:continue
        node=snap.snap((lat,lon))
        if node:stops[row['stop_id']]=(node,(lat,lon))
    trips=defaultdict(list)
    for row in csv.DictReader((raw/'gtfs/stop_times.txt').open(encoding='utf-8-sig')):
        stop=stops.get(row['stop_id'])
        if stop is None:continue
        trips[row['trip_id']].append((float(row['stop_sequence']),stop[0],stop[1],seconds(row.get('arrival_time')),seconds(row.get('departure_time')),number(row.get('shape_dist_traveled'))))
    with zipfile.ZipFile(raw/'gtfs.zip') as z:
        trip_shapes={}
        for row in ziprows(z,'trips.txt'):
            # The cached merged feed namespaces IDs with its source index.
            # gtfs.zip is source 0 and retains the original, unprefixed IDs.
            for tid in (row['trip_id'],'0:'+row['trip_id']):
                if tid in trips:trip_shapes[tid]=row.get('shape_id','')
        wanted=set(trip_shapes.values())-{''};shape_rows=defaultdict(list)
        for row in ziprows(z,'shapes.txt'):
            if row['shape_id'] in wanted:
                shape_rows[row['shape_id']].append((int(row['shape_pt_sequence']),float(row['shape_pt_lat']),float(row['shape_pt_lon']),number(row.get('shape_dist_traveled'))))
    shapes={sid:Shape(rows) for sid,rows in shape_rows.items() if len(rows)>1}
    observations=defaultdict(list)
    for tid,sequence in trips.items():
        sequence.sort(key=lambda r:r[0]);shape=shapes.get(trip_shapes.get(tid,''));prev=None;minimum=0.0
        for row in sequence:
            pos=shape.position(row[2],row[5],minimum) if shape else None
            if pos is not None:minimum=pos
            if prev is not None and prev[0][1]!=row[1]:
                prior,priorpos=prev;a,b=remap.get(prior[1],0),remap.get(row[1],0)
                duration=(row[3] if row[3]>=0 else row[4])-(prior[4] if prior[4]>=0 else prior[3])
                if (a,b) in expected and duration>0:
                    length=pos-priorpos if pos is not None and priorpos is not None and pos>priorpos else None
                    source='gtfs_shape' if length is not None else 'gtfs_stop_chord_fallback'
                    if length is None:length=hav(prior[2],row[2])
                    if length<=0:length=hav(meta['coords'][a-1],meta['coords'][b-1]);source='snapped_node_chord_fallback'
                    observations[a,b].append((duration,length,source))
            prev=(row,pos)
    rows=[];mismatches=[]
    for (a,b),cost in sorted(expected.items()):
        obs=observations.get((a,b),[])
        if obs:
            durations=sorted(x[0] for x in obs);median=durations[(len(durations)-1)//2]
            if max(1,round(100*median*meta['transit_value_of_time_eur_per_hour']/3600))!=cost:mismatches.append((a,b))
            cohort=[x for x in obs if x[0]==median]
            geometric=[x for x in cohort if x[2]=='gtfs_shape'];cohort=geometric or cohort
            cohort.sort(key=lambda x:(x[1],x[2]));_,length,source=cohort[(len(cohort)-1)//2]
        else:
            length=hav(meta['coords'][a-1],meta['coords'][b-1]);source='unmatched_service_node_chord_fallback'
        rows.append(dict(source=a+2*n,target=b+2*n,length_km=length,length_source=source))
    if mismatches:raise ValueError(f'{city}: {len(mismatches)} saved GTFS costs do not match reconstructed median-time services')
    output_root.mkdir(parents=True,exist_ok=True)
    with (output_root/(city+'.csv')).open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,['source','target','length_km','length_source']);writer.writeheader();writer.writerows(rows)
    report=dict(city=city,transit_arcs=len(rows),length_sources=dict(Counter(r['length_source'] for r in rows)),cost_validation_mismatches=len(mismatches),matched_trip_shapes=len(trip_shapes),loaded_shapes=len(shapes),length_selection='median-duration cohort; prefer geometry on ties; lower median length')
    (output_root/(city+'.json')).write_text(json.dumps(report,indent=2),encoding='utf-8');print(report,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--raw-root',type=Path,required=True);p.add_argument('--network-root',type=Path,required=True);p.add_argument('--output-root',type=Path,required=True);p.add_argument('cities',nargs='+');args=p.parse_args()
    for city in args.cities:recover(city,args.raw_root,args.network_root,args.output_root)
