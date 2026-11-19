from pathlib import Path
import csv,json,math,statistics
from collections import Counter
from paths import B, CODE
C=json.loads((CODE/'settings.json').read_text())
def scc_sizes(n,edges):
 adj=[[] for _ in range(n)];rev=[[] for _ in range(n)]
 for a,z in edges:adj[a].append(z);rev[z].append(a)
 seen=[False]*n;order=[]
 for root in range(n):
  if seen[root]:continue
  seen[root]=True;stack=[(root,0)]
  while stack:
   u,i=stack[-1]
   if i==len(adj[u]):order.append(u);stack.pop();continue
   stack[-1]=(u,i+1);v=adj[u][i]
   if not seen[v]:seen[v]=True;stack.append((v,0))
 seen=[False]*n;sizes=[]
 for root in reversed(order):
  if seen[root]:continue
  seen[root]=True;stack=[root];size=0
  while stack:
   u=stack.pop();size+=1
   for v in rev[u]:
    if not seen[v]:seen[v]=True;stack.append(v)
  sizes.append(size)
 return sorted(sizes,reverse=True)
def stats(city,k,arcs,scope,stations):
 n=4*k;degree=Counter(u for u,v in arcs);counts=Counter();travel={}
 for u,v in arcs:
  lu,lv=(u-1)//k,(v-1)//k
  if lu==lv:
   counts[['walk','car','transit','bike'][lu]]+=1
   travel.setdefault(lu,[]).append(((u-1)%k,(v-1)%k))
  else:counts['transfer']+=1
 sizes=scc_sizes(k,travel.get(0,[]));assert sizes==[k]
 return dict(city=city,scope=scope,physical_nodes=k,layered_vertices=n,eligible_sites=stations,total_arcs=len(arcs),walk_arcs=counts['walk'],car_arcs=counts['car'],transit_arcs=counts['transit'],bike_arcs=counts['bike'],transfer_arcs=counts['transfer'],mean_outdegree=len(arcs)/n,median_outdegree=statistics.median(degree.get(u,0) for u in range(1,n+1)),max_outdegree=max(degree.values()),isolated_outdegree_vertices=sum(u not in degree for u in range(1,n+1)),directed_density=len(arcs)/(n*(n-1)),directed_density_pct=100*len(arcs)/(n*(n-1)),walking_sccs=len(sizes))
rows=[]
for city in C['cities']:
 root=B/'inputs/city_networks_walk_connected/graphs'/city;m=json.loads((root/'metadata.json').read_text());k=m['nnodes']//4
 arcs={(int(x['source']),int(x['target'])) for x in csv.DictReader((root/'edges.csv').open())};sites=m['all_candidate_stations']
 for v in sites:arcs.add((v,v+3*k));arcs.add((v+3*k,v))
 rows.append(stats(city,k,arcs,'full_all_candidates_available',len(sites)))
 for target in C['targets']:
  path=B/'inputs/aggregation'/city/f'graph_{target}_edges.csv'
  if not path.exists():continue
  a={(int(x['source']),int(x['target'])) for x in csv.DictReader(path.open())};rows.append(stats(city,target,a,'aggregated_all_retained_available',target))
with (B/'data/network_topology.csv').open('w',newline='') as f:w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)
print(json.dumps([x for x in rows if x['scope'].startswith('full')],indent=2))
