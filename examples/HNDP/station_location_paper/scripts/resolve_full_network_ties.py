"""Assign minimum-cost complete routes: bicycle priority, then car.

Classification always uses ONE feasible optimal route's movement-mode mask;
it never unions modes from different alternatives. Within a priority tier,
prefer fewer non-walking modes, then numeric mode mask for determinism.
"""
from pathlib import Path
from fractions import Fraction
import csv,json,collections,hashlib
from paths import B, CODE
FAMILIES=['walking','car','transit','bike_sharing','multi_powered']
CITIES=['aachen','cologne','heidelberg','karlsruhe','nuremberg','dresden']
def frac(s):return Fraction(str(s).replace('//','/'))
def family(mask):
    nonwalk=mask&14
    return {0:'walking',2:'car',4:'transit',8:'bike_sharing'}.get(nonwalk,'multi_powered')
def select(bits):
    masks=[m for m in range(16) if bits&(1<<m)]
    if not masks:raise ValueError('Unreachable OD')
    choices=([m for m in masks if m&8] or [m for m in masks if m&2]
             or [m for m in masks if m&4] or masks)
    return min(choices,key=lambda m:((m&14).bit_count(),m&14,m))
def write(name,rows):
    with (B/'data'/name).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)
def main():
    # A bike-containing multimodal route stays multimodal. Different tied
    # bike-only and car-only routes must never be combined into a multimodal trip.
    assert family(select((1<<3)|(1<<9)))=='bike_sharing'
    assert family(select((1<<3)|(1<<13)))=='multi_powered'
    assert family(select((1<<9)|(1<<13)))=='bike_sharing'
    assert family(select((1<<1)|(1<<3)))=='car'
    profiles=[];classes=[];patterns=[];moves=[];audit={}
    for city in CITIES:
        folder=B/'raw/full_network_city_profiles'/city
        result=json.loads((folder/'result.json').read_text())
        bands={(r['class'],int(r['band'])):r for r in csv.DictReader((folder/'distance_bands.csv').open())}
        total=collections.defaultdict(Fraction);byclass=collections.defaultdict(lambda:collections.defaultdict(Fraction))
        ties=collections.defaultdict(Fraction);allocation=collections.defaultdict(Fraction)
        for r in csv.DictReader((folder/'band_mode_counts.csv').open()):
            band=bands[r['class'],int(r['band'])];bits=int(r['optimal_mode_set_bits'])
            conditional=frac(band['target_mass'])*Fraction(int(r['pairs']),int(band['original_pairs']))
            mass=frac(band['class_mass'])*conditional
            chosen=select(bits);fam=family(chosen)
            assert bits&(1<<chosen)
            if r['class']=='noncar':assert not chosen&2
            byclass[r['class']][fam]+=conditional;total[fam]+=mass
            old={family(m) for m in range(16) if bits&(1<<m)}
            if len(old)>1:
                key='+'.join(sorted(old));ties[key]+=mass;allocation[fam]+=mass
            else:assert fam==next(iter(old))
        assert sum(total.values())==1
        assert sum(ties.values())==frac(result['exact_shares']['tie'])
        for cl,shares in byclass.items():
            assert sum(shares.values())==1
            for fam in FAMILIES:
                classes.append(dict(city=city,availability_class=cl,family=fam,
                    conditional_share=float(shares[fam]),conditional_exact=str(shares[fam])))
        profiles.append(dict(city=city,physical_nodes=result['physical_nodes'],
            ordered_pairs_per_class=result['ordered_pairs_per_class'],stations=result['available_bike_stations'],
            car_probability_exact=result['car_probability_exact'],
            **{f+'_pct':100*float(total[f]) for f in FAMILIES},original_tie_pct=100*float(sum(ties.values()))))
        for signature,mass in sorted(ties.items(),key=lambda kv:-kv[1]):
            patterns.append(dict(city=city,tied_families=signature,demand_pct=100*float(mass),
                fraction_of_ambiguous_demand= float(mass/sum(ties.values())),exact_mass=str(mass)))
        for fam in FAMILIES:
            old=frac(result['exact_shares'][fam]);assert total[fam]==old+allocation[fam]
            moves.append(dict(city=city,assigned_family=fam,reassigned_demand_pct=100*float(allocation[fam]),
                reassigned_exact=str(allocation[fam])))
        audit[city]=dict(exact_shares={f:str(total[f]) for f in FAMILIES},
            input_histogram_sha256=hashlib.sha256((folder/'band_mode_counts.csv').read_bytes()).hexdigest(),
            counts_unchanged=True,passed=True)
        print(city, {f:round(100*float(total[f]),4) for f in FAMILIES},flush=True)
    write('full_network_city_profiles_resolved.csv',profiles)
    write('full_network_city_profiles_resolved_by_class.csv',classes)
    write('full_network_trip_tie_sources.csv',patterns)
    write('full_network_trip_tie_allocation.csv',moves)
    (B/'data/full_network_tie_resolution.json').write_text(json.dumps(dict(
        policy=__doc__,cities=audit,passed=True),indent=2),encoding='utf-8')
if __name__=='__main__':main()
