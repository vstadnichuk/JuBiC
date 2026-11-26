"""Figures of Section 6.2 (aggregation impact).

Reads data/aggregation/ (written by analyze_aggregation.py) and
data/station_full_network_design_replay.csv (written by replay_designs_full_network.jl) and writes
  figures/aggregation_paths_per_user.pdf   paths per user for 50-400 locations
  figures/aggregation_modal.pdf            modal split, aggregated (A) vs. full network (F)
"""
from pathlib import Path
import csv, json, math
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor, Color

from paths import B, CODE
F = B / "figures"
F.mkdir(parents=True, exist_ok=True)
D = B / "data/aggregation"
# Manuscript modal-split colors (plain fills, no patterns).
MC_PAPER = ["#e3c95b", "#8c6ba8", "#4895bd", "#e38735", "#379b86", "#e3e6e9"]
CITIES = json.loads((CODE / "settings.json").read_text())["cities"]
COLORS = ["#dd7725", "#286c9b", "#8d5f9d", "#409b88", "#bca642", "#c95f6b"]
MODES = ["walking", "car", "transit", "bike_sharing", "multi_powered"]
# Mode colors are independent of the city palette. Patterns also identify modes
# when photocopied or viewed without color.
MC = ["#F2E6B1", "#88CCEE", "#CC6677", "#EE7733", "#332288", "#EEEEEE"]

def modal_segment(can, x, y, w, h, index):
    """Vector fill and clipped pattern, shared by bars and legend swatches."""
    if h<=0:return
    can.saveState()
    can.setFillColor(HexColor(MC[index]));can.rect(x,y,w,h,stroke=0,fill=1)
    clip=can.beginPath();clip.rect(x,y,w,h);can.clipPath(clip,stroke=0,fill=0)
    can.setStrokeColor(HexColor('#FFFFFF' if index==4 else '#222222'))
    can.setFillColor(HexColor('#222222'));can.setLineWidth(.45)
    if index==0:  # horizontal
        for v in range(3,math.ceil(h),5):can.line(x,y+v,x+w,y+v)
    elif index in (1,2,4):
        for a in range(-math.ceil(h),math.ceil(w)+1,5):
            if index in (1,4):can.line(x+a,y,x+a+h,y+h)
            if index in (2,4):can.line(x+a,y+h,x+a+h,y)
    elif index==3:  # staggered dots
        for row,v in enumerate(range(2,math.ceil(h),4)):
            for u in range(2+(row%2)*2,math.ceil(w),4):can.circle(x+u,y+v,.55,stroke=0,fill=1)
    can.restoreState()
    can.setLineWidth(.3);can.setStrokeColor(HexColor('#555555'))
    can.rect(x,y,w,h,stroke=1,fill=0)

def read(name):
    p = D / name
    return list(csv.DictReader(p.open())) if p.exists() else []

def number(x, default=0):
    try: return float(x)
    except (ValueError, TypeError): return default

def label(can, x, y, text, size=8, bold=False, center=False):
    can.setFillColor(HexColor("#1f3445"))
    can.setFont("Helvetica-Bold" if bold else "Helvetica", size)
    (can.drawCentredString if center else can.drawString)(x, y, str(text))

def aggregation_paths_figure(cases):
    """Normalize by all OD users; distinguish city effects from quadratic OD growth."""
    from reportlab.graphics.shapes import Drawing, Line, String, Circle, Rect, Polygon
    from reportlab.graphics import renderPDF, renderSVG
    rows=[]
    for r in cases:
        if 'aggregation' not in r['groups'] or not r['path_count']:continue
        k,u,p=int(r['locations']),int(r['users']),int(r['path_count'])
        assert r['status']=='OPTIMAL' and u==2*k*(k-1) and p>=u
        rows.append(dict(city=r['city'],locations=k,users=u,paths=p,
            lower_bound_paths=u,paths_per_user=p/u,excess_paths_per_user=(p-u)/u,
            lower_bound_paths_per_user=1,candidates=int(r['candidates']),budget=int(r['budget'])))
    assert len(rows)==36
    with (D/'aggregation_paths_per_user.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,list(rows[0]));writer.writeheader();writer.writerows(rows)
    # Title and subtitle are given in the manuscript caption, not in the figure.
    d=Drawing(380,232);ink=HexColor('#1f3445')
    def text(x,y,s,size=8,bold=False,anchor='start'):
        d.add(String(x,y,str(s),fontName='Helvetica-Bold' if bold else 'Helvetica',
                     fontSize=size,fillColor=ink,textAnchor=anchor))
    def line(x1,y1,x2,y2,color='#dfe5e9',width=.5,dash=None):
        d.add(Line(x1,y1,x2,y2,strokeColor=HexColor(color),strokeWidth=width,strokeDashArray=dash))
    def marker(x,y,i,color):
        style=dict(strokeColor=HexColor(color),fillColor=HexColor(color),strokeWidth=1)
        if i==0:d.add(Circle(x,y,2.1,**style))
        elif i==1:d.add(Rect(x-2,y-2,4,4,**style))
        elif i==2:d.add(Polygon([x,y+2.8,x-2.5,y-2,x+2.5,y-2],**style))
        elif i==3:d.add(Polygon([x,y+2.8,x-2.5,y,x,y-2.8,x+2.5,y],**style))
        elif i==4:line(x-2.6,y,x+2.6,y,color,1.3);line(x,y-2.6,x,y+2.6,color,1.3)
        else:line(x-2,y-2,x+2,y+2,color,1.3);line(x-2,y+2,x+2,y-2,color,1.3)
    # Paths per user only; the total path count grows with the user count and is
    # reported in the CSV written above.
    x0,width,kmin,kmax,lo,hi=46,316,50,400,.99,1.40
    ticks,kticks=[1,1.1,1.2,1.3,1.4],[50,75,100,200,300,400]
    y0,h=65,144
    def xy(k,v):return x0+(k-kmin)/(kmax-kmin)*width,y0+(v-lo)/(hi-lo)*h
    text(x0,220,'Paths per user',8.5,True)
    for v in ticks:
        y=xy(kmin,v)[1]
        if v!=1:line(x0,y,x0+width,y)
        text(x0-6,y-2.5,f'{v:.1f}',7,anchor='end')
    for k in kticks:text(xy(k,lo)[0],50,k,7,anchor='middle')
    text(x0+width/2,35,'Retained locations',8,anchor='middle')
    for i,city in enumerate(CITIES):
        points=sorted((r['locations'],r['paths_per_user']) for r in rows if r['city']==city and kmin<=r['locations']<=kmax)
        for a,z in zip(points,points[1:]):line(*xy(*a),*xy(*z),COLORS[i],1.2)
        for k,v in points:marker(*xy(k,v),i,COLORS[i])
    # Lower bound of one path per user; drawn last so it stays visible when city counts are close to it.
    bound=[(k,1) for k in kticks]
    for a,z in zip(bound,bound[1:]):line(*xy(*a),*xy(*z),'#263746',1,[4,3])
    for i,city in enumerate(CITIES):
        x=12+i*61;line(x,13,x+10,13,COLORS[i],1.2);marker(x+5,13,i,COLORS[i]);text(x+15,10,city.title(),7.5)
    renderPDF.drawToFile(d,str(F/'aggregation_paths_per_user.pdf'))
    renderSVG.drawToFile(d,str(F/'aggregation_paths_per_user.svg'))

def modal_figure(cases, modes):
    optimal = {(r["city"],int(r["locations"])):r for r in cases
               if "aggregation" in r["groups"] and r["status"]=="OPTIMAL"}
    common=[k for k in [50,75,100,200,300,400] if all((c,k) in optimal for c in CITIES)]
    if not common:return None
    small=min(common);large=300 if 300 in common else max(common)
    shares={}
    for r in modes:
        key=(r["city"],int(r["locations"]))
        if r["candidates"]=="20" and r["budget"]=="10":
            shares.setdefault(key,{})[r["mode"]]=number(r["original_demand_pct"])
    replay={}
    p=B/"data/station_full_network_design_replay.csv"
    if p.exists():
        for r in csv.DictReader(p.open()):
            replay[r["city"],int(r["aggregation_locations"])]={k:number(v) for k,v in r.items() if k.endswith("_pct")}
    can=canvas.Canvas(str(F/"aggregation_modal.pdf"),pagesize=(510,430))
    for i,city in enumerate(CITIES):
        x0=34+(i%2)*255;y0=300-(i//2)*130;w=210;h=88
        label(can,x0,y0+h+12,city.title(),9,bold=True)
        for value in [0,50,100]:
            yy=y0+h*value/100;can.setStrokeColor(HexColor("#dfe5e9"));can.setLineWidth(.4);can.line(x0,yy,x0+w,yy)
            label(can,x0-20,yy-2,value,7)
        for j,(kind,k) in enumerate([("A",small),("A",large),("F",small),("F",large)]):
            xx=x0+12+j*49
            if kind=="A":
                r=optimal[city,k];values=[shares.get((city,k),{}).get(m,0) for m in MODES]+[number(r["intrazonal_pct"])]
            else:
                r=replay.get((city,k))
                if not r:label(can,xx+14,y0+h/2,"n/a",7,center=True);continue
                values=[r[m+"_pct"] for m in MODES]+[0]
            assert abs(sum(values)-100)<.01,(city,k,kind,values)
            current=y0
            for value,color in zip(values,MC_PAPER):
                can.setFillColor(HexColor(color));can.rect(xx,current,28,h*value/100,stroke=0,fill=1);current+=h*value/100
            label(can,xx+14,y0-12,f"{kind}{k}",7,center=True)
    for i,(name,color) in enumerate(zip(["Walk","Car","Transit","Sharing","Multi","Diagonal"],MC_PAPER)):
        xx=15+i*84;can.setFillColor(HexColor(color));can.rect(xx,10,8,6,stroke=0,fill=1);label(can,xx+12,9,name,7)
    can.save()
    return small,large

def main():
    cases=read("solver_cases.csv");modes=read("optimized_modal_shares.csv")
    aggregation_paths_figure(cases);sizes=modal_figure(cases,modes)
    (D/"displayed_aggregation_sizes.json").write_text(json.dumps(dict(sizes=sizes)))
    print("Vector result figures written; common displayed sizes:",sizes)

if __name__=="__main__":main()
