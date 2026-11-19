"""Six-city network maps for the appendix (map part of network_figures.py).

Same drawing as network_figures.py, but without the per-panel subtitle (location and site counts
are given in Table 1 of the manuscript); explanatory text belongs in the LaTeX caption.
Writes figures/network_maps.pdf.
"""
from pathlib import Path
import csv, json, math
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor
from paths import B, CODE; C = json.loads((CODE / 'settings.json').read_text())['cities']
out = B / 'figures'; out.mkdir(parents=True, exist_ok=True)
can = canvas.Canvas(str(out / 'network_maps.pdf'), pagesize=(500, 520), pageCompression=1)
for i, city in enumerate(C):
    meta = json.loads((B / 'inputs/city_networks_walk_connected/graphs' / city / 'metadata.json').read_text())
    coords = meta['coords']; n = len(coords); lat0 = sum(x[0] for x in coords) / n; lon0 = sum(x[1] for x in coords) / n
    points = [((x[1] - lon0) * 111.32 * math.cos(math.radians(lat0)), (x[0] - lat0) * 111.32) for x in coords]
    xmin, xmax = min(x[0] for x in points), max(x[0] for x in points); ymin, ymax = min(x[1] for x in points), max(x[1] for x in points)
    x0 = 20 + (i % 2) * 250; y0 = 378 - (i // 2) * 175; w, h = 210, 125
    scale = min(w / (xmax - xmin), h / (ymax - ymin)); ox = x0 + (w - scale * (xmax - xmin)) / 2; oy = y0 + (h - scale * (ymax - ymin)) / 2
    def pos(j):
        x, y = points[j]; return ox + (x - xmin) * scale, oy + (y - ymin) * scale
    can.setStrokeColor(HexColor('#b0bac4')); can.setLineWidth(.20); path = can.beginPath()
    for row in csv.DictReader((B / 'inputs/city_networks_walk_connected/graphs' / city / 'edges.csv').open()):
        u, v = int(row['source']), int(row['target'])
        if 1 <= u < v <= n: path.moveTo(*pos(u - 1)); path.lineTo(*pos(v - 1))
    can.drawPath(path)
    can.setFillColor(HexColor('#e08025'))
    for s in meta['all_candidate_stations']: can.circle(*pos(s - 1), .62, stroke=0, fill=1)
    can.setFillColor(HexColor('#183146')); can.setFont('Helvetica-Bold', 10); can.drawString(x0, y0 + h + 6, city.title())
    bar = 5 if xmax - xmin > 15 else 2
    can.setStrokeColor(HexColor('#183146')); can.setLineWidth(.8); can.line(x0, y0 - 11, x0 + bar * scale, y0 - 11)
    can.setFont('Helvetica', 7); can.drawString(x0, y0 - 21, f'{bar} km')
can.showPage(); can.save(); print('written', out / 'network_maps.pdf')
