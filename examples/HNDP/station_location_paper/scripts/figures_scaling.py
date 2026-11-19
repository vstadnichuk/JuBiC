"""Figures of Section 6.3 (increasing budget and potential location count).

Reads raw/scaling_cases/*/result.json and jobs_scaling.json directly and writes
figures/ and data/scaling/:
  candidate_runtime.pdf  grid of the 54 candidate/budget cases: total time and enumeration time
  station_scaling.pdf    20..300 potential stations at a 25% budget: enumeration time,
                         Gurobi time, paths per user
"""
from pathlib import Path
import csv, json, math
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor, Color

from paths import B, CODE
F = B / "figures"; F.mkdir(parents=True, exist_ok=True)
D = B / "data/scaling"; D.mkdir(parents=True, exist_ok=True)
CITIES = json.loads((CODE / "settings.json").read_text())["cities"]
COLORS = ["#dd7725", "#286c9b", "#8d5f9d", "#409b88", "#bca642", "#c95f6b"]
LIMIT = 3600
CANDIDATES = [20, 80, 140]   # candidate/budget study shown in the paper (equidistant steps of 60)

def load():
    jobs = json.loads((B / "jobs_scaling.json").read_text())
    rows = []
    for j in jobs:
        p = B / "raw/scaling_cases" / j["id"] / "result.json"
        r = json.loads(p.read_text()) if p.exists() else {}
        status = r.get("status", "skipped")
        row = dict(id=j["id"], city=j["city"], candidates=j["candidate_count"], budget=j["station_budget"], status=status)
        if r:
            e, b, m = r["enumeration_seconds"], r["build_seconds"], r["model_seconds"]
            row.update(total_seconds=m, enumeration_seconds=e, assembly_seconds=b - e, solve_seconds=m - b,
                       paths=r["path_count"], users=r["users"], paths_per_user=r["path_count"] / r["users"],
                       modeled_users=r["modeled_path_users"], variables=r["variables"], constraints=r["constraints"],
                       objective=r["objective"], bound=r.get("objective_bound"), gap=r.get("gap"),
                       peak_gib=r.get("peak_private_bytes", 0) / 2**30)
        rows.append(row)
    keys = sorted({k for r in rows for k in r}, key=lambda k: list(rows[-1]).index(k) if k in rows[-1] else 99)
    with (D / "fast_cases.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, keys); w.writeheader(); w.writerows(rows)
    return rows

def label(can, x, y, text, size=8, bold=False, center=False):
    can.setFillColor(HexColor("#1f3445")); can.setFont("Helvetica-Bold" if bold else "Helvetica", size)
    (can.drawCentredString if center else can.drawString)(x, y, str(text))

def shade(secs):
    frac = min(1, math.log10(1 + secs / 60) / math.log10(61))
    return Color(1 - .45 * frac, 1 - .28 * frac, 1 - .10 * frac)

def candidate_runtime(rows):
    cases = {(r["city"], r["candidates"], r["budget"]): r for r in rows}
    cols = [(s, round(s * f)) for s in CANDIDATES for f in [.25, .5, .75]]
    can = canvas.Canvas(str(F / "candidate_runtime.pdf"), pagesize=(510, 204))
    x0, y0, cw, ch = 75, 6, 46, 28
    for j, (s, b) in enumerate(cols): label(can, x0 + j * cw + cw / 2, y0 + 6 * ch + 7, f"{s}/{b}", 8, center=True)
    for i, city in enumerate(CITIES):
        yy = y0 + (5 - i) * ch; label(can, 3, yy + 10, city.title(), 8)
        for j, (s, b) in enumerate(cols):
            r = cases[(city, s, b)]
            if r["status"] == "OPTIMAL":
                can.setFillColor(shade(r["total_seconds"])); top = f"{r['total_seconds'] / 60:.1f}"
            else:  # time limit or failure; enumeration time is still shown if available
                can.setFillColor(HexColor("#e8b9b1")); top = "60+" if r["status"] == "INTERRUPTED" or "time" in r["status"] else r["status"][:6]
            can.rect(x0 + j * cw, yy, cw - 2, ch - 2, stroke=0, fill=1)
            label(can, x0 + j * cw + cw / 2, yy + 14, top, 8, bold=True, center=True)
            if r.get("enumeration_seconds") is not None:
                label(can, x0 + j * cw + cw / 2, yy + 5, f"{r['enumeration_seconds'] / 60:.1f}", 6.5, center=True)
    can.save()

def station_scaling(rows):
    from reportlab.graphics.shapes import Drawing, Line, String, Circle, Rect, Polygon
    from reportlab.graphics import renderPDF, renderSVG
    pts = [r for r in rows if r["budget"] == round(r["candidates"] * .25) and r["status"] != "skipped"]
    with (D / "station_scaling.csv").open("w", newline="") as f:
        keys = ["city", "candidates", "budget", "status", "total_seconds", "enumeration_seconds", "assembly_seconds",
                "solve_seconds", "paths", "paths_per_user", "modeled_users", "variables", "constraints", "gap", "peak_gib"]
        w = csv.DictWriter(f, keys, extrasaction="ignore"); w.writeheader(); w.writerows(sorted(pts, key=lambda r: (CITIES.index(r["city"]), r["candidates"])))
    d = Drawing(510, 212); ink = HexColor("#1f3445")   # title/subtitle are in the manuscript caption
    def text(x, y, s, size=8, bold=False, anchor="start"):
        d.add(String(x, y, str(s), fontName="Helvetica-Bold" if bold else "Helvetica", fontSize=size, fillColor=ink, textAnchor=anchor))
    def line(x1, y1, x2, y2, color="#dfe5e9", width=.5, dash=None):
        d.add(Line(x1, y1, x2, y2, strokeColor=HexColor(color), strokeWidth=width, strokeDashArray=dash))
    def marker(x, y, i, hollow=False):
        color = COLORS[i]
        style = dict(strokeColor=HexColor(color), fillColor=HexColor("#ffffff") if hollow else HexColor(color), strokeWidth=1)
        if i == 0: d.add(Circle(x, y, 2.1, **style))
        elif i == 1: d.add(Rect(x - 2, y - 2, 4, 4, **style))
        elif i == 2: d.add(Polygon([x, y + 2.8, x - 2.5, y - 2, x + 2.5, y - 2], **style))
        elif i == 3: d.add(Polygon([x, y + 2.8, x - 2.5, y, x, y - 2.8, x + 2.5, y], **style))
        elif i == 4: line(x - 2.6, y, x + 2.6, y, color, 1.3); line(x, y - 2.6, x, y + 2.6, color, 1.3)
        else: line(x - 2, y - 2, x + 2, y + 2, color, 1.3); line(x - 2, y + 2, x + 2, y - 2, color, 1.3)
    log = lambda v: math.log10(max(v, 1e-9))
    panels = [  # title, value, lo, hi, ticks, transform, tick format
        ("Enumeration time [min]", lambda r: r["enumeration_seconds"] / 60, 0, 30, [0, 10, 20, 30], lambda v: v, lambda v: f"{v:g}"),
        ("MIP solution time [min]", lambda r: max(r["solve_seconds"], .6) / 60, log(.01), log(100), [.01, .1, 1, 10, 60], log, lambda v: f"{v:g}"),
        ("Paths per user", lambda r: r["paths_per_user"], 1, 5, [1, 2, 3, 4, 5], lambda v: v, lambda v: f"{v:g}"),
    ]
    width, y0, h = 128, 62, 132
    for p, (title, val, lo, hi, ticks, tr, fmt) in enumerate(panels):
        x0 = 36 + p * 168
        def xy(s, v): return x0 + (s - 20) / 280 * width, y0 + (tr(v) - lo) / (hi - lo) * h
        text(x0, 202, title, 8.5, True)
        for v in ticks:
            y = xy(20, v)[1]; line(x0, y, x0 + width, y)
            text(x0 - 5, y - 2.5, fmt(v), 7, anchor="end")
        if p == 1:
            y = xy(20, 60)[1]; line(x0, y, x0 + width, y, "#263746", 1, [4, 3])
        for s in [20, 100, 200, 300]: text(xy(s, 1)[0], 48, s, 7, anchor="middle")
        text(x0 + width / 2, 33, "Potential stations", 8, anchor="middle")
        for i, city in enumerate(CITIES):
            cp = sorted((r for r in pts if r["city"] == city), key=lambda r: r["candidates"])
            xs = [xy(r["candidates"], val(r)) for r in cp]
            for a, z in zip(xs, xs[1:]): line(*a, *z, COLORS[i], 1.1)
            for r, (x, y) in zip(cp, xs): marker(x, y, i, hollow=r["status"] != "OPTIMAL")
    for i, city in enumerate(CITIES):
        x = 12 + i * 83; line(x, 13, x + 10, 13, COLORS[i], 1.2); marker(x + 5, 13, i); text(x + 15, 10, city.title(), 7.5)
    renderPDF.drawToFile(d, str(F / "station_scaling.pdf")); renderSVG.drawToFile(d, str(F / "station_scaling.svg"))

def main():
    rows = load(); candidate_runtime(rows); station_scaling(rows)
    print("written", F)

if __name__ == "__main__": main()
