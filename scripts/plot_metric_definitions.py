#!/usr/bin/env python3
"""How t_c and alpha depend on how they are read off the model trajectory (scripts/score_metric_definitions.jl output),
on the axes of Alventosa et al. (2023): centre of mass at z = R (instantaneous), and the experimental procedure on the
top of the drop crossing z = 2R with a 40-frame and an 80-frame parabolic fit at 10,000 fps.

    python3 scripts/plot_metric_definitions.py <metric_defs.csv> <out.svg>
"""
import csv, math, pathlib, sys

rows = list(csv.DictReader(open(sys.argv[1])))
out = pathlib.Path(sys.argv[2])
SERIES = [("com_instant", "centre of mass at z = R", "#2a78d6", ""),
          ("north_fitted", "top at z = 2R, 40-frame fit", "#eb6834", ""),
          ("north80", "top at z = 2R, 80-frame fit", "#eb6834", "6 4")]
COLS = [("alpha", "alpha_exp", "alpha_sd", "α"),
        ("tc", "tc_exp", "tc_sd", 't<tspan baseline-shift="sub" font-size="8">c</tspan>/t<tspan baseline-shift="sub" font-size="8">σ</tspan>')]
AX = {"water": (2.5, 0.5, [(0, 0.6, 0.2), (0, 6, 2)]), "oil": (8.0, 2.0, [(0, 0.4, 0.1), (0, 7, 2)])}
W, H = 520, 470
PL, PT, PW, PH, GX, GY = 52, 14, 190, 150, 62, 64
INK, MUTED = "#111", "#666"
o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
     f'font-family="Helvetica, Arial, sans-serif" font-size="11">', f'<rect width="{W}" height="{H}" fill="#fff"/>']
add = o.append


def ticks(lo, hi, step):
    t = lo
    while t <= hi + 1e-9:
        yield t
        t = round(t + step, 10)


def key(m, s):
    return f"{m}_{s}" if s != "north80" else f"{m}_north80"


for r_, fl in enumerate(("water", "oil")):
    xmax, xstep, yax = AX[fl]
    top = PT + r_ * (PH + GY)
    fr = [r for r in rows if r["fluid"] == fl]
    for c, (m, me, ms, ylab) in enumerate(COLS):
        x0 = PL + c * (PW + GX)
        ylo, yhi, ystep = yax[c]
        X = lambda w: x0 + w / xmax * PW
        Y = lambda v: top + PH - (v - ylo) / (yhi - ylo) * PH
        add(f'<path d="M {x0} {top} V {top+PH} H {x0+PW}" fill="none" stroke="{INK}"/>')
        for t in ticks(ylo, yhi, ystep):
            add(f'<line x1="{x0-4}" y1="{Y(t):.1f}" x2="{x0}" y2="{Y(t):.1f}" stroke="{INK}"/>'
                f'<text x="{x0-7}" y="{Y(t)+3.5:.1f}" text-anchor="end" fill="{MUTED}">{t:g}</text>')
        for t in ticks(0.0, xmax, xstep):
            add(f'<line x1="{X(t):.1f}" y1="{top+PH}" x2="{X(t):.1f}" y2="{top+PH+4}" stroke="{INK}"/>'
                f'<text x="{X(t):.1f}" y="{top+PH+16}" text-anchor="middle" fill="{MUTED}">{t:g}</text>')
        add(f'<text x="{x0+PW/2}" y="{top+PH+31}" text-anchor="middle" font-style="italic" fill="{INK}">We</text>')
        add(f'<text x="{x0-34}" y="{top+PH/2}" text-anchor="middle" font-size="13" fill="{INK}" '
            f'transform="rotate(-90 {x0-34} {top+PH/2})">{ylab}</text>')
        if c == 0:
            add(f'<text x="{x0+PW-4}" y="{top+12}" text-anchor="end" fill="{MUTED}">{fl}</text>')
        for s, label, col, dash in SERIES:
            pts = sorted((float(r["We"]), float(r[key(m, s)])) for r in fr if not math.isnan(float(r[key(m, s)])))
            d = f' stroke-dasharray="{dash}"' if dash else ""
            if len(pts) > 1:
                add('<path d="M ' + " L ".join(f"{X(w):.1f} {Y(v):.1f}" for w, v in pts) +
                    f'" fill="none" stroke="{col}" stroke-width="2"{d}/>')
        for r in fr:
            w, v, sd = float(r["We"]), float(r[me]), float(r[ms])
            add(f'<line x1="{X(w):.1f}" y1="{Y(v-sd):.1f}" x2="{X(w):.1f}" y2="{Y(v+sd):.1f}" stroke="{INK}"/>'
                f'<circle cx="{X(w):.1f}" cy="{Y(v):.1f}" r="2.6" fill="{INK}"/>')
lx, ly = PL - 40, H - 30
add(f'<circle cx="{lx+4}" cy="{ly-4}" r="2.6" fill="{INK}"/><text x="{lx+12}" y="{ly}" fill="{INK}">experiment</text>')
for k, (s, label, col, dash) in enumerate(SERIES):
    x, y = (lx + 90, ly) if k == 0 else (lx + (k - 1) * 230, ly + 18)
    if k == 0:
        x = lx + 90
    d = f' stroke-dasharray="{dash}"' if dash else ""
    add(f'<line x1="{x}" y1="{y-4}" x2="{x+22}" y2="{y-4}" stroke="{col}" stroke-width="2"{d}/>'
        f'<text x="{x+28}" y="{y}" fill="{INK}">{label}</text>')
add('</svg>')
out.write_text("\n".join(o))
print("wrote", out)
