#!/usr/bin/env python3
"""Model predictions over the 17 bath experiments, on the axes of Alventosa et al. (2023, figs. 6-7): panels
alpha, t_c/t_sigma, delta/R against We (linear), water on [0, 2.5], oil on [0, 8].

    python3 scripts/plot_score_overlay.py <dir with scripts/score_vs_experiment.jl csvs> <out.svg>

Experiment: black dots, +-1 sd. Each model configuration is a line through its predictions at the experimental We
(dash pattern distinguishes configurations that coincide). No plotting library: plain SVG.
"""
import csv, pathlib, sys

src, out = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
CONFIGS = [  # file stem, legend label, colour, dash (None: markers only)
    ("crossing_minnorm_geometric", "model (default)", "#2a78d6", ""),
    ("tangency_minnorm_geometric", "tangency edge rule", "#eb6834", None),
]
METRICS = [("alpha", "alpha_exp", "alpha_sd"), ("tc", "tc_exp", "tc_sd"), ("delta", "delta_exp", "delta_sd")]
AXES = {"water": (2.5, 0.5, [(0, 0.6, 0.2), (0, 6, 2), (0, 1.5, 0.5)]),
        "oil": (8.0, 2.0, [(0, 0.4, 0.1), (0, 7, 2), (0, 2.5, 0.5)])}
YLAB = {"alpha": "α", "delta": "δ/R",
        "tc": 't<tspan baseline-shift="sub" font-size="8">c</tspan>/t<tspan baseline-shift="sub" font-size="8">σ</tspan>'}
rows = {}
for stem, *_ in CONFIGS:
    with open(src / f"{stem}.csv") as f:
        rows[stem] = list(csv.DictReader(f))

W, H = 760, 470
PL, PT, PW, PH, GX, GY = 52, 14, 190, 150, 62, 64
INK, MUTED = "#111", "#666"
o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
     f'font-family="Helvetica, Arial, sans-serif" font-size="11">', f'<rect width="{W}" height="{H}" fill="#fff"/>']
add = o.append


def diamond(x, y, col):
    return f'<path d="M {x:.1f} {y-4:.1f} l 4 4 l -4 4 l -4 -4 z" fill="none" stroke="{col}" stroke-width="1.5"/>'


def ticks(lo, hi, step):
    t = lo
    while t <= hi + 1e-9:
        yield t
        t = round(t + step, 10)


for r_, fl in enumerate(("water", "oil")):
    xmax, xstep, yax = AXES[fl]
    top = PT + r_ * (PH + GY)
    exp = [r for r in rows[CONFIGS[0][0]] if r["fluid"] == fl]
    for c, (m, me, ms) in enumerate(METRICS):
        x0 = PL + c * (PW + GX)
        ylo, yhi, ystep = yax[c]
        X = lambda w: x0 + w / xmax * PW
        Y = lambda v: top + PH - (v - ylo) / (yhi - ylo) * PH
        add(f'<path d="M {x0} {top} V {top+PH} H {x0+PW}" fill="none" stroke="{INK}" stroke-width="1"/>')
        for t in ticks(ylo, yhi, ystep):
            add(f'<line x1="{x0-4}" y1="{Y(t):.1f}" x2="{x0}" y2="{Y(t):.1f}" stroke="{INK}"/>')
            add(f'<text x="{x0-7}" y="{Y(t)+3.5:.1f}" text-anchor="end" fill="{MUTED}">{t:g}</text>')
        for t in ticks(0.0, xmax, xstep):
            add(f'<line x1="{X(t):.1f}" y1="{top+PH}" x2="{X(t):.1f}" y2="{top+PH+4}" stroke="{INK}"/>')
            add(f'<text x="{X(t):.1f}" y="{top+PH+16}" text-anchor="middle" fill="{MUTED}">{t:g}</text>')
        add(f'<text x="{x0+PW/2}" y="{top+PH+31}" text-anchor="middle" font-style="italic" fill="{INK}">We</text>')
        add(f'<text x="{x0-34}" y="{top+PH/2}" text-anchor="middle" font-size="13" fill="{INK}" '
            f'transform="rotate(-90 {x0-34} {top+PH/2})">{YLAB[m]}</text>')
        if c == 0:
            add(f'<text x="{x0+PW-4}" y="{top+12}" text-anchor="end" fill="{MUTED}">{fl}</text>')
        for stem, label, col, dash in CONFIGS:
            pts = sorted((float(r["We"]), float(r[m])) for r in rows[stem] if r["fluid"] == fl and r["status"] == "ok")
            if dash is None:
                o.extend(diamond(X(w), Y(v), col) for w, v in pts)
            elif len(pts) > 1:
                d = f' stroke-dasharray="{dash}"' if dash else ""
                add('<path d="M ' + " L ".join(f"{X(w):.1f} {Y(v):.1f}" for w, v in pts) +
                    f'" fill="none" stroke="{col}" stroke-width="2"{d}/>')
        for r in exp:
            w, v, s = float(r["We"]), float(r[me]), float(r[ms])
            add(f'<line x1="{X(w):.1f}" y1="{Y(v-s):.1f}" x2="{X(w):.1f}" y2="{Y(v+s):.1f}" stroke="{INK}" stroke-width="1"/>')
            add(f'<circle cx="{X(w):.1f}" cy="{Y(v):.1f}" r="2.6" fill="{INK}"/>')
lx, ly = PL, H - 12
add(f'<circle cx="{lx+4}" cy="{ly-4}" r="2.6" fill="{INK}"/><text x="{lx+12}" y="{ly}" fill="{INK}">experiment</text>')
lx += 90
for stem, label, col, dash in CONFIGS:
    if dash is None:
        add(diamond(lx + 10, ly - 4, col))
    else:
        d = f' stroke-dasharray="{dash}"' if dash else ""
        add(f'<line x1="{lx}" y1="{ly-4}" x2="{lx+22}" y2="{ly-4}" stroke="{col}" stroke-width="2"{d}/>')
    add(f'<text x="{lx+28}" y="{ly}" fill="{INK}">{label}</text>')
    lx += 28 + 6.2 * len(label) + 18
add('</svg>')
out.write_text("\n".join(o))
print("wrote", out)
