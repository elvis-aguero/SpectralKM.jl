#!/usr/bin/env python3
"""
Writes assets/schematic.svg, the schematic of FORMULATION.md (figure 1). Not to scale.

Panel (a): the undeformed drop above the finite cylindrical bath, with the two coordinate systems. Geometry and
symbols only: the boundary conditions are equations (3.4a-d) in the text, not figure annotations. Panel (b): a deformed drop pressed into the surface, with the
pressed surface S(t), its projection C(t) of radius r_c, the contact angle theta_c and the height z_cm.

    python3 derivations/make_schematic.py            # rewrites assets/schematic.svg

The drop in (b) is xi(theta) = 1 + beta_2 P_2(cos theta) + beta_4 P_4(cos theta) outside the patch, drawn from
the same series as eq. (3.28); inside the patch it coincides with the bath surface, as the kinematic match requires.
"""
import math, pathlib

W, H = 1130, 520
out = []
add = out.append

import re as _re

def _markup(s, size):
    """Turn a_x and a_{xy} into real subscripts (tspans with a baseline shift), so no raw underscores are drawn."""
    d = round(0.28 * size, 1)
    parts = _re.split(r"_(\{[^}]*\}|\w)", s)           # text, sub, text, sub, ..., text
    out = parts[0]
    for k in range(1, len(parts), 2):
        out += f'<tspan dy="{d}" font-size="{round(0.72 * size, 1)}" font-style="normal">{parts[k].strip("{}")}</tspan>'
        out += f'<tspan dy="-{d}">{parts[k + 1]}</tspan>'
    return out

def txt(x, y, s, size=17, anchor="middle", style="italic", weight="normal", fill="#111", rot=None):
    t = f' transform="rotate({rot} {x} {y})"' if rot is not None else ""
    s = _markup(s, size)
    add(f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" text-anchor="{anchor}" font-style="{style}" '
        f'font-weight="{weight}" fill="{fill}" stroke="#fff" stroke-width="4" paint-order="stroke" stroke-linejoin="round"{t}>{s}</text>')

def line(x1, y1, x2, y2, w=1.6, col="#111", dash=None, marker=None, mstart=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    m = f' marker-end="url(#{marker})"' if marker else ""
    ms = f' marker-start="url(#{mstart})"' if mstart else ""
    add(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{col}" stroke-width="{w}"{d}{m}{ms}/>')

def path(pts, w=2, col="#111", fill="none", dash=None, close=False):
    d = "M " + " L ".join(f"{x:.1f} {y:.1f}" for x, y in pts) + (" Z" if close else "")
    da = f' stroke-dasharray="{dash}"' if dash else ""
    add(f'<path d="{d}" stroke="{col}" stroke-width="{w}" fill="{fill}"{da} stroke-linejoin="round"/>')

add(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
    f'font-family="Times New Roman, Times, serif">')
add('<defs>'
    '<marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerUnits="userSpaceOnUse" markerWidth="11" markerHeight="11" orient="auto-start-reverse">'
    '<path d="M0 0 L10 5 L0 10 z" fill="#111"/></marker>'
    '<marker id="ahr" viewBox="0 0 10 10" refX="9" refY="5" markerUnits="userSpaceOnUse" markerWidth="11" markerHeight="11" orient="auto-start-reverse">'
    '<path d="M0 0 L10 5 L0 10 z" fill="#c0392b"/></marker>'
    '</defs>')
add(f'<rect width="{W}" height="{H}" fill="#ffffff"/>')

BATH, DROP, WALL, RED = "#cfe3f5", "#fbf3dc", "#222", "#c0392b"

def box(ox, ysurf, ybot, x0=50, x1=470, fill=True):
    """Bath cylinder cross-section (both halves), walls, bottom. fill=False when the surface is not flat."""
    if fill:
        add(f'<rect x="{ox+x0}" y="{ysurf}" width="{x1-x0}" height="{ybot-ysurf}" fill="{BATH}"/>')
    line(ox + x0, ysurf - 40, ox + x0, ybot, w=5, col=WALL)
    line(ox + x1, ysurf - 40, ox + x1, ybot, w=5, col=WALL)
    line(ox + x0, ybot, ox + x1, ybot, w=5, col=WALL)

P2 = lambda x: (3 * x * x - 1) / 2
P4 = lambda x: (35 * x**4 - 30 * x**2 + 3) / 8

# =============================== panel (a): before contact ===============================
ox, ys, yb = 10, 300, 410
box(ox, ys, yb)
xa = ox + 260
line(xa, 40, xa, yb, w=1.2, dash="10 4 2 4")                          # axis of symmetry
line(ox + 50, ys, ox + 470, ys, w=1.4, col="#555", dash="7 5")        # z = 0
R = 85
cy = ys - R - 2          # the 4 px stroke is centred on the circle, so lift it 2 px: the drop touches, never enters
add(f'<circle cx="{xa}" cy="{cy}" r="{R}" fill="{DROP}" stroke="#111" stroke-width="4"/>')
add(f'<circle cx="{xa}" cy="{cy}" r="3" fill="#111"/>')
# spherical frame: downward axis, ray at theta, arc
th = math.radians(50)
ex, ey = xa + R * math.sin(th), cy + R * math.cos(th)
line(xa, cy, xa, cy + R, w=1.2, dash="4 3")
line(xa, cy, ex, ey, w=1.6, marker="ah")
add(f'<path d="M {xa} {cy+40} A 40 40 0 0 0 {xa+40*math.sin(th):.1f} {cy+40*math.cos(th):.1f}" fill="none" stroke="#111" stroke-width="1.4"/>')
txt(xa + 21, cy + 62, "θ", 19)
txt(ex + 30, ey + 10, "ξ(θ, t)", 17, anchor="start")
txt(xa + 8, cy - 8, "COM", 13, anchor="start", style="normal")
txt(xa - 34, cy - 36, "R", 18)
line(xa, cy, xa - R * math.cos(math.radians(40)), cy - R * math.sin(math.radians(40)), w=1.2, col="#555")
txt(xa - 14, cy + 34, "(ρ, σ, ν)", 16, anchor="end")
# cylindrical frame at the pole of the undisturbed surface
line(xa, ys, xa, ys - 34, w=1.8, marker="ah"); line(xa, ys, xa + 40, ys, w=1.8, marker="ah")
txt(xa - 8, ys - 24, "z", 17, anchor="end"); txt(xa + 44, ys + 18, "r", 17, anchor="start")
# velocity and gravity
line(xa + 118, 120, xa + 118, 170, w=2.0, marker="ah"); txt(xa + 130, 150, "V₀", 18, anchor="start")
line(ox + 18, 70, ox + 18, 130, w=2.0, marker="ah"); txt(ox + 32, 106, "g", 18, anchor="start")
# bath annotations: boundary conditions
txt(ox + 110, ys + 26, "(ρ, σ, ν)", 17)
txt(ox + 58, ys - 8, "z = 0", 16, anchor="start")
txt(ox + 410, ys - 8, "η(r, t)", 16, anchor="start")
# dimensions
line(ox + 26, ys, ox + 26, yb, w=1.3, marker="ah", mstart="ah"); txt(ox + 14, (ys + yb) / 2 + 5, "h₀", 17, anchor="end")
line(xa, yb + 34, ox + 470, yb + 34, w=1.3, marker="ah", mstart="ah"); txt((xa + ox + 470) / 2, yb + 56, "b", 17)
line(xa, yb + 20, xa, yb + 40, w=1, col="#555"); line(ox + 470, yb + 20, ox + 470, yb + 40, w=1, col="#555")
txt(ox + 260, H - 14, "(a)  before contact", 18, style="normal", weight="bold")

# =============================== panel (b): during impact ===============================
ox, ys, yb = 570, 270, 425
x0w, x1w = ox + 50, ox + 470
xa = ox + 260
R = 108
b2, b4, thc = -0.20, 0.03, 0.62
xi = lambda t: R * (1 + b2 * P2(math.cos(t)) + b4 * P4(math.cos(t)))
cy = ys - 24                                                          # z_cm = 24 px above the surface
pt = lambda t: (xa + xi(t) * math.sin(t), cy + xi(t) * math.cos(t))
pole, edge = (xa, cy + xi(0.0)), pt(thc)
inside = [pt(thc * k / 40) for k in range(41)]                        # patch: pole -> right edge
rest = [pt(thc + (math.pi - thc) * k / 140) for k in range(141)]       # drop outside the patch: right edge -> top
mir = lambda pts: [(2 * xa - x, y) for x, y in pts]
# Bath surface outside the patch leaves TANGENTIALLY to the drop (the model's smooth-pasting condition) and then
# flattens: its slope decays while the drop's grows, so the bath stays below the drop and never crosses it.
h_ = 1e-4
dx = (pt(thc + h_)[0] - pt(thc - h_)[0]); dy = (pt(thc + h_)[1] - pt(thc - h_)[1])
slope = -dy / dx                                                       # rise (px up) per px outward at the edge
Lm = (edge[1] - ys) / slope                                            # y_b(s) = ys + (y_edge - ys) exp(-s/Lm)
sright = [(edge[0] + sx, ys + (edge[1] - ys) * math.exp(-sx / Lm)) for sx in [k * (x1w - edge[0]) / 120 for k in range(121)]]
surface = mir(sright[::-1]) + mir(inside[::-1]) + inside + sright      # left wall -> ... -> right wall
# bath fill: only BELOW its own surface
add('<path d="M ' + " L ".join(f"{x:.1f} {y:.1f}" for x, y in surface) + f' L {x1w} {yb} L {x0w} {yb} Z" fill="{BATH}"/>')
box(ox, ys, yb, fill=False)
line(xa, 40, xa, yb, w=1.2, dash="10 4 2 4")
line(x0w, ys, x1w, ys, w=1.4, col="#555", dash="7 5")
# drop polygon closed along the pressed surface (not along a chord), so it never overlaps the bath polygon
drop = mir(rest[::-1]) + mir(inside[::-1]) + inside + rest
add('<path d="M ' + " L ".join(f"{x:.1f} {y:.1f}" for x, y in drop) + f' Z" fill="{DROP}" stroke="none"/>')
# geometric audit: outside the patch the drop must be strictly above the bath (SVG y grows downward)
def bath_y(x):
    d = abs(x - xa)
    return ys + (edge[1] - ys) * math.exp(-(d - (edge[0] - xa)) / Lm)
worst = min(bath_y(x) - y for x, y in rest if abs(x - xa) >= (edge[0] - xa) + 1.0 and x <= x1w and x >= x0w)
assert worst > -0.5, f"the drop enters the liquid by {-worst:.1f} px"
assert all(x0w < x < x1w for x, _ in drop), "the drop leaves the cylinder"
path(mir(sright[::-1]), w=2.8, col="#2c6aa0")
path(sright, w=2.8, col="#2c6aa0")
path(mir(rest[::-1]) + rest, w=4, col="#111")                          # drop outline outside the patch
path(mir(inside[::-1]) + inside, w=7, col=RED)                         # pressed surface S(t)
add(f'<circle cx="{xa}" cy="{cy}" r="3" fill="#111"/>')
rc = edge[0] - xa
# C(t): projection of S onto z = 0, drawn on the undisturbed level
line(xa - rc, ys, xa + rc, ys, w=5, col="#666", dash="12 6")
line(xa - rc, ys, xa - rc, edge[1], w=1, col="#aaa"); line(xa + rc, ys, xa + rc, edge[1], w=1, col="#aaa")
txt(xa - R - 18, ys - 30, "C(t)", 16, anchor="end")
line(xa - R - 14, ys - 26, xa - rc - 6, ys - 4, w=1.2, marker="ah")
# r_c dimension below the pole
yd = pole[1] + 44
line(xa, yd, xa + rc, yd, w=1.4, marker="ah", mstart="ah"); txt(xa + rc / 2, yd + 20, "r_c(t)", 16)
line(xa + rc, edge[1] + 3, xa + rc, yd + 6, w=1, col="#888"); line(xa, pole[1] + 3, xa, yd + 6, w=1, col="#888")
# z_cm
line(xa - 30, ys, xa - 30, cy, w=1.4, marker="ah"); txt(xa - 38, cy - 4, "z_{cm}", 16, anchor="end")
# spherical frame: theta_c and xi
line(xa, cy, xa, cy + xi(0.0), w=1.2, dash="4 3")
line(xa, cy, edge[0], edge[1], w=1.4)
add(f'<path d="M {xa} {cy+52} A 52 52 0 0 0 {xa+52*math.sin(thc):.1f} {cy+52*math.cos(thc):.1f}" fill="none" stroke="#111" stroke-width="1.4"/>')
txt(xa + 8, cy + 70, "θ_c", 16, anchor="start")
t2 = 1.35
x2, y2 = pt(t2)
line(xa, cy, x2, y2, w=1.6, marker="ah"); txt(x2 + 10, y2 - 6, "ξ(θ, t)", 16, anchor="start")
txt(xa - 30, cy - 46, "(ρ, σ, ν)", 17)
# labels on the bath side: one label for S(t) and the pressure, to the left of the patch, with a leader that stops short
lx = xa - rc - 30
txt(lx, pole[1] + 18, "S(t)", 17, anchor="end", fill=RED)
line(lx + 6, pole[1] + 10, xa - rc + 16, edge[1] + 7, w=1.2, col=RED, marker="ahr")
txt(ox + 400, ys - 36, "η(r, t)", 16, anchor="start", fill="#2c6aa0")
txt(ox + 260, H - 14, "(b)  during contact", 18, style="normal", weight="bold")
add('</svg>')

pathlib.Path(__file__).resolve().parent.parent.joinpath("assets", "schematic.svg").write_text("\n".join(out))
print("wrote assets/schematic.svg")
