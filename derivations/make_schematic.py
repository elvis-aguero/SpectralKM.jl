#!/usr/bin/env python3
"""
Writes assets/schematic.svg, the single schematic of FORMULATION.md (figure 1). Not to scale.

One panel, no repetition: a deformed drop pressed into the surface of a bath in a finite cylinder. It carries every
symbol the text uses: the cylinder (b, h0), the two coordinate systems ((r, z) on the undisturbed surface z = 0,
(xi, theta) at the drop's centre of mass), the interfaces (eta, xi), and the contact quantities (pressed surface S(t),
its projection C(t) of radius r_c, the contact angle theta_c, the height z_cm). Geometry and symbols only: the boundary
conditions are equations (3.4a-d) in the text.

The drop is xi(theta) = 1 + beta_2 P_2(cos theta) + beta_4 P_4(cos theta) outside the patch, as in eq. (3.28), and
coincides with the bath surface on the patch, as the kinematic match requires. The bath leaves the patch tangent to the
drop and then flattens, so the drop never enters the liquid; the script asserts this.

    python3 derivations/make_schematic.py        # rewrites assets/schematic.svg
    python3 derivations/audit_schematic.py       # checks that no label touches a line or another label
"""
import math, pathlib, re as _re

W, H = 600, 500
out = []
add = out.append


def _markup(s, size):
    """a_x and a_{xy} become real subscripts (tspans with a baseline shift): no raw underscores are drawn."""
    d = round(0.28 * size, 1)
    parts = _re.split(r"_(\{[^}]*\}|\w)", s)           # text, sub, text, sub, ..., text
    out_ = parts[0]
    for k in range(1, len(parts), 2):
        out_ += f'<tspan dy="{d}" font-size="{round(0.72 * size, 1)}" font-style="normal">{parts[k].strip("{}")}</tspan>'
        out_ += f'<tspan dy="-{d}">{parts[k + 1]}</tspan>'
    return out_


def txt(x, y, s, size=17, anchor="middle", style="italic", weight="normal", fill="#111"):
    add(f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" text-anchor="{anchor}" font-style="{style}" '
        f'font-weight="{weight}" fill="{fill}" stroke="#fff" stroke-width="4" paint-order="stroke" '
        f'stroke-linejoin="round">{_markup(s, size)}</text>')


def line(x1, y1, x2, y2, w=1.6, col="#111", dash=None, marker=None, mstart=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    m = f' marker-end="url(#{marker})"' if marker else ""
    ms = f' marker-start="url(#{mstart})"' if mstart else ""
    add(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{col}" stroke-width="{w}"{d}{m}{ms}/>')


def path(pts, w=2, col="#111"):
    add('<path d="M ' + " L ".join(f"{x:.1f} {y:.1f}" for x, y in pts) + f'" stroke="{col}" stroke-width="{w}" fill="none" stroke-linejoin="round"/>')


add(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
    f'font-family="Times New Roman, Times, serif">')
add('<defs><marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerUnits="userSpaceOnUse" markerWidth="11" '
    'markerHeight="11" orient="auto-start-reverse"><path d="M0 0 L10 5 L0 10 z" fill="#111"/></marker>'
    '<marker id="ahr" viewBox="0 0 10 10" refX="9" refY="5" markerUnits="userSpaceOnUse" markerWidth="11" '
    'markerHeight="11" orient="auto-start-reverse"><path d="M0 0 L10 5 L0 10 z" fill="#c0392b"/></marker></defs>')
add(f'<rect width="{W}" height="{H}" fill="#ffffff"/>')

BATH, DROP, WALL, RED = "#cfe3f5", "#fbf3dc", "#222", "#c0392b"
P2 = lambda x: (3 * x * x - 1) / 2
P4 = lambda x: (35 * x**4 - 30 * x**2 + 3) / 8

# ---------------------------------------------------------------- geometry
x0w, x1w, xa = 60, 540, 300            # left wall, right wall, axis
ys, yb = 270, 430                      # undisturbed surface z = 0, bottom z = -h0
R = 108
b2, b4, thc = -0.20, 0.03, 0.62
xi = lambda t: R * (1 + b2 * P2(math.cos(t)) + b4 * P4(math.cos(t)))
cy = ys - 24                           # centre of mass, z_cm = 24 px
pt = lambda t: (xa + xi(t) * math.sin(t), cy + xi(t) * math.cos(t))
pole, edge = (xa, cy + xi(0.0)), pt(thc)
inside = [pt(thc * k / 40) for k in range(41)]                       # patch: pole -> right edge
rest = [pt(thc + (math.pi - thc) * k / 140) for k in range(141)]      # outside the patch: right edge -> top
mir = lambda pts: [(2 * xa - x, y) for x, y in pts]
h_ = 1e-4                                                             # slope of the drop at the patch edge
slope = -(pt(thc + h_)[1] - pt(thc - h_)[1]) / (pt(thc + h_)[0] - pt(thc - h_)[0])
Lm = (edge[1] - ys) / slope                                           # bath: tangent at the edge, then flattens
sright = [(edge[0] + s_, ys + (edge[1] - ys) * math.exp(-s_ / Lm)) for s_ in [k * (x1w - edge[0]) / 120 for k in range(121)]]
surface = mir(sright[::-1]) + mir(inside[::-1]) + inside + sright
rc = edge[0] - xa

# ---------------------------------------------------------------- bath, drop
add('<path d="M ' + " L ".join(f"{x:.1f} {y:.1f}" for x, y in surface) + f' L {x1w} {yb} L {x0w} {yb} Z" fill="{BATH}"/>')
line(x0w, ys - 40, x0w, yb, w=5, col=WALL); line(x1w, ys - 40, x1w, yb, w=5, col=WALL); line(x0w, yb, x1w, yb, w=5, col=WALL)
line(xa, 40, xa, yb, w=1.2, dash="10 4 2 4")                          # axis of symmetry
line(x0w, ys, x1w, ys, w=1.4, col="#555", dash="7 5")                 # z = 0
drop = mir(rest[::-1]) + mir(inside[::-1]) + inside + rest            # closed along the patch, not a chord
add('<path d="M ' + " L ".join(f"{x:.1f} {y:.1f}" for x, y in drop) + f' Z" fill="{DROP}" stroke="none"/>')
worst = min((ys + (edge[1] - ys) * math.exp(-(abs(x - xa) - rc) / Lm)) - y
            for x, y in rest if abs(x - xa) >= rc + 1.0 and x0w <= x <= x1w)
assert worst > -0.5, f"the drop enters the liquid by {-worst:.1f} px"
assert all(x0w < x < x1w for x, _ in drop), "the drop leaves the cylinder"
path(mir(sright[::-1]), w=2.8, col="#2c6aa0"); path(sright, w=2.8, col="#2c6aa0")
path(mir(rest[::-1]), w=4); path(rest, w=4)                     # drop outside the patch: two paths, no chord across S(t)
path(mir(inside[::-1]) + inside, w=7, col=RED)                         # pressed surface S(t)
add(f'<circle cx="{xa}" cy="{cy}" r="3" fill="#111"/>')

# ---------------------------------------------------------------- contact quantities
line(xa - rc, ys, xa + rc, ys, w=5, col="#666", dash="12 6")           # C(t): S projected on z = 0
line(xa - rc, ys, xa - rc, edge[1], w=1, col="#aaa"); line(xa + rc, ys, xa + rc, edge[1], w=1, col="#aaa")
txt(xa - R - 18, ys - 30, "C(t)", 16, anchor="end"); line(xa - R - 14, ys - 26, xa - rc - 6, ys - 4, w=1.2, marker="ah")
yd = pole[1] + 44                                                      # r_c below the pole
line(xa, yd, xa + rc, yd, w=1.4, marker="ah", mstart="ah"); txt(xa + rc / 2, yd + 20, "r_c", 16)
line(xa + rc, edge[1] + 3, xa + rc, yd + 6, w=1, col="#888"); line(xa, pole[1] + 3, xa, yd + 6, w=1, col="#888")
line(xa - 30, ys, xa - 30, cy, w=1.4, marker="ah"); txt(xa - 38, cy - 4, "z_{cm}", 16, anchor="end")
lx = xa - rc - 30                                                      # S(t) with a leader that stops short
txt(lx, pole[1] + 18, "S(t)", 17, anchor="end", fill=RED)
line(lx + 6, pole[1] + 10, xa - rc + 16, edge[1] + 7, w=1.2, col=RED, marker="ahr")

# ---------------------------------------------------------------- spherical frame at the centre of mass
line(xa, cy, xa, cy + xi(0.0), w=1.2, dash="4 3")
line(xa, cy, edge[0], edge[1], w=1.4)
add(f'<path d="M {xa} {cy+50} A 50 50 0 0 0 {xa+50*math.sin(thc):.1f} {cy+50*math.cos(thc):.1f}" fill="none" stroke="#111" stroke-width="1.4"/>')
txt(xa + 8, cy + 70, "θ_c", 16, anchor="start")
t2 = 1.35
x2, y2 = pt(t2)
line(xa, cy, x2, y2, w=1.6, marker="ah"); txt(x2 + 10, y2 - 6, "ξ(θ, t)", 16, anchor="start")
txt(xa - 30, cy - 46, "(ρ, σ, ν)", 17)

# ---------------------------------------------------------------- bath labels, frame, dimensions
txt(x0w + 70, yb - 40, "(ρ, σ, ν)", 17)
txt(x1w - 14, ys - 36, "η(r, t)", 16, anchor="end", fill="#2c6aa0")
txt(x0w + 8, ys - 8, "z = 0", 16, anchor="start")
line(xa, 100, xa, 62, w=1.8, marker="ah"); txt(xa + 10, 80, "z", 17, anchor="start")          # frame (r, z): z along the axis
line(x1w - 90, ys, x1w - 40, ys, w=1.8, marker="ah"); txt(x1w - 65, ys - 10, "r", 17)        # r along z = 0
line(x0w - 28, ys, x0w - 28, yb, w=1.3, marker="ah", mstart="ah"); txt(x0w - 38, (ys + yb) / 2 + 5, "h₀", 17, anchor="end")
line(xa, yb + 34, x1w, yb + 34, w=1.3, marker="ah", mstart="ah"); txt((xa + x1w) / 2, yb + 56, "b", 17)
line(xa, yb + 20, xa, yb + 40, w=1, col="#555"); line(x1w, yb + 20, x1w, yb + 40, w=1, col="#555")
line(x0w + 24, 70, x0w + 24, 130, w=2, marker="ah"); txt(x0w + 38, 106, "g", 18, anchor="start")  # gravity
line(xa + 130, 70, xa + 130, 120, w=2, marker="ah"); txt(xa + 142, 100, "V₀", 18, anchor="start")  # impact speed

add('</svg>')
pathlib.Path(__file__).resolve().parent.parent.joinpath("assets", "schematic.svg").write_text("\n".join(out))
print("wrote assets/schematic.svg")
