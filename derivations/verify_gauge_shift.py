"""The pressure shift p* = p + C0 applied to film AND ambient everywhere (FORMULATION.md, track A3). Symbolic checks."""
import sympy as sp

r, b, k, x, C, h0 = sp.symbols("r b k x C h_0", positive=True)

# ---- (A) bath: a uniform pressure C over the whole surface 0<=r<=b drives nothing
I = sp.integrate(sp.besselj(0, k * r) * r, (r, 0, b))
assert sp.simplify(I - b * sp.besselj(1, k * b) / k) == 0
print("(A) int_0^b J0(k r) r dr = b J1(k b)/k  -> 0 on every eigenvalue J1(k_m b)=0 (m>=1);"
      " piston m=0 forcing ~ k0 tanh(k0 h0) = 0 at k0=0:", sp.limit(k * sp.tanh(k * h0), k, 0))

# ---- (B) drop, Alventosa's projection Pi_l = int_{-1}^{1} p P_l dx: uniform C excites only l=0, whose forcing (2l+1) l = 0
for l in range(0, 7):
    Pi = sp.integrate(C * sp.legendre(l, x), (x, -1, 1))
    assert (l == 0 and Pi == 2 * C) or (l > 0 and Pi == 0)
print("(B) Alventosa projection: uniform C on the whole sphere -> Pi_l = 2C delta_l0, no dynamics. GAUGE INVARIANT.")

# ---- (C) drop, OUR projection b_l = int p x P_l w dx, f = 2 int p w dx, with w = x on the sphere (r^2 = 1 - x^2)
w = -sp.Rational(1, 2) * sp.diff(1 - x**2, x)
assert sp.simplify(w - x) == 0
f = 2 * sp.integrate(C * w, (x, -1, 1))
bl = {l: sp.integrate(C * x * sp.legendre(l, x) * w, (x, -1, 1)) for l in range(0, 7)}
forcing = {l: sp.simplify(-bl[l] * (2 * l + 1) * l) for l in bl}
print("(C) OUR projection over the whole sphere: f =", f, "; b_l =", bl, "; mode forcing -(2l+1) l b_l =", forcing)
assert f == 0 and forcing[2] != 0
print("    -> uniform C would drive the l=2 mode: our x*w convention (exact only on a polar patch, where w>0)"
      " is NOT gauge invariant beyond the patch. The shift must not be applied through it.")

# ---- (D) what the shift does to the contact law. q = p_film - p_amb is what both interfaces feel (all dynamics depend on q only).
#      Signorini in the gauge-invariant variable:  q >= 0, G >= 0, q G = 0   (G = opening, >= 0 separated)
#      Shift everything by C:  p* = p_film + C, p*_amb = C.  The same law reads  p* >= C, G >= 0, (p* - C) G = 0 :
#      the cone moves with the ambient; it does not widen. Imposing p* >= 0 instead is the law q >= -C:
q, G = sp.symbols("q G")
ps = q + C
law_shifted_signorini = [sp.Ge(ps, C), sp.Eq((ps - C) * G, 0)]          # identical to q>=0, qG=0
assert sp.simplify(law_shifted_signorini[0].lhs - law_shifted_signorini[0].rhs - q) == 0
print("(D) Signorini after the gauge shift: p* >= C and (p* - C) G = 0, i.e. exactly q >= 0, qG = 0 again.")
print("    Imposing p* >= 0 means q >= -C: a NEW physical parameter (an absolute-pressure floor), not a gauge.")
print("    On a separated surface q = 0, i.e. p* = C != 0 with G > 0, so p* G = 0 fails: the law is")
print("    {G = 0, q >= -C} U {G > 0, q = 0}; pair (G=0,q=-C) vs (G=1,q=0) gives (q_a-q_b)(G_a-G_b) =",
      (-C - 0) * (0 - 1), "> 0: non-monotone, so no convex QP (Rockafellar Thm 24.8).")

# ---- (E) physical size of the floor if p* is ABSOLUTE pressure (C = p_atm), in units of 2 sigma / R, R = 0.35 mm
for name, sig in (("water", 0.0722), ("oil 5 cSt", 0.0205)):
    print(f"(E) {name}: p_atm / (2 sigma/R) = {101325 / (2 * sig / 3.5e-4):.0f}")
