"""Symbolic half of the edge-pinning check (numerical half: derivations/verify_edge_pinning_model.jl).
Model problem for the edge-pinned kinematic match: a static, axisymmetric capillary-gravity membrane
    eta'' + eta'/r - eta/l^2 = p(r)   on 0 <= r < oo,  p = 0 for r > a,  eta -> 0 as r -> oo,
indented by a body z_d(r). Contact: eta = z_d on r <= a (kinematic match), eta < z_d for r > a (no overlap),
eta'(a) = z_d'(a) (tangency). Lengths scaled so that l = 1 below (r -> r/l)."""
import sympy as sp, mpmath as mp

r, s, a, Q = sp.symbols("r s a Q", positive=True)
g = sp.Function("g")                     # pressure, arbitrary
I0, K0, I1, K1 = (lambda x: sp.besseli(0, x)), (lambda x: sp.besselk(0, x)), (lambda x: sp.besseli(1, x)), (lambda x: sp.besselk(1, x))

# (1) Green representation of the decaying solution, and that it solves the ODE (needs the Wronskian I0 K1 + I1 K0 = 1/x)
eta_in = -K0(r) * sp.Integral(I0(s) * g(s) * s, (s, 0, r)) - I0(r) * sp.Integral(K0(s) * g(s) * s, (s, r, a))
ode = sp.diff(eta_in, r, 2) + sp.diff(eta_in, r) / r - eta_in - g(r)
ode = sp.simplify(ode.doit())
W = sp.symbols("W")
ode_W = sp.simplify(ode.subs(sp.besseli(0, r) * sp.besselk(1, r), 1 / r - sp.besseli(1, r) * sp.besselk(0, r)))
print("(1) ODE residual of the Green representation, after the Wronskian:", sp.simplify(ode_W))
mp.mp.dps = 40
for x in (mp.mpf("0.01"), mp.mpf("0.4"), mp.mpf(3), mp.mpf(17)):
    assert abs(mp.besseli(0, x) * mp.besselk(1, x) + mp.besseli(1, x) * mp.besselk(0, x) - 1 / x) < mp.mpf(10) ** -35
print("    Wronskian I0 K1 + I1 K0 = 1/x checked to 1e-35 at x = 0.01, 0.4, 3, 17 (DLMF 10.28.2)")
assert sp.simplify(ode_W) == 0

# (2) the edge data depend on the pressure ONLY through Q = int_0^a I0(s) p(s) s ds
eta_a = sp.simplify(eta_in.subs(r, a).doit())
deta_a = sp.simplify(sp.diff(eta_in, r).doit().subs(r, a))
Qint = sp.Integral(I0(s) * g(s) * s, (s, 0, a))
print("(2) eta(a)  =", sp.simplify(eta_a))
print("    eta'(a) =", sp.simplify(deta_a))
assert sp.simplify(eta_a + K0(a) * Qint) == 0
assert sp.simplify(deta_a - K1(a) * Qint) == 0          # d/dr K0 = -K1
# outside the patch the solution is eta = -Q K0(r): a ONE-parameter family

# (3) the two edge conditions E = eta(a) - z_d(a) = 0 and T = eta'(a) - z_d'(a) = 0, written in Q:
zd = sp.Function("z_d")
E = -K0(a) * Q - zd(a)
T = K1(a) * Q - sp.diff(zd(a), a)
Qsol = sp.solve(E, Q)[0]
edge_eq = sp.simplify(T.subs(Q, Qsol))
print("(3) E = 0 and T = 0  <=>  ", sp.Eq(sp.simplify(edge_eq * K0(a)), 0), "   (Q eliminated)")
print("    -> the contact radius a* is fixed by the edge data alone, whatever the interior pressure is.")
print("    -> ANY discretisation that imposes E(a) = 0 exactly and picks a by T(a) = 0 returns a* exactly, at every N")
print("       (provided its linear system is solvable). A Galerkin match that imposes E only weakly does not.")

# (4) the exact solution for a sphere z_d = -D + 1 - sqrt(1 - r^2) (radius 1; l scaled out as rho = 1/l)
rho, D = sp.symbols("rho D", positive=True)          # drop radius in units of l
zs = lambda rr: -D + (1 - sp.sqrt(1 - (rr * rho) ** 2)) / rho   # in l-units
Dstar = sp.solve(sp.Eq(sp.diff(zs(a), a) * K0(a), -K1(a) * zs(a)), D)[0]
print("(4) sphere: D*(a) =", sp.simplify(Dstar))
