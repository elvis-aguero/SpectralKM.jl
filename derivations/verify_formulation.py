#!/usr/bin/env python3
"""
Symbolic verification of FORMULATION.md.

Reads the equations out of the markdown (labelled display equations and the inline definitions), translates them
to SymPy with a purpose-built LaTeX reader, and checks that they follow from one another and from first
principles. If an equation in the markdown is edited so that it no longer holds (a sign, a factor of two), the
corresponding check FAILS: the checks are bound to the text, not to a retyped copy of it.

    pip install sympy mpmath
    python3 derivations/verify_formulation.py [path/to/FORMULATION.md]       (exit status 1 on any failure)

What "verified" means here, by kind:
  [SYM]  exact symbolic identity or derivation (SymPy), no numerics.
  [NUM]  numerical check at a sample state (mpmath, 30 digits). Evidence, not proof.
  [BIND] a statement in the text is located and compared with the verified expression.

What is NOT covered: the physical modelling assumptions (linearisation, thin air film, single-point-disk contact
set), the solver (Levenberg-Marquardt, selector robustness), and agreement with experiment. Those are tested by
the Julia test suite and the experiment comparison, not here.

READING CONVENTION. Typeset physics writes F/2\\pi\\sigma R for F/(2\\pi\\sigma R). The reader therefore binds
juxtaposition tighter than '/', i.e. a/b c = a/(b c), and the checks that depend on it say so.
"""
import re
import sys
import pathlib
import itertools
import mpmath as mp
import sympy as sp

HERE = pathlib.Path(__file__).resolve().parent
MD_PATH = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else HERE.parent / "FORMULATION.md"
TEXT = MD_PATH.read_text()

# ----------------------------------------------------------------------------------------------------------
# Extraction from the markdown
# ----------------------------------------------------------------------------------------------------------
def display_equations(text):
    eqs = {}
    for m in re.finditer(r"\$\$\n(.*?)\n\$\$", text, re.S):
        body = m.group(1).strip()
        t = re.search(r"\\qquad\\text\{\(([^)]*)\)\}\s*$", body)
        if t:
            eqs[t.group(1)] = re.sub(r"[.,;]\s*$", "", body[: t.start()].strip())
    return eqs


def inline_math(text):
    return re.findall(r"(?<!\$)\$([^$\n]+)\$(?!\$)", text)


EQ = display_equations(TEXT)
INL = inline_math(TEXT)


def eq(label):
    if label not in EQ:
        raise KeyError(f"equation ({label}) not found in {MD_PATH.name}")
    return EQ[label]


def inl(prefix, must=None, mustnot=None):
    for s in INL:
        s = s.strip()
        if s.startswith(prefix) and (must is None or must in s) and (mustnot is None or mustnot not in s):
            return s
    raise KeyError(f"inline definition starting {prefix!r} not found in {MD_PATH.name}")


# ----------------------------------------------------------------------------------------------------------
# A small LaTeX -> SymPy reader for the vocabulary the markdown uses
# ----------------------------------------------------------------------------------------------------------
_SYMS = {}


def S(name, real_only=False):
    """One Symbol per name across the whole script, so parsed text and hand derivations share objects."""
    if name not in _SYMS:
        _SYMS[name] = sp.Symbol(name, real=True) if real_only else sp.Symbol(name, positive=True)
    return _SYMS[name]


_STATE = re.compile(r"__d{1,2}$")  # \dot / \ddot symbols are real, not positive


def _strip(tex):
    for junk in (r"\,", r"\;", r"\!", r"\ ", r"\Big", r"\big", r"\left", r"\right", r"\quad"):
        tex = tex.replace(junk, " " if junk == r"\ " else "")
    tex = tex.replace(r"\tfrac", r"\frac").replace(r"\dfrac", r"\frac")
    tex = re.sub(r"\\frac(\d|[A-Za-z])(\d|[A-Za-z])", r"\\frac{\1}{\2}", tex)  # \frac2b
    return tex


def _brace(tex, i):
    """tex[i] == '{': return (inner, index after the matching '}')."""
    depth, j = 0, i
    while j < len(tex):
        depth += tex[j] == "{"
        depth -= tex[j] == "}"
        if depth == 0:
            return tex[i + 1 : j], j + 1
        j += 1
    raise ValueError("unbalanced braces in " + tex)


def _expand_macros(tex):
    """\\frac{A}{B} -> (A)/(B), \\sqrt{A} -> \\sqrt(A), ^{A} -> ^(A), _{A} kept for the symbol tokenizer."""
    out, i = "", 0
    while i < len(tex):
        if tex.startswith(r"\frac", i):
            a, i = _brace(tex, i + 5)
            b, i = _brace(tex, i)
            out += "((" + _expand_macros(a) + ")/(" + _expand_macros(b) + "))"
        elif tex.startswith(r"\sqrt", i):
            a, i = _brace(tex, i + 5)
            out += r"\sqrt(" + _expand_macros(a) + ")"
        elif tex[i] == "^" and i + 1 < len(tex) and tex[i + 1] == "{":
            a, i = _brace(tex, i + 1)
            out += "^(" + _expand_macros(a) + ")"
        else:
            out += tex[i]
            i += 1
    return out


_TOK = re.compile(
    r"""\s*(?:
    (?P<num>\d+(?:\.\d+)?)
  | (?P<func>J_[01]\(|\\tanh\(|\\cosh\(|\\sinh\(|\\sqrt\(|\\sin\(|\\cos\()
  | (?P<sym>(?:\\ddot\s*|\\dot\s*|\\bar\s*)?(?:\\mathfrak\s*[A-Za-z]|\\mathrm\{[A-Za-z]+\}|\\[A-Za-z]+|[A-Za-z])
            (?:_\{[^}]*\}|_[A-Za-z0-9])?)
  | (?P<op>[+\-/^()])
    )""",
    re.X,
)

_FUNC = {"J_0(": "J0(", "J_1(": "J1(", r"\tanh(": "tanh(", r"\cosh(": "cosh(", r"\sinh(": "sinh(",
         r"\sqrt(": "sqrt(", r"\sin(": "sin(", r"\cos(": "cos("}


def _symname(tok):
    pre = ""
    m = re.match(r"\\(ddot|dot|bar)\s*", tok)
    if m:
        tok = tok[m.end():]
        pre = m.group(1)
    tok = re.sub(r"\\mathfrak\s*", "frak_", tok)
    tok = re.sub(r"\\mathrm\{([A-Za-z]+)\}", r"\1", tok)
    tok = tok.replace("\\", "").replace("{", "").replace("}", "")
    return {"": tok, "dot": tok + "__d", "ddot": tok + "__dd", "bar": tok.split("_")[0] + "bar_" + "_".join(tok.split("_")[1:])
            if "_" in tok else tok + "bar"}[pre]


def _tokenize(tex):
    toks, pos = [], 0
    tex = tex.strip()
    while pos < len(tex):
        m = _TOK.match(tex, pos)
        if not m or m.end() == pos:
            raise ValueError(f"cannot tokenize {tex[pos:pos+25]!r} in {tex!r}")
        pos = m.end()
        if m.group("num"):
            toks.append(("num", m.group("num")))
        elif m.group("func"):
            toks.append(("func", _FUNC[m.group("func")]))
        elif m.group("sym"):
            toks.append(("sym", _symname(m.group("sym"))))
        else:
            toks.append(("op", m.group("op")))
    return toks


class _Parser:
    """sum := term (('+'|'-') term)* ; term := ['-'] seg ('/' seg)* ; seg := factor+ (juxtaposition)."""

    def __init__(self, toks):
        self.t, self.i = toks, 0

    def peek(self):
        return self.t[self.i] if self.i < len(self.t) else (None, None)

    def take(self):
        tok = self.t[self.i]
        self.i += 1
        return tok

    def sum(self):
        s = self.term()
        while self.peek() in (("op", "+"), ("op", "-")):
            s += self.take()[1] + self.term()
        return s

    def term(self):
        neg = ""
        if self.peek() in (("op", "-"), ("op", "+")):
            neg = self.take()[1]
        s = neg + "(" + self.seg() + ")"
        while self.peek() == ("op", "/"):
            self.take()
            s += "/(" + self.seg() + ")"
        return s

    def seg(self):
        fs = [self.factor()]
        while self.peek()[0] in ("num", "sym", "func") or self.peek() == ("op", "("):
            fs.append(self.factor())
        return "*".join(fs)

    def factor(self):
        a = self.atom()
        if self.peek() == ("op", "^"):
            self.take()
            neg = ""
            if self.peek() == ("op", "-"):
                self.take()
                neg = "-"
            a += "**(" + neg + self.atom() + ")"
        return a

    def atom(self):
        kind, v = self.take()
        if kind == "num":
            return v
        if kind == "sym":
            return v
        if kind == "func":
            inner = self.sum()
            assert self.take() == ("op", ")")
            return v + inner + ")"
        if (kind, v) == ("op", "("):
            inner = self.sum()
            assert self.take() == ("op", ")")
            return "(" + inner + ")"
        raise ValueError(f"unexpected token {kind}:{v}")


def _ns(names):
    ns = {n: S(n, real_only=bool(_STATE.search(n))) for n in names}
    ns.update(J0=lambda x: sp.besselj(0, x), J1=lambda x: sp.besselj(1, x), tanh=sp.tanh, cosh=sp.cosh,
              sinh=sp.sinh, sqrt=sp.sqrt, sin=sp.sin, cos=sp.cos)
    return ns


def tex2expr(tex):
    """Translate one side of an equation (no '=') to a SymPy expression."""
    toks = _tokenize(_expand_macros(_strip(tex)))
    p = _Parser(toks)
    s = p.sum()
    if p.i != len(toks):
        raise ValueError(f"trailing tokens in {tex!r}")
    names = {v for k, v in toks if k == "sym"}
    return sp.sympify(s, locals=_ns(names), evaluate=True)


def tex2eq(tex):
    parts = tex.split("=")
    if len(parts) != 2:
        raise ValueError(f"expected one '=' in {tex!r}")
    return tex2expr(parts[0]), tex2expr(parts[1])


def residual(tex):
    l, r = tex2eq(tex)
    return sp.expand(l - r)


def reduce_bessel(e):
    """Rewrite higher-order Bessel functions that SymPy's derivatives produce in terms of J0, J1."""
    e = e.replace(lambda x: isinstance(x, sp.besselj) and x.args[0] == 2,
                  lambda x: 2 / x.args[1] * sp.besselj(1, x.args[1]) - sp.besselj(0, x.args[1]))
    return e


# ----------------------------------------------------------------------------------------------------------
# Check registry
# ----------------------------------------------------------------------------------------------------------
RESULTS = []


def check(name, kind, where):
    def deco(fn):
        def run():
            try:
                note = fn()
                RESULTS.append((True, name, kind, where, note or ""))
            except Exception as e:  # noqa: BLE001  (a failed binding is a failed check)
                RESULTS.append((False, name, kind, where, f"{type(e).__name__}: {e}"))
        run.__name__ = fn.__name__
        CHECKS.append(run)
        return run
    return deco


CHECKS = []

# common symbols
t, r, z = sp.symbols("t r z", real=True)
k, h, nu, sig, rho, g, b, R, tsig, Bo, Oh = (S(n) for n in
                                             ("k", "h", "nu", "sigma", "rho", "g", "b", "R", "t_sigma", "Bo", "Oh"))


# ---- 1. Bath: separation of variables in the finite cylinder -> (3.9), (3.10) ----------------------------
@check("Bath mode: Laplace + bottom + kinematic + dynamic conditions give (3.10), up to the O(nu^2) term the text drops",
       "SYM", "(3.1)-(3.4), (3.9), (3.10)")
def bath_mode():
    a = sp.Function("a")(t)
    d = sp.Function("d")(t)
    J0 = sp.besselj(0, k * r)
    Z = sp.cosh(k * (z + h)) / (k * sp.sinh(k * h))
    phi = (a.diff(t) + 2 * nu * k**2 * a) * J0 * Z
    lapH = lambda f: f.diff(r, 2) + f.diff(r) / r
    # (3.1) Laplace
    lap = sp.simplify(reduce_bessel(lapH(phi) + phi.diff(z, 2)))
    assert lap == 0, f"Laplace residual {lap}"
    # (3.4) bottom: dphi/dz = 0 at z = -h
    assert sp.simplify(phi.diff(z).subs(z, -h)) == 0
    # wall: d/dr J0(k r) = -k J1(k r), so J1(k b) = 0 is the no-flux condition at every depth
    assert sp.simplify(sp.diff(J0, r) + k * sp.besselj(1, k * r)) == 0
    # (3.2) kinematic: eta_t = phi_z + 2 nu lap(eta) at z=0, with eta = a J0(k r)
    eta = a * J0
    kin = sp.simplify(reduce_bessel(eta.diff(t) - (phi.diff(z).subs(z, 0) + 2 * nu * lapH(eta))))
    assert kin == 0, f"kinematic residual {kin}"
    # (3.3) dynamic: phi_t = -g eta - 2 nu phi_zz + sigma lap(eta)/rho - p/rho, p = d(t) J0(k r), at z=0
    dyn = phi.diff(t) + g * eta + 2 * nu * phi.diff(z, 2) - sig * lapH(eta) / rho + d * J0 / rho
    coeff = sp.simplify(reduce_bessel(dyn.subs(z, 0)) / J0)
    E = sp.expand(sp.simplify(coeff * k * sp.tanh(k * h)))
    # compare with the markdown's (3.10); derivatives a, a', a'' <-> symbols
    A0, A1, A2, D = (S("a_m", True), S("a_m__d", True), S("a_m__dd", True), S("d_m"))
    R310 = residual(eq("3.10"))
    R310 = R310.xreplace({A0: a, A1: a.diff(t), A2: a.diff(t, 2), D: d, S("k_m"): k, S("h_0"): h})
    diff = sp.simplify(sp.expand(E - R310))
    assert sp.simplify(diff - 4 * nu**2 * k**4 * a) == 0, f"unexpected remainder {diff}"
    return "remainder is exactly 4 nu^2 k^4 a (weak-viscosity term)"


@check("Fourier-Bessel basis: norm, orthogonality, volume neutrality, normaliser in (3.21)", "SYM", "(3.21), text of 3.1")
def bessel_basis():
    kk, kn = sp.symbols("k_m k_n", positive=True)
    # d/dr [ r^2/2 (J0^2 + J1^2) ] = r J0^2  =>  int_0^b J0^2 r dr = b^2/2 (J0(kb)^2 + J1(kb)^2)
    F = r**2 / 2 * (sp.besselj(0, kk * r) ** 2 + sp.besselj(1, kk * r) ** 2)
    assert sp.simplify(reduce_bessel(F.diff(r)) - r * sp.besselj(0, kk * r) ** 2) == 0
    # Lommel: d/dr[ r (J0(kn r) dJ0(kk r)/dr - J0(kk r) dJ0(kn r)/dr) ] = (kn^2 - kk^2) r J0(kk r) J0(kn r)
    Jm, Jn = sp.besselj(0, kk * r), sp.besselj(0, kn * r)
    W = r * (Jn * Jm.diff(r) - Jm * Jn.diff(r))
    assert sp.simplify(reduce_bessel(W.diff(r)) - (kn**2 - kk**2) * r * Jm * Jn) == 0
    # boundary term vanishes when J0'(k b) = 0 for both (J0' = -k J1)
    # volume: d/dr [ r J1(k r)/k ] = r J0(k r)
    assert sp.simplify(reduce_bessel((r * sp.besselj(1, kk * r) / kk).diff(r)) - r * Jm) == 0
    # normaliser: markdown prefactor of (3.21) equals 1/norm when J1(k b) = 0
    lhs = EQ["3.21"].split("=", 1)[1].split(r"\int")[0]
    pref = tex2expr(lhs)
    J0b = sp.besselj(0, S("k_m") * b)
    norm = b**2 / 2 * J0b**2
    assert sp.simplify(pref - 1 / norm) == 0, f"prefactor {pref} != 1/norm"
    return "prefactor 2/(b J0(k_m b))^2 = 1/norm"


@check("Dispersion: (3.10) with nu=0, p=0 reproduces the capillary-gravity relation", "SYM", "(3.10)")
def dispersion():
    R310 = residual(eq("3.10"))
    om = sp.symbols("omega", positive=True)
    A0, A1, A2 = S("a_m", True), S("a_m__d", True), S("a_m__dd", True)
    sub = {S("d_m"): 0, nu: 0}
    c = sp.expand(R310.xreplace(sub))
    # a = exp(i om t): a'' = -om^2 a  => coefficient of a: -om^2 + (sigma k^2/rho + g) k tanh(k h) = 0
    char = sp.simplify(c.xreplace({A2: -om**2, A1: 0, A0: 1}))
    target = -om**2 + (sig * S("k_m") ** 2 / rho + g) * S("k_m") * sp.tanh(S("k_m") * S("h_0"))
    assert sp.simplify(char - target) == 0
    return "omega^2 = (g + sigma k^2/rho) k tanh(k h0)"


# ---- 2. Non-dimensionalisation ---------------------------------------------------------------------------
def _scale_state(sym_to_nd, names_scales):
    return sym_to_nd


@check("Non-dimensional bath mode (3.29) is (3.10) under the stated scales (length R, time t_sigma, pressure 2 sigma/R)",
       "SYM+BIND", "(3.10) -> (3.29), (3.17)")
def nondim_bath():
    assert r"\frac{2\sigma}{R}" in EQ["3.17"], "pressure scale 2 sigma / R not present in (3.17)"
    dim = residual(eq("3.10"))
    scale = {
        S("a_m", True): R * S("a_m", True),
        S("a_m__d", True): R / tsig * S("a_m__d", True),
        S("a_m__dd", True): R / tsig**2 * S("a_m__dd", True),
        S("k_m"): S("k_m") / R,
        S("h_0"): S("h_0") * R,
        nu: Oh * R**2 / tsig,
        sig: rho * R**3 / tsig**2,
        g: Bo * R / tsig**2,
        S("d_m"): 2 * (rho * R**3 / tsig**2) / R * S("c_m"),
    }
    dim_nd = dim.xreplace(scale) * tsig**2 / R
    gam = tex2eq(inl(r"\gamma_m="))[1]
    om2 = tex2eq(inl(r"\omega_m^2="))[1]
    nd = residual(eq("3.29")).xreplace({S("gamma_m"): gam, S("omega_m"): sp.sqrt(om2)})
    assert sp.simplify(sp.expand(dim_nd - nd)) == 0, sp.simplify(sp.expand(dim_nd - nd))
    return "gamma_m = 2 Oh k^2, omega_m^2 = (k^2 + Bo) k tanh(k h0), forcing -2 c_m k tanh(k h0)"


@check("Non-dimensional droplet mode (3.30) is (3.14) with Lamb's alpha_l, omega_l; Pi_l scaled by 2 sigma/R",
       "SYM+BIND", "(3.14) -> (3.30)")
def nondim_drop():
    dim = residual(eq("3.14"))
    al = tex2eq(inl(r"\alpha_l="))[1]
    ol = tex2eq(inl(r"\omega_l^2=", must=r"\sigma"))[1]
    l_ = S("l")
    dim = dim.xreplace({S("alpha_l"): al, S("omega_l"): sp.sqrt(ol)})
    dim = dim.xreplace({S("delta_1l"): 0, S("mu"): rho * nu})  # l >= 2
    scale = {
        S("beta_l", True): R * S("beta_l", True),
        S("beta_l__d", True): R / tsig * S("beta_l__d", True),
        S("beta_l__dd", True): R / tsig**2 * S("beta_l__dd", True),
        nu: Oh * R**2 / tsig,
        sig: rho * R**3 / tsig**2,
        S("Pi_l"): 2 * (rho * R**3 / tsig**2) / R * S("b_l"),
    }
    dim_nd = sp.expand(dim.xreplace(scale) * tsig**2 / R)
    gam = tex2eq(inl(r"\gamma_l="))[1]
    om2 = tex2eq(inl(r"\omega_l^2=", mustnot=r"\sigma"))[1]
    nd = residual(eq("3.30")).xreplace({S("gamma_l"): gam, S("omega_l"): sp.sqrt(om2)})
    assert sp.simplify(dim_nd - nd) == 0, sp.simplify(dim_nd - nd)
    return "gamma_l = Oh (2l+1)(l-1), omega_l^2 = l(l-1)(l+2), forcing -(2l+1) l b_l"


@check("Centre of mass: z_cm'' = F/m - g (the beta_1 = -z_cm mode) becomes (3.31) with m = 4/3 pi rho R^3, F = 2 pi sigma R f",
       "SYM+BIND", "text of 3.3, (3.31)")
def com_law():
    # the equation as written in the markdown (section 3.3)
    law = next(s for s in INL if s.strip().startswith(r"\ddot z_{cm}=") and "F/m" in s)
    l, rr = tex2eq(law)
    F, m_ = S("F"), S("m")
    zdd = S("z_cm__dd", True)
    dim = sp.expand(l - rr)
    # the signed physical statement: Alventosa et al. (3.25) is for beta_1 = -z_cm:  beta1'' = -F/m + g
    ref = zdd - (F / m_ - g)  # => z'' = F/m - g
    assert sp.simplify(dim - ref) == 0, \
        f"text says z_cm'' = {sp.solve(dim, zdd)[0]}, but beta_1'' = -F/m + g with beta_1 = -z_cm gives z_cm'' = F/m - g"
    # definitions present in the text
    assert any(s.strip().replace(" ", "").startswith(r"m=\tfrac43\pi\rho R^3".replace(" ", "")) for s in INL), \
        "droplet mass m = 4/3 pi rho R^3 is not defined in the markdown"
    mm = sp.Rational(4, 3) * sp.pi * rho * R**3
    FF = 2 * sp.pi * sig * R * S("f")
    nd = dim.xreplace({F: FF, m_: mm}).xreplace({zdd: R / tsig**2 * zdd, sig: rho * R**3 / tsig**2, g: Bo * R / tsig**2})
    nd = sp.expand(nd * tsig**2 / R)
    ref31 = residual(eq("3.31")).xreplace({S("z_cm__dd", True): zdd})
    assert sp.simplify(nd - ref31) == 0, f"{sp.simplify(nd)} vs {ref31}"
    assert r"f=F/2\pi\sigma R".replace(" ", "") in "".join(q.replace(" ", "") for q in INL), "f = F/(2 pi sigma R) not stated"
    return "z_cm'' = 3/2 f - Bo"


@check("Net force: f = 2 int p w dx with w = -1/2 d(r^2)/dx equals 2 int p r dr (and F = 2 pi sigma R f)",
       "SYM", "(3.5), 3.3")
def force_weight():
    x, xc = sp.symbols("x x_c", positive=True)
    # sample deformed droplet: xi = 1 + b2 P2(x) + b3 P3(x); p = polynomial in x
    b2, b3 = sp.Rational(1, 50), sp.Rational(-1, 80)
    xi = 1 + b2 * sp.legendre(2, x) + b3 * sp.legendre(3, x)
    r2 = xi**2 * (1 - x**2)
    w = -sp.diff(r2, x) / 2
    p = 3 + 2 * x - x**2
    lhs = sp.integrate(sp.expand(p * w), (x, xc, 1))
    # direct: substitute r = r(x): int_0^{r_c} p r dr = int_{x_c}^{1} p (-r dr/dx) dx
    rx = sp.sqrt(r2)
    rhs = sp.integrate(sp.simplify(p * (-rx * sp.diff(rx, x))), (x, xc, 1))
    assert sp.simplify(lhs - rhs) == 0
    # dimensional factor: F = int p_dim dA = 2 pi R^2 (2 sigma/R) int p r dr = 2 pi sigma R (2 int p r dr)
    f_ = 2 * sp.pi * R**2 * (2 * sig / R) / (2 * sp.pi * sig * R)
    assert sp.simplify(f_ - 2) == 0
    assert "f=2\\int_{x_c}^1p\\,w\\,dx" in "".join(s for s in INL).replace(" ", "") or \
        "f=2\\int_{x_c}^1p\\,w\\,dx" in TEXT.replace(" ", ""), "f = 2 int p w dx not found in text"
    return "f = 2 int_{x_c}^1 p w dx = 2 int p r dr"


# ---- 3. Time discretisation ------------------------------------------------------------------------------
@check("BDF2: variable-step coefficients (a, b, c) are exact for 1, t, t^2; markdown a = (1+2s)/(1+s)", "SYM+BIND", "text of 3.6")
def bdf2_coeffs():
    s, dl = sp.symbols("s delta", positive=True)
    A, B, C = sp.symbols("A B C")
    ts = [dl, 0, -dl / s]  # t_{k+1}, t_k, t_{k-1}
    eqs = []
    for f, fp in ((lambda tt: 1, 0), (lambda tt: tt, 1), (lambda tt: tt**2, 2 * dl)):
        eqs.append(sp.Eq((A * f(ts[0]) + B * f(ts[1]) + C * f(ts[2])) / dl, fp))
    sol = sp.solve(eqs, [A, B, C], dict=True)[0]
    assert sp.simplify(sol[A] - (1 + 2 * s) / (1 + s)) == 0
    assert sp.simplify(sol[B] + (1 + s)) == 0 and sp.simplify(sol[C] - s**2 / (1 + s)) == 0
    a_md = tex2eq(inl(r"\mathfrak a="))[1]
    assert sp.simplify(a_md - sol[A].subs(sp.Symbol("delta"), 1) ) == 0 or sp.simplify(a_md.subs(S("s"), s) - sol[A]) == 0
    return "a = (1+2s)/(1+s), b = -(1+s), c = s^2/(1+s)"


@check("Affine response (3.33): eliminating the velocity from BDF2 gives x^{k+1} = history + kappa q^{k+1}", "SYM+BIND", "(3.33)")
def affine():
    x1, y1, Hx, Hy, q = sp.symbols("x1 y1 H_x H_y q")
    a_, dl, gm, om, F = sp.symbols("a delta gamma omega F", positive=True)
    # a x1 + b x0 + c x_-1 = delta y1 ; a y1 + b y0 + c y_-1 = delta(-2 gamma y1 - omega^2 x1 + F q)
    sol = sp.solve([a_ * x1 - dl * y1 - Hx, a_ * y1 + dl * (2 * gm * y1 + om**2 * x1 - F * q) - Hy], [x1, y1], dict=True)[0]
    kappa_derived = sp.simplify(sp.diff(sol[x1], q))
    md = EQ["3.33"].split(r"\qquad")[1]
    kmd = tex2eq(md)[1]
    sub = {S("delta"): dl, S("frak_a"): a_, S("gamma"): gm, S("omega"): om, S("F"): F}
    assert sp.simplify(kmd.xreplace(sub) - kappa_derived) == 0, (kmd, kappa_derived)
    # the history does not multiply q: x1 = (history terms) + kappa q
    assert sp.simplify(sp.diff(sol[x1], Hx).diff(q)) == 0
    return "kappa = delta^2 F/(a^2 + 2 gamma a delta + delta^2 omega^2)"


@check("Sign structure of the slopes: kappa_m < 0 (m>=1), kappa_0 = 0, lambda_l < 0, kappa_cm > 0, D < 0", "SYM",
       "(3.33), (3.35), A3")
def slope_signs():
    dl, a_ = sp.symbols("delta frak_a", positive=True)
    kap = tex2eq(EQ["3.33"].split(r"\qquad")[1])[1]
    gam_m = tex2eq(inl(r"\gamma_m="))[1]
    om2_m = tex2eq(inl(r"\omega_m^2="))[1]
    Fm = sp.simplify(tex2eq(eq("3.29"))[1] / S("c_m"))
    kk = S("k_m")
    km = kap.xreplace({S("gamma"): gam_m, S("omega"): sp.sqrt(om2_m), S("F"): Fm})
    km = sp.simplify(km)
    num, den = sp.fraction(sp.together(km))
    assert den.is_positive, f"denominator not manifestly positive: {den}"
    assert sp.simplify(num / (-dl**2 * 2 * kk * sp.tanh(kk * S("h_0")))) == 1 or num.is_negative, num
    assert sp.limit(sp.tanh(kk * S("h_0")) * kk, kk, 0) == 0  # F_0 = 0 -> kappa_0 = 0
    Fl = sp.simplify(tex2eq(eq("3.30"))[1] / S("b_l"))
    lam = kap.xreplace({S("gamma"): tex2eq(inl(r"\gamma_l="))[1], S("omega"): sp.sqrt(tex2eq(inl(r"\omega_l^2=", mustnot=r"\sigma"))[1]),
                        S("F"): Fl})
    ln, ld = sp.fraction(sp.together(lam))
    lp = sp.Symbol("l_minus_2", nonnegative=True)  # droplet modes start at l = 2
    assert ln.is_negative and sp.expand(ld.xreplace({S("l"): lp + 2})).is_positive
    Fc = sp.diff(tex2eq(eq("3.31"))[1], S("f"))
    kcm = kap.xreplace({S("gamma"): 0, S("omega"): 0, S("F"): Fc})
    assert sp.simplify(kcm - 3 * dl**2 / (2 * a_**2)) == 0 and kcm.is_positive
    return "all four signs follow from F<0 (bath, drop), F=0 (piston), F>0 (centre of mass)"


# ---- 4. Pressure basis -----------------------------------------------------------------------------------
@check("Shifted-Legendre pressure basis: exactly degree n in x, orthogonal, Gram condition number 2N+1 (Hilbert for monomials)",
       "SYM", "(3.17)")
def pressure_basis():
    x, xc, psi = sp.symbols("x x_c psi")
    pt = lambda n: sp.legendre(n, 2 * (x - xc) / (1 - xc) - 1)
    for n in range(0, 9):
        poly = sp.Poly(sp.together(pt(n)).as_numer_denom()[0], x)
        assert poly.degree() == n, (n, poly.degree())
    for n, m in itertools.combinations_with_replacement(range(0, 6), 2):
        val = sp.integrate(sp.legendre(n, 2 * psi - 1) * sp.legendre(m, 2 * psi - 1), (psi, 0, 1))
        assert val == (sp.Rational(1, 2 * n + 1) if n == m else 0)
    for N in range(1, 8):
        gram = sp.diag(*[sp.Rational(1, 2 * n + 1) for n in range(N + 1)])
        assert max(gram.diagonal()) / min(gram.diagonal()) == 2 * N + 1
    hil = sp.Matrix(4, 4, lambda i, j: sp.Rational(1, i + j + 1))  # monomials psi^n, N=3
    cond = float(mp.cond(mp.matrix(hil.tolist())))
    assert cond > 1e4, cond
    return f"monomial Gram at N=3 has cond {cond:.3g} (> 1e4); shifted Legendre has cond 2N+1"


@check("Quadrature: the droplet and centre-of-mass integrands are polynomial in x (degrees reported)", "SYM", "text of 3.3")
def quad_degrees():
    x = sp.symbols("x")
    N_, L_ = 3, 6
    be = sp.symbols("beta2:%d" % (L_ + 1))
    xi = 1 + sum(be[l - 2] * sp.legendre(l, x) for l in range(2, L_ + 1))
    w = -sp.diff(xi**2 * (1 - x**2), x) / 2
    p = sum(sp.Symbol(f"c{n}") * x**n for n in range(N_ + 1))
    deg = lambda e: sp.Poly(sp.expand(e), x).degree()
    d_com = deg(p * w)
    d_bl = max(deg(p * x * sp.legendre(l, x) * w) for l in range(2, L_ + 1))
    d_zd = deg(xi * x * sp.legendre(N_, x) * w)
    assert d_com == N_ + 2 * L_ + 1
    return f"N=3, L=6: COM integrand degree {d_com} (= N+2L+1); b_l integrand {d_bl}; z_d-Galerkin {d_zd} (= N+3L+2)"


# ---- 5. Contact closure ----------------------------------------------------------------------------------
@check("Tangency is degenerate at the poles: C(theta) is even, so dC/dtheta = 0 at theta = 0 and pi for any state", "SYM",
       "A (tangency remark)")
def tangency_degenerate():
    th = sp.symbols("theta", real=True)
    Mn, Ln = 3, 3
    am = sp.symbols("a0:%d" % Mn)
    be = sp.symbols("b2:%d" % (Ln + 1))
    km = sp.symbols("k0:%d" % Mn, positive=True)
    zc = sp.Symbol("z_cm")
    xi = 1 + sum(be[l - 2] * sp.legendre(l, sp.cos(th)) for l in range(2, Ln + 1))
    rr = xi * sp.sin(th)
    eta = sum(am[m] * sp.besselj(0, km[m] * rr) for m in range(Mn))
    C = eta - zc + xi * sp.cos(th)
    dC = sp.diff(C, th)
    assert sp.simplify(dC.subs(th, 0)) == 0
    assert sp.simplify(dC.subs(th, sp.pi)) == 0
    assert sp.simplify(C.subs(th, -th) - C) == 0  # even
    return "C(-theta) = C(theta)"


@check("Convex form: with C = C_free + A p and Q = -A symmetric, the gradient of (3.26b) is -C, so its KKT conditions are (3.26a)",
       "SYM", "(3.26a), (3.26b)")
def kkt():
    n = 3
    Qs = sp.Matrix(n, n, lambda i, j: sp.Symbol(f"q{min(i,j)}{max(i,j)}"))
    pv = sp.Matrix(n, 1, lambda i, j: sp.Symbol(f"p{i}"))
    Cf = sp.Matrix(n, 1, lambda i, j: sp.Symbol(f"C{i}"))
    J = (pv.T * Qs * pv)[0] / 2 - (pv.T * Cf)[0]
    grad = sp.Matrix([sp.diff(J, pv[i]) for i in range(n)])
    Cexpr = Cf - Qs * pv  # C = C_free + A p, A = -Q
    assert sp.simplify(grad + Cexpr) == sp.zeros(n, 1)
    # the text states (3.26b)
    flat = EQ["3.26b"].replace(" ", "").replace("\\ ", "")
    for piece in (r"\min_{p\ge0}", r"\tfrac12\langle p,(-\mathcal A)p\rangle_w", r"-\langle p,C_{\rm free}\rangle_w"):
        assert piece.replace(" ", "") in flat, f"(3.26b) in the text does not contain {piece}: it is not (1/2)<p,Qp> - <p,C_free> with Q = -A"
    return "grad J = -C  =>  p>0: C=0;  p=0: C<=0 (i.e. -C>=0)"


# ---- 6. Wall conditions ----------------------------------------------------------------------------------
def _indexed(expr, M, skip0=False):
    """Sum expr over m = 0..M-1, where symbols whose name contains '_m' become per-mode symbols;
    mode 0 is the piston: k_0 = 0, kappa_0 = 0, J0(k_0 b) = 1."""
    total = 0
    for i in range(M):
        sub = {}
        for s_ in expr.free_symbols:
            if "_m" in s_.name:
                sub[s_] = S(s_.name.replace("_m", f"_m{i}"))
        e = expr.xreplace(sub)
        if i == 0:
            e = e.xreplace({S("kappa_m0"): 0})
        total += e
    return total


@check("Pinned wall (3.35): Lambda enforces sum a_m J0(k_m b) = 0 exactly; line-force coefficient 2 Lambda/(b J0) matches (3.21)",
       "SYM+BIND", "(3.35), (3.21)")
def pinning():
    clean = lambda q: re.sub(r"[.,;]\s*$", "", q.strip())
    a_eq, lam_eq, d_eq = (clean(q) for q in EQ["3.35"].split(r"\qquad"))
    M = 5
    Jsym = sp.besselj(0, S("k_m") * b)
    Jgen = S("Jb_m")

    def idx(e, i):
        e = e.xreplace({Jsym: Jgen})
        e = e.xreplace({s_: S(s_.name.replace("_m", f"_m{i}")) for s_ in e.free_symbols if "_m" in s_.name})
        if i == 0:  # piston mode: k_0 = 0, J0 = 1, kappa_0 = 0
            e = e.xreplace({S("Jb_m0"): sp.Integer(1), S("kappa_m0"): sp.Integer(0), S("k_m0"): sp.Integer(0)})
        return e

    am = tex2eq(a_eq)[1]
    mnum = re.search(r"\\frac\{\\sum_m(.*)\}\{D\}", lam_eq)
    assert mnum, "cannot read the numerator of Lambda in (3.35)"
    num = tex2expr(mnum.group(1))
    assert lam_eq.replace(" ", "").startswith(r"\Lambda=-"), "Lambda should carry a leading minus sign"
    dfac_tex, dsum_tex = d_eq.split("=", 1)[1].split(r"\sum_m")
    Dval = tex2expr(dfac_tex) * sum(idx(tex2expr(dsum_tex), i) for i in range(M))
    Lam = -sum(idx(num, i) for i in range(M)) / Dval
    cons = sum(idx(am, i).xreplace({S("Lambda"): Lam}) * (idx(Jsym, i) if i else 1) for i in range(M))
    cons = sp.simplify(cons.xreplace({S(f"Jb_m{i}"): S(f"Jb_m{i}") for i in range(M)}))
    assert cons == 0, f"constraint residual {cons}"
    # p = Lambda delta(r-b): c_m = (prefactor of (3.21)) * Lambda * b J0(k_m b)  vs the line-force term of (3.35)
    pref = tex2expr(EQ["3.21"].split("=", 1)[1].split(r"\int")[0])
    cLam = pref * S("Lambda") * b * Jsym
    line = tex2expr(r"(2/b)\,\Lambda/J_0(k_mb)")
    assert sp.simplify(cLam - line) == 0
    assert sp.simplify(am.coeff(S("Lambda"))) is not None
    return "sum_m a_m J0(k_m b) = 0 for 5 modes (piston included); c_m^Lambda = 2 Lambda/(b J0(k_m b))"


@check("Pinned wall: the rank-one correction keeps the response operator self-adjoint in the pairing w_m ~ J0(k_m b)^2",
       "SYM", "(3.35), C2")
def pinning_selfadjoint():
    M = 4
    kap = [S(f"kappa_m{i}") for i in range(M)]
    J = [S(f"Jb_m{i}") for i in range(M)]
    c = [sp.Symbol(f"c{i}") for i in range(M)]
    D = 2 / b * sum(kap)
    Lam = -sum(kap[i] * c[i] * J[i] for i in range(M)) / D
    a = [kap[i] * c[i] + 2 / b * kap[i] * Lam / J[i] for i in range(M)]
    B = sp.Matrix(M, M, lambda i, j: sp.diff(a[i], c[j]))
    WB = sp.Matrix(M, M, lambda i, j: J[i] ** 2 * B[i, j])
    assert sp.simplify(WB - WB.T) == sp.zeros(M, M)
    # and the constraint direction is annihilated: sum_m J_m a_m = 0 for every c
    assert sp.simplify(sum(J[i] * a[i] for i in range(M))) == 0
    return "W B is symmetric (4 modes, generic kappa, J)"


# ---- 7. Compliance operator: symmetry and sign at a sample state ----------------------------------------
@check("Compliance operator: the Galerkin matrix <P_n, A P_n'>_w built from the markdown's definitions is symmetric and "
       "negative semidefinite at a sample state", "NUM+BIND", "A3, 3.6")
def compliance():
    mp.mp.dps = 30
    flat = TEXT.replace(" ", "")
    for needle in (r"c_m=\frac{2}{(bJ_0(k_mb))^2}\int p\,J_0(k_mr)\,r\,dr", r"b_l=\int_{x_c}^1p\,x\,P_l\,w\,dx",
                   r"f=2\int_{x_c}^1p\,w\,dx", r"w=-\tfrac12\,d(r^2)/dx"):
        assert needle.replace(" ", "") in flat, f"definition {needle} not found in the text; the sample-state check encodes it"
    x = sp.symbols("x")
    xc = sp.Rational(93, 100)
    N_, M_, L_, bb, hh = 2, 5, 5, 6, 3
    OhV, BoV, dlt, aV = mp.mpf("0.006"), mp.mpf("0.017"), mp.mpf("1e-3"), mp.mpf("1.5")
    betas = [sp.Rational(1, 100), sp.Rational(-1, 200), sp.Rational(1, 500), sp.Rational(-1, 1000)]
    xi = 1 + sum(betas[l - 2] * sp.legendre(l, x) for l in range(2, L_ + 1))
    r2 = sp.expand(xi**2 * (1 - x**2))
    wexpr = -sp.diff(r2, x) / 2
    rexpr = sp.sqrt(r2)
    wf = sp.lambdify(x, wexpr, "mpmath")
    rf = sp.lambdify(x, rexpr, "mpmath")
    ks = [mp.mpf(0)] + [mp.besseljzero(1, j) / bb for j in range(1, M_)]
    psi = lambda xv: (xv - mp.mpf(xc.p) / xc.q) / (1 - mp.mpf(xc.p) / xc.q)
    xcv = mp.mpf(xc.p) / xc.q
    Pt = lambda n, xv: mp.legendre(n, 2 * psi(xv) - 1)
    # slopes from the markdown's (3.33), (3.29), (3.30), (3.31)
    kap_expr = tex2eq(EQ["3.33"].split(r"\qquad")[1])[1]
    fam = {}
    for fname, eq_lab, cname in (("bath", "3.29", "c_m"), ("drop", "3.30", "b_l")):
        l_, r_ = tex2eq(eq(eq_lab))
        fam[fname] = sp.simplify(r_ / S(cname))
    Fcm = sp.diff(tex2eq(eq("3.31"))[1], S("f"))
    gm = tex2eq(inl(r"\gamma_m="))[1]; om2m = tex2eq(inl(r"\omega_m^2="))[1]
    gl = tex2eq(inl(r"\gamma_l="))[1]; om2l = tex2eq(inl(r"\omega_l^2=", mustnot=r"\sigma"))[1]
    sub0 = {S("delta"): sp.Float(str(dlt)), S("frak_a"): sp.Float(str(aV)), S("Oh"): sp.Float(str(OhV)), S("Bo"): sp.Float(str(BoV)),
            S("h_0"): hh}
    def slope(gam, om2, F, extra):
        e = kap_expr.xreplace({S("gamma"): gam, S("omega"): sp.sqrt(om2), S("F"): F}).xreplace(sub0).xreplace(extra)
        return mp.mpf(str(sp.N(e, 30)))
    kappa = [slope(gm, om2m, fam["bath"], {S("k_m"): sp.Float(str(ks[m]), 30)}) if m else mp.mpf(0) for m in range(M_)]
    lam = {l: slope(gl, om2l, fam["drop"], {S("l"): l}) for l in range(2, L_ + 1)}
    kcm = slope(sp.Integer(0), sp.Integer(0), Fcm, {})
    assert all(kappa[m] < 0 for m in range(1, M_)) and kappa[0] == 0 and all(v < 0 for v in lam.values()) and kcm > 0
    J0b = [mp.besselj(0, ks[m] * bb) for m in range(M_)]
    quad = lambda g_: mp.quad(g_, [xcv, 1])
    def A_apply(n):  # returns function x -> (A P~_n)(x)
        p = lambda xv: Pt(n, xv)
        cm = [2 / (bb * J0b[m]) ** 2 * quad(lambda xv: p(xv) * mp.besselj(0, ks[m] * rf(xv)) * wf(xv)) for m in range(M_)]
        bl = {l: quad(lambda xv: p(xv) * xv * mp.legendre(l, xv) * wf(xv)) for l in range(2, L_ + 1)}
        f = 2 * quad(lambda xv: p(xv) * wf(xv))
        def Ap(xv):
            eta = sum(kappa[m] * cm[m] * mp.besselj(0, ks[m] * rf(xv)) for m in range(M_))
            zd = xv * sum(lam[l] * bl[l] * mp.legendre(l, xv) for l in range(2, L_ + 1))
            return eta - kcm * f + zd
        return Ap
    Aps = [A_apply(n) for n in range(N_ + 1)]
    G = mp.matrix(N_ + 1, N_ + 1)
    for n in range(N_ + 1):
        for n2 in range(N_ + 1):
            G[n, n2] = quad(lambda xv: Pt(n, xv) * Aps[n2](xv) * wf(xv))
    scale = max(abs(G[i, j]) for i in range(N_ + 1) for j in range(N_ + 1))
    asym = max(abs(G[i, j] - G[j, i]) for i in range(N_ + 1) for j in range(N_ + 1)) / scale
    ev = mp.eigsy((G + G.T) / 2, eigvals_only=True)
    assert asym < mp.mpf("1e-20"), f"asymmetry {mp.nstr(asym, 3)}"
    assert max(ev) < mp.mpf("1e-20") * scale, f"largest eigenvalue {mp.nstr(max(ev), 3)}"
    return f"N=2, M=L=5, x_c=0.93: relative asymmetry {mp.nstr(asym, 2)}, eigenvalues {[mp.nstr(e, 3) for e in ev]}"


# ----------------------------------------------------------------------------------------------------------
def main():
    for run in CHECKS:
        run()
    width = max(len(r[1]) for r in RESULTS)
    ok = all(r[0] for r in RESULTS)
    print(f"Verification of {MD_PATH.name}   ({len(EQ)} labelled display equations, {len(INL)} inline expressions read)\n")
    for passed, name, kind, where, note in RESULTS:
        print(f"[{'PASS' if passed else 'FAIL'}] [{kind}] {name}\n        binds: {where}\n        {note}\n")
    npass = sum(r[0] for r in RESULTS)
    print(f"{npass}/{len(RESULTS)} checks passed")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
