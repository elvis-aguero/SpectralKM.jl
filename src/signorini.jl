# THE UNILATERAL (SIGNORINI) CONTACT CLOSURE, `selector = :signorini` (FORMULATION.md, track A3).
#
# EXPERIMENTAL. Not the default; every claim below is about what this file does, and each is tied to
# a test (test/test_signorini.jl) or to a measurement in scripts/signorini_study.jl. Nothing here is
# borrowed from the sister repository DropRebound.jl, which solves a different problem (a rigid
# substrate) with a different operator; its conclusions are not assumed to transfer.
#
# THE PROBLEM. With the pressure p(x) >= 0 on 0 <= theta <= theta_max and the contact mismatch
# C = C_free + A p, the unilateral conditions are
#
#       p >= 0,    C <= 0,    p C = 0,
#
# so the contact set {p > 0} is an output of the solve, not an input as in :crossing / :feasible.
# -A is symmetric negative semidefinite in <u,v>_w = int u v w dx (FORMULATION.md, verified by
# derivations/verify_formulation.py), so these are the KKT conditions of the convex programme
#
#       min_{p >= 0}  1/2 <p,Qp> - <p,C_free>,        Q = -A.
#
# DISCRETISATION. Pressure p = sum_i p_i phi_i with piecewise-linear "hat" functions in theta on a
# uniform grid of `p.signorini_nodes + 1` nodes on [0, theta_max]; p >= 0 everywhere iff every p_i >= 0,
# so the cone is represented exactly. Testing against the same hats gives the symmetric matrix
#
#       Q = -( H' diag(kappa .* bath_norm) H  -  2 kappa_cm e e'  +  T' diag(lambda) T ),
#       H[m,i] = int phi_i J0(k_m r) w dx,  T[l,i] = int phi_i x P_l w dx,  e_i = int phi_i w dx,
#
# and g_i = int phi_i C_free w dx. A nodal moment is then c_m = bath_norm_m (H p)_m, b_l = (T p)_l,
# f = 2 e'p, exactly as in the polynomial closures.
#
# WHAT IS AND IS NOT CLAIMED. The truncated operator has small numerical rank (it follows M and L, not
# the node count). Nodes beyond that rank do not add information, and the QP minimiser is then a sparse
# comb of narrow hats: p, its maximum and its support are NOT converged outputs of this closure.
# Only the pressure MOMENTS (c_m, b_l, f), which are all the dynamics see, are meaningful, and they are
# unique where p is not (Q = B' D B with D strictly signed). `theta_c` returned here is therefore the
# moment-weighted extent, never the edge of the support.
#
# GEOMETRY. r(theta) and w(x) depend on beta, which depends on p, so the frozen-geometry QP is wrapped
# in a Picard iteration. It is capped and tolerance-limited because the active set is discontinuous in
# beta: it contracts geometrically while the active set is fixed, and can stall at an active-set noise
# floor otherwise. The number of iterations and the final change are returned.

"""
    SignoriniGrid

Per-step node data for the signorini closure: the theta grid, the quadrature nodes inside every hat
element, and everything that does not depend on the pressure.
"""
struct SignoriniGrid
    n::Int                          # number of elements; n+1 hats
    thmax::Float64
    theta::Vector{Float64}          # quadrature abscissae, length n*ng
    wth::Vector{Float64}            # quadrature weights in theta
    phi_left::Vector{Float64}       # hat left of the abscissa's element (value at the abscissa)
    elem::Vector{Int}               # element index (1..n); left hat is elem, right hat is elem+1
    x::Vector{Float64}              # cos(theta)
    sinth::Vector{Float64}          # sin(theta): dx = sin(theta) dtheta
    P::Vector{Vector{Float64}}      # Legendre tables at x
    dP::Vector{Vector{Float64}}
end

const _SIG_NG = 6                   # Gauss points per hat element
const _SIG_GRADE = 1.6              # grid grading: theta_j = theta_max (j/n)^_SIG_GRADE

"""
    signorini_grid(thmax, n, L; grade) -> SignoriniGrid

Graded grid `theta_j = thmax (j/n)^grade`: the first contact patch is of angular size ~ sqrt(2 v dt) and
is narrower than a uniform element, so the grid is refined toward the pole. Hat functions are piecewise
linear in theta on these (unequal) elements.
"""
function signorini_grid(thmax::Float64, n::Int, L::Int; grade::Float64=_SIG_GRADE)
    gn, gw = gauss_legendre_nodes(_SIG_NG)
    node = [thmax * (j / n)^grade for j in 0:n]
    NQ = n * _SIG_NG
    theta = Vector{Float64}(undef, NQ); wth = similar(theta); phil = similar(theta); elem = Vector{Int}(undef, NQ)
    q = 0
    for e in 1:n
        a, h = node[e], node[e+1] - node[e]
        for j in 1:_SIG_NG
            q += 1
            s = (1 + gn[j]) / 2                 # position within the element in [0,1]
            theta[q] = a + s * h
            wth[q] = gw[j] * h / 2
            phil[q] = 1 - s                     # left hat; right hat is s
            elem[q] = e
        end
    end
    x = cos.(theta)
    P, dP = legendre_tables(x, L)
    return SignoriniGrid(n, thmax, theta, wth, phil, elem, x, sin.(theta), P, dP)
end

"""Largest `theta_max <= asin(0.9998)` such that `w(x) > 0` (droplet not past its widest point) on
`[0, theta_max]` for the geometry `beta`."""
function signorini_theta_max(beta::AbstractVector, p::Params; nsample::Int=600)
    cap = asin(0.9998)
    tm = cap
    for th in range(1e-3, cap; length=nsample)
        if w_of_x(beta, cos(th), p.L) <= 1e-3
            tm = th
            break
        end
    end
    return tm
end

"""
    signorini_assemble(grid, beta, kappa, bath_norm, lambda, kappa_cm, alpha, gam, mu, p)
        -> (Q, g, H, T, e)

Frozen-geometry assembly at the geometry `beta` (the unknown pressure enters only through `Q`).
`gam` is the pressure-free droplet state, so `C_free = sum alpha_m J0(k_m r) - mu + x (1 + sum gam_l P_l)`.
"""
function signorini_assemble(grid::SignoriniGrid, beta::AbstractVector, kappa, alpha, lambda, gam,
                            kappa_cm, mu, p::Params)
    n1 = grid.n + 1
    NQ = length(grid.theta)
    M, L = p.M, p.L
    xi, r, w = geom_at_nodes(beta, grid.x, grid.P, grid.dP, L)
    # measure for dx = sin(theta) dtheta: the integrals are int (...) w sin(theta) dtheta
    meas = [grid.wth[q] * w[q] * grid.sinth[q] for q in 1:NQ]
    H = zeros(M + 1, n1)
    T = zeros(L + 1, n1)
    e = zeros(n1)
    Cfree_w = zeros(n1)
    J0tab = [besselj0(p.k[m+1] * r[q]) for m in 0:M, q in 1:NQ]
    for q in 1:NQ
        i = grid.elem[q]; pl = grid.phi_left[q]; pr = 1 - pl
        mq = meas[q]
        for m in 0:M
            v = J0tab[m+1, q] * mq
            H[m+1, i] += pl * v
            H[m+1, i+1] += pr * v
        end
        xq = grid.x[q]
        for l in 2:L
            v = xq * grid.P[q][l+1] * mq
            T[l+1, i] += pl * v
            T[l+1, i+1] += pr * v
        end
        e[i] += pl * mq; e[i+1] += pr * mq
        # C_free: history-only mismatch at frozen geometry
        eta_free = 0.0
        for m in 0:M
            eta_free += alpha[m+1] * J0tab[m+1, q]
        end
        xif = 1.0
        for l in 2:L
            xif += gam[l+1] * grid.P[q][l+1]
        end
        cf = eta_free - mu + xq * xif
        Cfree_w[i] += pl * cf * mq; Cfree_w[i+1] += pr * cf * mq
    end
    G = H' * Diagonal(kappa .* p.bath_norm) * H
    G .+= T' * Diagonal(lambda) * T
    G .-= 2 * kappa_cm * (e * e')
    Q = -Symmetric((G + G') / 2)
    return Matrix(Q), Cfree_w, H, T, e
end

"""
    nnqp(Q, g; rank_tol) -> (p, nactive)

`min_{p >= 0} 1/2 p'Qp - g'p` for symmetric positive semidefinite `Q`, solved as nonnegative least
squares on the eigen-factorisation `Q = V diag(lam) V'`: with `A = diag(sqrt lam) V'` and
`b = diag(1/sqrt lam) V' g` the objective is `1/2 |A p - b|^2` up to a constant, and Lawson-Hanson
solves that with a pivoted-QR least-squares step on the free set. Eigenvalues below
`rank_tol * max(lam)` are dropped, i.e. directions along which the operator responds with less than
that fraction of its largest response are not allowed to influence the solution.

WHY NOT SOLVE `Q_FF s = g_F` DIRECTLY. `Q` is numerically rank-deficient (its numerical rank follows the
truncations M and L, not the node count), and `g` has components along the near-null directions. A
direct solve with a tiny ridge amplifies those into huge alternating-sign pressures and, through them,
into wrong forces (observed: f ~ 10 against a physical ~ 0.5 at M=30, L=40, n=96). The truncated
factorisation is the minimum-norm minimiser on the resolvable subspace, which is all the moments see.
"""
function nnqp(Q::AbstractMatrix, g::AbstractVector; rank_tol::Float64=1e-12, maxiter::Int=0, _equilibrated::Bool=false)
    n = length(g)
    maxiter = maxiter == 0 ? 6n + 20 : maxiter
    # Jacobi equilibration: with y_i = sqrt(Q_ii) p_i the diagonal is 1, so a narrow hat is not discarded
    # merely for being narrow (its raw energy scales with its width). Nonnegativity of p is unchanged.
    d = _equilibrated ? ones(n) : [Q[i, i] > 0 ? sqrt(Q[i, i]) : 1.0 for i in 1:n]
    if !_equilibrated && !all(==(1.0), d)
        y, na = nnqp(Matrix(Q) ./ (d * d'), g ./ d; rank_tol=rank_tol, maxiter=maxiter, _equilibrated=true)
        return y ./ d, na
    end
    E = eigen(Symmetric(Matrix(Q)))
    lmax = maximum(E.values)
    keep = findall(>(rank_tol * lmax), E.values)
    isempty(keep) && return zeros(n), 0
    A = Diagonal(sqrt.(E.values[keep])) * E.vectors[:, keep]'          # r x n
    b = (E.vectors[:, keep]' * g) ./ sqrt.(E.values[keep])
    pv = zeros(n)
    free = falses(n)
    tol = 1e-12 * max(norm(b), floatmin())
    for _ in 1:maxiter
        wv = A' * (b .- A * pv)                                        # = g - Q p on the resolvable subspace
        cand = [i for i in 1:n if !free[i] && wv[i] > tol]
        isempty(cand) && break
        free[argmax([i in cand ? wv[i] : -Inf for i in 1:n])] = true
        for _ in 1:(2n + 5)
            F = findall(free)
            s = zeros(n)
            s[F] = A[:, F] \ b                                         # min-norm least squares on the free set
            if all(s[F] .> 0)
                pv = s
                break
            end
            neg = [i for i in F if s[i] <= 0]
            alpha = minimum(pv[i] / (pv[i] - s[i]) for i in neg if pv[i] - s[i] > 0; init=1.0)
            pv = pv .+ alpha .* (s .- pv)
            for i in F
                pv[i] <= 1e-14 * max(1.0, maximum(abs.(pv))) && (free[i] = false; pv[i] = 0.0)
            end
            any(free) || break
        end
    end
    pv = max.(pv, 0.0)
    return pv, count(>(0), pv)
end

"""Counters of `select_signorini` calls and rejections, for diagnosing step-size cascades."""
const SIGNORINI_STATS = Dict(:calls => 0, :rej_picard => 0, :rej_w => 0)

"""
    select_signorini(hist, dt, p; picard_tol, picard_maxit) -> payload | nothing

One unilateral contact step. Returns the same payload fields as the other selectors
(`chat` holds the NODAL pressure here, followed by nothing; `theta_c` is the moment-weighted extent),
plus `info` with the solver audit. Returns the zero-pressure ("contact lost") payload with
`theta_c = 0` when the solve finds `p == 0`.
"""
function select_signorini(hist::SimHistory, dt::Float64, p::Params;
                          picard_tol::Float64=1e-7, picard_maxit::Int=1, accept_tol::Float64=1e-2)
    SIGNORINI_STATS[:calls] += 1
    p.wall === :free || throw(ArgumentError("selector=:signorini supports wall=:free only (the clamped-wall rank-one term is not in Q)"))
    kappa, alpha, lambda, gam, kappa_cm, mu = step_affine(hist, dt, p)
    # Geometry is frozen at the second-order extrapolation of the two previous time levels (the order of
    # BDF2), so each step is ONE convex programme. Iterating the geometry to a fixed point (picard_maxit > 1)
    # is available but is not the default: the active set is discontinuous in beta, and the iteration was
    # observed to chatter (and, with rejection on non-convergence, to cascade into dt -> 0) at some (M, L, n).
    dtprev = hist.curr.dt
    ratio = dtprev > 0 ? dt / dtprev : 0.0
    beta = hist.curr.drop.beta .+ ratio .* (hist.curr.drop.beta .- hist.prev.drop.beta)
    thmax = signorini_theta_max(beta, p)
    n = p.signorini_nodes
    grid = signorini_grid(thmax, n, p.L)
    pn = zeros(n + 1)
    dbeta = Inf; nact = 0; it = 0
    cm = zeros(p.M + 1); bl = zeros(p.L + 1); f = 0.0
    for outer in 1:picard_maxit
        it = outer
        Q, g, H, T, e = signorini_assemble(grid, beta, kappa, alpha, lambda, gam, kappa_cm, mu, p)
        pn, nact = nnqp(Q, g)
        cm = p.bath_norm .* (H * pn)
        bl = T * pn
        f = 2 * dot(e, pn)
        beta_new = gam .+ lambda .* bl
        dbeta = maximum(abs.(beta_new .- beta))
        beta = beta_new
        (nact == 0 || dbeta < picard_tol) && break
    end
    # A step whose geometry has not settled, or whose geometry has left the region where the area weight
    # is positive, is REJECTED (the caller halves dt), not accepted with a stale state: the active set is
    # discontinuous in beta, so an unconverged Picard iterate is not "nearly" a solution.
    if nact > 0 && dbeta > accept_tol
        SIGNORINI_STATS[:rej_picard] += 1
        return nothing
    end
    if nact == 0 || f <= 0
        am = copy(alpha); zcm = mu
        return (; chat=pn, result=nothing, am, beta=copy(gam), zcm, f=0.0, cm=zeros(p.M + 1), bl=zeros(p.L + 1),
                theta_c=0.0, picard_its=it, dbeta=dbeta, nactive=0, thmax=thmax, max_overlap=0.0)
    end
    am = alpha .+ kappa .* cm
    zcm = mu + kappa_cm * f
    # moment-weighted extent: theta at which the pressure-weighted mean of r^2 is attained
    xi, r, w = geom_at_nodes(beta, grid.x, grid.P, grid.dP, p.L)
    if !(minimum(w) > 0)
        SIGNORINI_STATS[:rej_w] += 1
        return nothing
    end
    # pressure at quadrature points from the nodal hats
    pq = [grid.phi_left[q] * pn[grid.elem[q]] + (1 - grid.phi_left[q]) * pn[grid.elem[q]+1] for q in eachindex(grid.theta)]
    mw = [grid.wth[q] * w[q] * grid.sinth[q] * pq[q] for q in eachindex(grid.theta)]
    rc2 = sum(mw .* r .^ 2) / sum(mw)
    theta_c = asin(clamp(sqrt(rc2), 1e-6, 0.9998))
    # audit: the true (nonlinear) gap on the final geometry, relative to the free overlap
    ths = range(1e-3, thmax; length=300)
    Cfinal = [C_at_theta(am, beta, zcm, th, p) for th in ths]
    Cfree0 = [C_at_theta(alpha, gam, mu, th, p) for th in ths]
    max_overlap = max(maximum(Cfinal), 0.0)
    return (; chat=pn, result=nothing, am, beta, zcm, f, cm, bl, theta_c, picard_its=it, dbeta=dbeta,
            nactive=nact, thmax=thmax, max_overlap=max_overlap, free_overlap=max(maximum(Cfree0), 0.0))
end
