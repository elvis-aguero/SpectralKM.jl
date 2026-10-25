# THE EDGE-PINNED TANGENCY CLOSURE, `selector = :tangency` (FORMULATION.md, track A4). EXPERIMENTAL.
#
# THE TWO EDGE CONDITIONS. The kinematic match of Agüero et al. (2026) imposes, at the contact edge,
# both that the surfaces meet, C(theta_c) = 0, and that they meet tangentially, dC/dtheta(theta_c) = 0,
# together with no overlap beyond the edge. They impose the match pointwise at every node of a grid, so
# C(theta_c) = 0 holds by construction and tangency is the one condition left to choose the contact
# radius (by the smallest tangency error among non-overlapping candidates).
#
# The Galerkin match of this package (residual.jl) imposes N+1 MOMENTS of C on the patch, not C itself,
# so C(theta_c) is not zero in general. Measured on the water impact We = 1.0958, the root of
# dC/dtheta then has C(theta_c) > 0, i.e. it overlaps at the edge and is inadmissible
# (+1.2e-3, +2e-4, +1.5e-4 at three steps), and the rule "smallest tangency error among admissible
# candidates" collapses onto the smallest admissible angle, which is selector = :feasible.
#
# THE FIX. Replace the highest Galerkin moment by the edge condition itself:
#
#     R_n = int_{x_c}^1 C P~_n w dx = 0,  n = 0..N-1,       R_N = C(theta_c) * int_{x_c}^1 w dx = 0,
#
# (the edge row is scaled by the patch measure so that every row carries the same units), and choose
# theta_c by dC/dtheta(theta_c) = 0 among candidates with no overlap beyond the edge.
#
# WHY THIS IS A SENSIBLE TRY, NOT WISHFUL THINKING. For the static axisymmetric capillary-gravity
# membrane indented by a sphere, the exterior solution is a one-parameter family -Q K0(r/l), so the edge
# value and slope depend on the pressure ONLY through Q = int I0(s/l) p s ds. Eliminating Q, the two
# edge conditions fix the contact radius by themselves, whatever the interior pressure:
# derivations/verify_edge_pinning.py proves this symbolically, and derivations/verify_edge_pinning_model.jl checks numerically that the
# edge-pinned scheme returns the exact radius at every N (error ~1e-14), whereas the plain Galerkin
# scheme's tangency root overlaps at the edge at every finite N and converges only as N grows. The
# dynamic problem is not one-parameter outside the patch, so exactness does not transfer; the structure
# (two edge conditions, one unknown radius, both imposed) does.

"""
    residual_pinned(chat, xc, q, kappa, alpha, lambda, gam, kappa_cm, mu, p) -> Vector

The edge-pinned inner system (see the file header): Galerkin moments `n = 0..N-1` and the edge row.
"""
function residual_pinned(chat::AbstractVector, xc, q, kappa::Vector{Float64},
    alpha::Vector{Float64}, lambda::Vector{Float64}, gam::Vector{Float64},
    kappa_cm::Float64, mu::Float64, p::Params)
    x, wq, P, dP, Ptil = q
    am, beta, zcm, _, _, _ = unpack_state(chat, xc, q, kappa, alpha, lambda, gam, kappa_cm, mu, p)
    xi, r, w = geom_at_nodes(beta, x, P, dP, p.L)
    T = promote_type(eltype(chat), typeof(xc))
    R = zeros(T, p.N + 1)
    meas = zero(T)
    @inbounds for i in eachindex(x)
        eta = zero(T)
        for m in 0:p.M
            eta += am[m+1] * besselj0(p.k[m+1] * r[i])
        end
        gapterm = (eta - zcm + xi[i] * x[i]) * w[i] * wq[i]
        for n in 0:p.N-1
            R[n+1] += gapterm * Ptil[i][n+1]
        end
        meas += w[i] * wq[i]
    end
    R[p.N+1] = C_at_theta(am, beta, zcm, acos(xc), p) * meas
    return R
end

"""Inner solve at fixed `theta_c` for the edge-pinned system; same return shape as `inner_solve`."""
function inner_solve_pinned(hist::SimHistory, dt::Float64, theta_c::Float64,
    chat_guess::Vector{Float64}, p::Params)
    kappa, alpha, lambda, gam, kappa_cm, mu = step_affine(hist, dt, p)
    xc = cos(theta_c)
    q = contact_quad(xc, p)
    R = chat -> residual_pinned(chat, xc, q, kappa, alpha, lambda, gam, kappa_cm, mu, p)
    result = minnorm_gauss_newton(R, chat_guess)
    am, beta, zcm, f, cm, bl = unpack_state(result.X, xc, q, kappa, alpha, lambda, gam, kappa_cm, mu, p)
    return result.X, result, am, beta, zcm, f, cm, bl
end

"""
    max_overlap_beyond(am, beta, zcm, theta_c, p; nsample=400) -> Float64

Largest `C` on `(theta_c, pi/2]`, sampled. Positive means the surfaces overlap outside the patch.
"""
function max_overlap_beyond(am, beta, zcm, theta_c, p::Params; nsample::Int=400)
    m = -Inf
    for th in range(theta_c, π / 2; length=nsample + 1)[2:end]
        m = max(m, C_at_theta(am, beta, zcm, th, p))
    end
    return m
end

# One candidate: the solved state at theta, its tangency residual and its admissibility. A candidate is
# admissible when the solve is usable, r is monotone on the patch, and nothing overlaps beyond the
# edge by more than `ovl_tol` (the edge itself has C = 0 by construction, so a strict test would reject
# every candidate on round-off).
function _tangency_candidate(hist, dt, theta, guess, p; ovl_tol=1e-9)
    chat, result, am, beta, zcm, f, cm, bl = inner_solve_pinned(hist, dt, theta, guess, p)
    usable = result.status == Converged
    T = usable ? tangency_residual(am, beta, zcm, theta, p) : NaN
    ok = usable && check_monotone_r(beta, cos(theta), p.L) &&
         max_overlap_beyond(am, beta, zcm, theta, p) <= ovl_tol
    return (; theta, T, ok, payload=(; chat, result, am, beta, zcm, f, cm, bl, theta_c=theta))
end

"""
    select_theta_c_tangency(hist, dt, theta_c_prev, chat_guess, p) -> payload | nothing

`theta_c` is the root of the tangency residual `T = dC/dtheta(theta_c)` of the edge-pinned solution.
With `C(theta_c) = 0` imposed, `T > 0` means the bath rises above the drop just outside the edge (overlap,
inadmissible) and `T < 0` that it falls away, so the root is the lower boundary of the admissible set.

Search: start at `theta_c_prev`; if that candidate is inadmissible, walk UP geometrically until one is
admissible; if it is admissible, walk DOWN until one is not (or the floor is reached). The two ends
bracket the boundary, which is refined by regula falsi (Illinois) on `T` when `T` changes sign across
the bracket and by bisection on admissibility otherwise. The admissible end is returned. If the floor is
reached while still admissible, the floor candidate is returned: the patch is then below what the modes
resolve (see FORMULATION.md). Agüero et al.'s fallback, the smallest `|T|`, is not used: `T` vanishes
at `theta = 0` for every state by axisymmetry, so `|T|` is biased towards the smallest patch.
"""
function select_theta_c_tangency(hist::SimHistory, dt::Float64, theta_c_prev::Float64,
                                 chat_guess::Vector{Float64}, p::Params; onset::Bool=false,
                                 theta_tol::Float64=1e-6, maxwalk::Int=24, maxit::Int=40)
    lo_cap, hi_cap = asin(THETA_C_FLOOR_R), THETA_C_SAT
    cand(th) = _tangency_candidate(hist, dt, th, chat_guess, p)
    th0 = clamp(theta_c_prev, lo_cap, hi_cap)
    c0 = cand(th0)
    bad, good = nothing, nothing
    if c0.ok
        good = c0
        th = th0
        for _ in 1:maxwalk
            th <= lo_cap && return good.payload           # admissible down to the floor
            th = max(th / 1.25, lo_cap)
            c = cand(th)
            if c.ok
                good = c
            else
                bad = c; break
            end
        end
        bad === nothing && return good.payload
    else
        bad = c0
        th = th0
        for _ in 1:maxwalk
            th >= hi_cap && return nothing                # nothing admissible below saturation
            th = min(th * 1.25, hi_cap)
            c = cand(th)
            if c.ok
                good = c; break
            else
                bad = c
            end
        end
        good === nothing && return nothing
    end
    # Refine. Classify by the sign of T (T > 0: overlap just beyond the edge), not by admissibility alone:
    # admissibility uses a small overlap tolerance, so judging by it would stop up to ~1e-3 rad short of the
    # root on the overlapping side. A point with T <= 0 that is still inadmissible (an overlap far from the
    # edge, or an unusable solve) also counts as the bad side.
    a, b = bad, good
    side = 0
    Ta, Tb = a.T, b.T
    for _ in 1:maxit
        abs(b.theta - a.theta) < theta_tol && break
        mid = if isfinite(Ta) && isfinite(Tb) && Ta * Tb < 0
            (a.theta * Tb - b.theta * Ta) / (Tb - Ta)
        else
            (a.theta + b.theta) / 2
        end
        lo_, hi_ = minmax(a.theta, b.theta)
        (isfinite(mid) && lo_ < mid < hi_) || (mid = (a.theta + b.theta) / 2)
        c = cand(mid)
        if c.ok && isfinite(c.T) && c.T <= 0
            b = c; Tb = c.T; side == 1 && (Ta /= 2); side = 1
        else
            a = c; Ta = c.T; side == -1 && (Tb /= 2); side = -1
        end
    end
    return (a.ok && isfinite(a.T) && abs(a.T) < abs(b.T)) ? a.payload : b.payload
end
