# Damped Newton solve for R(X)=0 (design doc §subsec:newton, eq:newton-step), with the
# Jacobian obtained by ForwardDiff.jacobian rather than a hand-derived analytic
# Jacobian (km-viscous-drop's DropSolver approach, deliberately not carried over: see
# design doc's own disclosure that R involves J0, oscillatory, composed with an
# argument affine in X — no global convergence guarantee, so this solver never throws
# and reports enough diagnostics for the caller to judge the outcome, rather than
# asserting success).
#
# EMPIRICAL FINDING beyond what the design doc's own audit anticipated: the Jacobian
# of R is not merely "possibly multi-rooted" (the disclosed J0-oscillation risk) but
# routinely severely ill-conditioned even at modest N — condition numbers ~1e17-1e18
# were observed on ordinary test cases, with singular values spanning ~13 orders of
# magnitude WITHIN the pressure-coefficient (ĉ_n) subspace alone, independent of the
# rescaled ψ basis (pressure.jl) and independent of the timestep. This is a real
# near-rank-deficiency, not just a poor-scaling artifact of the ψ=(x-xc)/(1-xc)
# reparametrization (which fixed a SEPARATE, and now negligible, n-dependent decay).
# A plain `J \ R` Newton step is essentially numerical noise on such a system (the
# achievable relative accuracy of a linear solve is ~cond(J)*eps, i.e. `>> 1` here).
#
# Fix: solve the Newton step by Levenberg-Marquardt-damped least squares,
# `(JᵀJ + μI) ΔX = JᵀR(X)`, with μ adapted each iteration (standard, well-established
# practice for exactly this situation — a Jacobian that is singular or near-singular
# at points of interest — rather than assuming an ordinary Newton step is well-posed).
# This trades quadratic for linear local convergence but is robust to the observed
# near-singularity; the damping μ shrinks toward 0 as the iteration succeeds, so
# well-conditioned regions still converge fast.

@enum NewtonStatus Converged MaxIterExceeded Stalled

"""
    NewtonResult

- `X`: final iterate (converged solution if `status==Converged`, best attempt otherwise).
- `status`: `Converged`, `MaxIterExceeded`, or `Stalled` (LM damping `μ` was driven to
  its ceiling without residual convergence — a stronger warning sign than simply
  running out of iterations).
- `resid_norm_hist`, `step_norm_hist`: per-iteration `‖R(X)‖`, `‖ΔX‖` — kept as part of
  the saved output (the design doc's own disclosed convergence risk, and the
  ill-conditioning found during implementation, both call for this, not just a debug
  log discarded after the run).
"""
struct NewtonResult
    X::Vector{Float64}
    status::NewtonStatus
    resid_norm_hist::Vector{Float64}
    step_norm_hist::Vector{Float64}
    n_iter::Int
end

"""
    newton_solve(R, X0; tol_abs=1e-7, tol_rel=1e-6, maxiter=50, mu0=1e-6, mu_max=1e10) -> NewtonResult

Default tolerances are deliberately looser than the naive `1e-10` one might reach for:
given the confirmed `~1e17-1e18` Jacobian condition numbers (module note), driving the
residual below `~1e-7` is not reliably achievable in double precision regardless of
algorithm, and demanding it produces spurious `Stalled` outcomes on otherwise-good
solutions (most acutely right at contact onset, where the true contact patch has zero
width and grows like `√t` — a genuine degeneracy, not a solver failure).

Levenberg-Marquardt-damped Newton's method on `R(X)=0`, `R: R^n -> R^n`, from warm
start `X0`. Never throws. Convergence criterion: `‖R(X)‖ < max(tol_abs, tol_rel*‖R(X0)‖)`
(mixed absolute/relative, since the warm start can be excellent mid-bounce or poor at
contact onset). At each iterate, solves `(JᵀJ + μI)ΔX = JᵀR(X)` rather than `JΔX=R(X)`
directly — necessary because `J` was empirically found to be severely ill-conditioned
(see module note), not merely occasionally singular. `μ` is decreased (toward a floor)
after a successful (residual-reducing) step and increased (up to `mu_max`) after a
failed one, the standard adaptive LM strategy.
"""
function newton_solve(R, X0::AbstractVector{Float64}; tol_abs::Float64=1e-7,
    tol_rel::Float64=1e-6, maxiter::Int=50, mu0::Float64=1e-6, mu_max::Float64=1e10,
    mu_floor::Float64=1e-12)
    X = copy(X0)
    R0 = R(X)
    resid0 = norm(R0)
    tol = max(tol_abs, tol_rel * resid0)

    resid_hist = Float64[resid0]
    step_hist = Float64[]

    resid_norm = resid0
    mu = mu0
    for iter in 1:maxiter
        if resid_norm < tol
            return NewtonResult(X, Converged, resid_hist, step_hist, iter - 1)
        end
        Rval = R(X)
        if !all(isfinite, Rval)
            return NewtonResult(X, Stalled, resid_hist, step_hist, iter - 1)
        end
        J = ForwardDiff.jacobian(R, X)
        if !all(isfinite, J)
            return NewtonResult(X, Stalled, resid_hist, step_hist, iter - 1)
        end
        JtJ = J' * J
        JtR = J' * Rval
        n = length(X)

        accepted = false
        step_norm = 0.0
        Xtrial = X
        resid_trial = resid_norm
        local_tries = 0
        while !accepted && local_tries < 30
            deltaX = try
                (JtJ + mu * I(n)) \ JtR
            catch e
                e isa Union{LinearAlgebra.SingularException,ArgumentError} || rethrow()
                mu = min(mu * 4, mu_max)
                local_tries += 1
                continue
            end
            if !all(isfinite, deltaX)
                mu = min(mu * 4, mu_max)
                local_tries += 1
                continue
            end
            Xtrial = X .- deltaX
            resid_trial = norm(R(Xtrial))
            if !isfinite(resid_trial)
                mu = min(mu * 4, mu_max)
                local_tries += 1
                continue
            end
            if resid_trial < resid_norm
                accepted = true
                step_norm = norm(deltaX)
                mu = max(mu / 3, mu_floor)
            else
                mu = min(mu * 4, mu_max)
                local_tries += 1
                if mu >= mu_max
                    step_norm = norm(deltaX)
                    break
                end
            end
        end

        push!(step_hist, step_norm)
        X = Xtrial
        resid_norm = resid_trial
        push!(resid_hist, resid_norm)

        if !accepted && mu >= mu_max
            return NewtonResult(X, Stalled, resid_hist, step_hist, iter)
        end
    end

    if resid_norm < tol
        return NewtonResult(X, Converged, resid_hist, step_hist, maxiter)
    end
    return NewtonResult(X, MaxIterExceeded, resid_hist, step_hist, maxiter)
end

"""
    minnorm_gauss_newton(R, X0; rtol_ladder=(1e-10, 1e-8, 1e-6, 1e-4, 1e-3), tol=1e-11, maxiter=30) -> NewtonResult

Gauss-Newton for the inner contact system, with each iterate the MINIMUM-NORM solution of the
linearised problem, `X_{k+1} = J_r^+ (J X_k - R(X_k))`, where `J_r^+` is the pseudo-inverse of
`J = dR/dX` keeping singular values above `rtol_svd * sigma_max`.

Why not `newton_solve`. That routine accepts an iterate once `‖R‖ < 1e-7` in ABSOLUTE terms. The
inner residual is a gap moment integrated over the contact patch, so it scales with the patch
area: on the first contact steps the starting guess (`[1e-3, 0, ...]` at onset) is already below
`1e-7` and was returned unchanged as `Converged`, and on later early steps the
Levenberg-Marquardt damping (also absolute, `mu0 = 1e-6` against `J'J ~ 1e-12`) left the iteration
at its cap with the residual barely reduced. Measured on the water impact We = 1.0958 (M = L = 60,
N = 3): 2 of 217 contact steps returned the guess and 38 hit the iteration cap; at two of those
the accepted net force had the wrong sign relative to the actual solution of the same equations.

Why minimum norm. On small patches `J` is genuinely rank-deficient to working precision (measured
condition numbers 1e8-1e14 for theta_c < 0.05 at M = L = 60): the truncated bath and drop modes
cannot resolve pressure structure finer than their shortest wavelengths, so some combinations of
the coefficients change no pressure moment the dynamics see. Those combinations are not
determined by the equations, and the minimum-norm iterate sets them to zero instead of carrying
whatever the starting guess or round-off put there. Because the inner map is affine up to the
mild geometry coupling, the iterate is independent of `X0` and converges in a few steps.

The truncation level: the most truncated level on `rtol_ladder` that satisfies every equation, else
the most truncated level that converged on its retained subspace. A level that keeps nearly-null
directions asks for coefficients of order 1e6 (measured at theta_c = 0.046) or worse and is never
preferred over one that already solves the system. Above theta_c ~ 0.1 every level gives the
same answer (residual 1e-14 relative, identical f); below it the answer depends on the level,
which is the identifiability limit of the truncated modes, not a solver defect.

Convergence is judged relative to the problem's own scale `s = max(‖R(0)‖, ‖R(X0)‖)` (the gap
moments with no pressure), not in absolute terms: `‖U_r' R‖ <= tol * s` on the retained subspace and
a step below `1e-9` relative. `resid_norm_hist` records the full residual, which keeps a component in
the discarded directions; that component cannot be removed by any bounded pressure.
"""
function minnorm_gauss_newton(R, X0::AbstractVector{Float64};
                              rtol_ladder=(1e-3, 1e-4, 1e-6, 1e-8, 1e-10), full_tol::Float64=1e-10, kwargs...)
    # Most truncated level first. A level that satisfies EVERY equation (full residual <= full_tol * scale) ends the
    # search: keeping more directions cannot improve it and only admits nearly-null ones. If no level satisfies all of
    # them (a patch below what the modes resolve), the most truncated level that converged on its own subspace is
    # returned. Going the other way -- least truncated first -- accepted, at We = 0.0231, theta_c = 0.049, a
    # "converged" iterate with net force 7.8e30.
    # The minimum-norm iterate does not depend on the start except through the mild geometry coupling, so a start
    # that the geometry cannot evaluate (w changes sign, R non-finite) is replaced by zero.
    all(isfinite, R(X0)) || (X0 = zero(X0))
    s = max(norm(R(zero(X0))), norm(R(X0)))
    first_ok = nothing
    local res
    for rt in rtol_ladder
        res = _minnorm_gn(R, X0; rtol_svd=rt, kwargs...)
        res.status == Converged || continue
        first_ok === nothing && (first_ok = res)
        res.resid_norm_hist[end] <= full_tol * s && return res
    end
    return first_ok === nothing ? res : first_ok
end

function _minnorm_gn(R, X0::AbstractVector{Float64}; rtol_svd::Float64=1e-10,
                     tol::Float64=1e-11, maxiter::Int=30)
    X = copy(X0)
    Rv = R(X)
    s = max(norm(R(zero(X))), norm(Rv))
    resid_hist = Float64[norm(Rv)]
    step_hist = Float64[]
    s == 0 && return NewtonResult(zero(X), Converged, resid_hist, step_hist, 0)
    for iter in 1:maxiter
        J = ForwardDiff.jacobian(R, X)
        (all(isfinite, J) && all(isfinite, Rv)) || return NewtonResult(X, Stalled, resid_hist, step_hist, iter - 1)
        F = svd(J)
        keep = F.S .> rtol_svd * F.S[1]
        Ur, Vr, Sr = F.U[:, keep], F.V[:, keep], F.S[keep]
        Xn = Vr * ((Ur' * (J * X .- Rv)) ./ Sr)
        Rn = R(Xn)
        all(isfinite, Rn) || return NewtonResult(X, Stalled, resid_hist, step_hist, iter - 1)
        step = norm(Xn .- X)
        push!(step_hist, step); push!(resid_hist, norm(Rn))
        X, Rv = Xn, Rn
        if norm(Ur' * Rn) <= tol * s && step <= 1e-9 * max(norm(X), 1e-300)
            return NewtonResult(X, Converged, resid_hist, step_hist, iter)
        end
    end
    return NewtonResult(X, MaxIterExceeded, resid_hist, step_hist, maxiter)
end
