# Numerical half of the model problem in derivations/verify_edge_pinning.py (static capillary-gravity membrane, l = 1,
# indented by a sphere of radius 1/rho). Compare plain Galerkin + tangency root against edge-pinned Galerkin +
# tangency root, at every N, against the exact contact radius a*.
using SpectralKM, SpecialFunctions, LinearAlgebra, Printf
const S = SpectralKM
gl(n, lo, hi) = (let (x, w) = S.gauss_legendre_nodes(n); ((hi - lo) / 2 .* x .+ (hi + lo) / 2, (hi - lo) / 2 .* w) end)
legP(n, t) = S.legendre_P_table(n, t)            # P_0..P_n at t
zd(r, rho, D) = -D + (1 - sqrt(1 - (rho * r)^2)) / rho
dzd(r, rho) = rho * r / sqrt(1 - (rho * r)^2)
Dstar(a, rho) = a * rho * besselk(0, a) / (sqrt(1 - a^2 * rho^2) * besselk(1, a)) - sqrt(1 - a^2 * rho^2) / rho + 1 / rho

function solve_patch(a, N, rho, D; pinned::Bool, nq=48, ni=48)
    phi(s) = legP(N, 2 * (1 - s^2 / a^2) - 1)[1:N+1]          # edge s=a is psi=0
    ro, wo = gl(nq, 0.0, a)
    eta = zeros(nq, N + 1)                                      # eta_m(r_i)
    for i in 1:nq
        s1, w1 = gl(ni, 0.0, ro[i]); s2, w2 = gl(ni, ro[i], a)
        A1 = sum(w1[j] * besseli(0, s1[j]) * s1[j] .* phi(s1[j]) for j in 1:ni)
        A2 = sum(w2[j] * besselk(0, s2[j]) * s2[j] .* phi(s2[j]) for j in 1:ni)
        eta[i, :] = -besselk(0, ro[i]) .* A1 .- besseli(0, ro[i]) .* A2
    end
    sq, wq = gl(ni, 0.0, a)
    Qm = sum(wq[j] * besseli(0, sq[j]) * sq[j] .* phi(sq[j]) for j in 1:ni)
    G = zeros(N + 1, N + 1); rhs = zeros(N + 1)
    for i in 1:nq
        ph = phi(ro[i])
        G .+= wo[i] * ro[i] .* (ph * eta[i, :]')
        rhs .+= wo[i] * ro[i] .* ph .* zd(ro[i], rho, D)
    end
    if pinned
        G[end, :] = -besselk(0, a) .* Qm; rhs[end] = zd(a, rho, D)
    end
    c = G \ rhs
    Qc = dot(Qm, c)
    E = -besselk(0, a) * Qc - zd(a, rho, D)
    T = besselk(1, a) * Qc - dzd(a, rho)
    Cin = maximum(abs(dot(eta[i, :], c) - zd(ro[i], rho, D)) for i in 1:nq)
    (; E, T, c, cond=cond(G), pedge=sum(c .* legP(N, -1.0)[1:N+1]), Cin)
end

function root_T(N, rho, D, astar; pinned)
    lo, hi = 0.5astar, 1.5astar
    Tlo = solve_patch(lo, N, rho, D; pinned).T; Thi = solve_patch(hi, N, rho, D; pinned).T
    Tlo * Thi < 0 || return nothing
    for _ in 1:60
        mid = (lo + hi) / 2; Tm = solve_patch(mid, N, rho, D; pinned).T
        (Tm * Tlo > 0) ? (lo = mid; Tlo = Tm) : (hi = mid)
    end
    (lo + hi) / 2
end

for (rho, astar) in ((0.3, 0.12), (0.3, 0.03), (1.0, 0.4))
    D = Dstar(astar, rho)
    @printf("\nrho = R/l = %.2f, exact a* = %.4f (a*/R = %.3f), D = %.6e\n", rho, astar, astar * rho, D)
    @printf("%3s | %12s %11s %10s %9s | %12s %11s %9s %9s\n", "N", "plain a_T", "a_T/a*-1", "E(a_T)", "cond", "pinned a_T", "a_T/a*-1", "cond", "maxC_in")
    for N in (0, 1, 2, 3, 4, 6, 8, 10)
        aP = root_T(N, rho, D, astar; pinned=false)
        aQ = N == 0 ? nothing : root_T(N, rho, D, astar; pinned=true)
        sP = aP === nothing ? nothing : solve_patch(aP, N, rho, D; pinned=false)
        sQ = aQ === nothing ? nothing : solve_patch(aQ, N, rho, D; pinned=true)
        @printf("%3d | %12s %11s %10s %9s | %12s %11s %9s %9s\n", N,
                aP === nothing ? "no root" : @sprintf("%.8f", aP), aP === nothing ? "" : @sprintf("%+.2e", aP / astar - 1),
                sP === nothing ? "" : @sprintf("%+.2e", sP.E), sP === nothing ? "" : @sprintf("%.1e", sP.cond),
                aQ === nothing ? "-" : @sprintf("%.8f", aQ), aQ === nothing ? "" : @sprintf("%+.2e", aQ / astar - 1),
                sQ === nothing ? "" : @sprintf("%.1e", sQ.cond), sQ === nothing ? "" : @sprintf("%.1e", sQ.Cin))
    end
end

# Assertions (the table above is the evidence; these make the script fail if it stops holding).
let rho = 1.0, astar = 0.4, D = Dstar(0.4, 1.0)
    for N in (1, 3, 6)
        aQ = root_T(N, rho, D, astar; pinned=true)
        @assert abs(aQ / astar - 1) < 1e-10 "edge-pinned scheme does not return the exact radius at N=$N"
    end
    for N in (0, 1, 2, 3)
        aP = root_T(N, rho, D, astar; pinned=false)
        sP = solve_patch(aP, N, rho, D; pinned=false)
        @assert sP.E > 0 "plain Galerkin tangency root expected to overlap at the edge at N=$N"
    end
    println("assertions passed: edge-pinned exact at N = 1, 3, 6; plain Galerkin root overlaps at the edge at N = 0..3")
end
