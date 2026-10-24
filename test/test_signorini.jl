# Tests for the unilateral (Signorini) contact closure, selector = :signorini (src/signorini.jl).
#
# What is tested is deliberately limited to what the closure can honestly claim: that the discrete
# operator is the one the polynomial route would assemble, that the QP solver solves the stated
# problem, and that the PRESSURE MOMENTS (not the pointwise pressure, which is not a converged
# output) are consistent and converge. Whether the closure agrees with experiment is a measurement,
# not a test: see scripts/signorini_study.jl.

@testset "signorini closure" begin

    @testset "nnqp solves min_{p>=0} 1/2 p'Qp - g'p (KKT, against brute-force active-set enumeration)" begin
        for trial in 1:6
            n = 6
            B = [sin(1.7i + 2.3j + 0.9trial) + 0.3cos(0.7i * j + trial) for i in 1:n+2, j in 1:n]   # deterministic data
            Q = B' * B + 1e-3 * I                                                                   # symmetric positive definite
            g = [cos(1.3i + 2.1trial) for i in 1:n]
            p, _ = nnqp(Q, g; rank_tol=0.0)
            @test all(p .>= 0)
            grad = Q * p .- g
            @test all(grad .>= -1e-9 * maximum(abs.(g)))               # dual feasibility
            @test maximum(abs.(p .* grad)) <= 1e-9 * maximum(abs.(g))  # complementarity
            # brute force: best objective over all 2^n active sets with a nonnegative free-set solution
            best = Inf
            for mask in 0:(2^n - 1)
                F = [i for i in 1:n if (mask >> (i - 1)) & 1 == 1]
                s = zeros(n)
                isempty(F) || (s[F] = Q[F, F] \ g[F])
                all(s .>= -1e-12) || continue
                best = min(best, 0.5 * s' * Q * s - g' * s)
            end
            @test isapprox(0.5 * p' * Q * p - g' * p, best; atol=1e-9)
        end
        # semidefinite Q (rank 2) with a positive-force constraint, insensitive to the truncation tolerance
        v1 = collect(1.0:8.0); v2 = cos.(1.0:8.0)
        Qs = v1 * v1' + v2 * v2'; gs = Qs * abs.(cos.(0.7 .* (1:8)))
        pa, _ = nnqp(Qs, gs; rank_tol=1e-13); pb, _ = nnqp(Qs, gs; rank_tol=1e-9)
        @test norm(Qs * pa - Qs * pb) <= 1e-6 * norm(Qs * pa)          # the gap response is unique
    end

    @testset "assembled operator: symmetric, positive semidefinite, and consistent with the polynomial moments" begin
        p = Params(We=1.0958, Bo=0.017, Oh=0.006, b=6.0, h0=3.0, M=30, L=40, N=3, nq=60,
                   selector=:signorini, signorini_nodes=24, check_budget=false)
        lvl = initial_level(p)
        hist = SimHistory(lvl, lvl)
        dt = 1e-3
        kappa, alpha, lambda, gam, kappa_cm, mu = step_affine(hist, dt, p)
        beta = copy(gam)
        for n in (24, 48)
            grid = signorini_grid(0.9, n, p.L)
            Q, g, H, T, e = signorini_assemble(grid, beta, kappa, alpha, lambda, gam, kappa_cm, mu, p)
            @test norm(Q - Q') <= 1e-12 * norm(Q)
            @test minimum(eigvals(Symmetric(Q))) >= -1e-10 * maximum(eigvals(Symmetric(Q)))
        end
        # moments of a smooth pressure on the whole grid agree with the polynomial route, at O(h^2)
        thmax = 0.9; xc = cos(thmax)
        chat = [1.0, -0.4, 0.2, 0.1]
        errs = Float64[]
        for n in (12, 24, 48)
            grid = signorini_grid(thmax, n, p.L)
            Q, g, H, T, e = signorini_assemble(grid, beta, kappa, alpha, lambda, gam, kappa_cm, mu, p)
            pn = [pressure_poly(chat, xc, cos(thmax * (j / n)^SpectralKM._SIG_GRADE)) for j in 0:n]   # the graded grid's own nodes
            f_hat = 2 * dot(e, pn)
            q = contact_quad(xc, p)
            x, wq, P, dP, _ = q
            xi, r, w = geom_at_nodes(beta, x, P, dP, p.L)
            f_poly = 2 * sum(pressure_poly_raw(chat, xc, x[i]) * w[i] * wq[i] for i in eachindex(x))
            cm_hat = p.bath_norm .* (H * pn)
            cm_poly = c_m_all(chat, xc, x, wq, r, w, p.k, p.bath_norm)
            push!(errs, max(abs(f_hat - f_poly) / abs(f_poly), maximum(abs.(cm_hat[2:6] .- cm_poly[2:6])) / maximum(abs.(cm_poly[2:6]))))
        end
        @test errs[3] < errs[2] < errs[1]
        @test errs[3] < errs[1] / 8                        # second order: two halvings cut the error by ~16
        @test errs[3] < 2e-3
    end

    @testset "no contact: free flight above the bath gives p = 0" begin
        p = Params(We=1.0958, Bo=0.017, Oh=0.006, b=6.0, h0=3.0, M=20, L=20, N=3, nq=40,
                   selector=:signorini, signorini_nodes=24, check_budget=false)
        lvl = Level(BathModeState(p.M), DropModeState(p.L), COMState(2.0, -sqrt(p.We)), 0.0, 0.0, nothing)
        out = select_signorini(SimHistory(lvl, lvl), 1e-3, p)
        @test out.nactive == 0 && out.f == 0.0 && out.theta_c == 0.0
    end

    @testset "only the free wall is supported; the clamped-wall term is not in Q" begin
        p = Params(We=1.0, Bo=0.017, Oh=0.006, b=6.0, h0=3.0, M=10, L=10, N=3, nq=20, wall=:clamped,
                   selector=:signorini, signorini_nodes=12, check_budget=false)
        lvl = initial_level(p)
        @test_throws ArgumentError select_signorini(SimHistory(lvl, lvl), 1e-3, p)
    end

    @testset "coupled impact (M=L=60): momentum identity, convergence in n, agreement with :crossing" begin
        run1(n; sel=:signorini) = begin
            p = Params(We=1.0958, Bo=0.017, Oh=0.006, b=6.0, h0=3.0, M=60, L=60, N=3, nq=40,
                       selector=sel, signorini_nodes=n, check_budget=false)
            levels, diag, phases = run_simulation(p; t_end=8.0)
            (p=p, levels=levels, diag=diag, phases=phases, ts=[l.t for l in levels])
        end
        a = run1(24); b = run1(32)
        # discrete momentum identity dv = 3/2 * impulse - Bo * dt, as for the polynomial closures
        ic = findall(==(InContact), a.phases)
        i0, i1 = first(ic) - 1, last(ic) + 1
        rows = [d for d in a.diag if haskey(d, :f)]
        imp = sum(0.5 * (rows[i].f + rows[i-1].f) * (rows[i].t - rows[i-1].t) for i in 2:length(rows))
        dv = a.levels[i1].com.v - a.levels[i0].com.v
        @test isapprox(dv, 1.5 * imp - a.p.Bo * (a.levels[i1].t - a.levels[i0].t); rtol=0.01, atol=0.005)
        # the pressure is non-negative by construction at every step (the point of the closure)
        @test all(d.nactive >= 1 for d in rows)
        # the observables are converged in n (node count) to well inside the experimental scatter
        tca = threshold_contact_time(a.ts, a.levels); tcb = threshold_contact_time(b.ts, b.levels)
        ala = threshold_coefficient_of_restitution(a.ts, a.levels); alb = threshold_coefficient_of_restitution(b.ts, b.levels)
        @test isapprox(tca, tcb; rtol=5e-3)
        @test isapprox(ala, alb; rtol=1e-2)
        @test isapprox(max_penetration_depth(a.levels, a.p.L), max_penetration_depth(b.levels, b.p.L); rtol=5e-3)
        # and agree with the default edge rule on the observables (they differ in the pressure, not here)
        c = run1(2; sel=:crossing)
        @test isapprox(tca, threshold_contact_time(c.ts, c.levels); rtol=1e-2)
        # the true nonlinear overlap left after the solve is a small fraction of the free overlap
        worst = maximum(d.max_overlap for d in rows)
        @test worst < 1e-2
    end

    @testset "a rejection cascade raises instead of continuing a garbage trajectory" begin
        SpectralKM.SIGNORINI_STATS[:consecutive] = 0
        lvl = initial_level(Params(We=1.0, Bo=0.01, Oh=0.01, b=6.0, h0=3.0, M=4, L=4, N=3, nq=10, check_budget=false))
        hist = SimHistory(lvl, lvl)
        for _ in 1:SpectralKM.SIGNORINI_MAX_CONSECUTIVE_REJECTIONS
            @test SpectralKM._signorini_reject(hist, 1e-3, NamedTuple[]) === nothing
        end
        @test_throws ErrorException SpectralKM._signorini_reject(hist, 1e-3, NamedTuple[])
        SpectralKM.SIGNORINI_STATS[:consecutive] = 0
    end
end
