# Tests for the edge-pinned tangency closure, selector = :tangency (src/tangency.jl). What is tested is the
# construction (the edge row holds after the solve; the root is the admissibility boundary), not agreement with
# experiment, which is a measurement.

@testset "tangency closure" begin
    p = Params(We=1.0958, Bo=0.017, Oh=0.006, b=6.0, h0=3.0, M=40, L=80, N=3, nq=80,
               selector=:tangency, check_budget=false)
    levels, diag, phases = run_simulation(Params(We=1.0958, Bo=0.017, Oh=0.006, b=6.0, h0=3.0, M=40, L=80, N=3, nq=80,
                                                 selector=:crossing, check_budget=false); t_end=0.6)
    ic = findall(==(InContact), phases)
    k = ic[end]                                       # a step well into contact
    hist = SimHistory(levels[k-2], levels[k-1]); dt = levels[k].dt

    @testset "edge-pinned solve: C(theta_c) = 0 to round-off, N moments satisfied" begin
        for th in (0.4, 0.6, 0.8)
            chat, res, am, beta, zcm, f = SpectralKM.inner_solve_pinned(hist, dt, th, zeros(p.N + 1), p)
            @test res.status == Converged
            @test abs(C_at_theta(am, beta, zcm, th, p)) < 1e-12
        end
    end

    @testset "selected edge: tangent (|T| small), admissible, and T changes sign across it" begin
        th_prev = levels[k-1].X[end]
        pay = select_theta_c_tangency(hist, dt, th_prev, zeros(p.N + 1), p)
        @test pay !== nothing
        th = pay.theta_c
        @test max_overlap_beyond(pay.am, pay.beta, pay.zcm, th, p) <= 1e-9
        Tc = tangency_residual(pay.am, pay.beta, pay.zcm, th, p)
        lo = SpectralKM._tangency_candidate(hist, dt, th - 0.02, zeros(p.N + 1), p)
        hi = SpectralKM._tangency_candidate(hist, dt, th + 0.02, zeros(p.N + 1), p)
        @test lo.T > 0 && hi.T < 0                    # overlapping below the edge, separating above
        @test abs(Tc) < 1e-3 * max(abs(lo.T), abs(hi.T))
        @test !lo.ok && hi.ok
    end
end
