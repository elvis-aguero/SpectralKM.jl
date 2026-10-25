@testset "newton.jl" begin
    @testset "method of manufactured solutions — linear system" begin
        Atrue = [2.0 0.3; -0.1 1.5]
        btrue = [0.4, -0.2]
        R(X) = Atrue * X - btrue
        result = newton_solve(R, [0.0, 0.0])
        @test result.status == Converged
        @test result.X ≈ Atrue \ btrue atol = 1e-6
    end

    @testset "method of manufactured solutions — mildly nonlinear" begin
        Xtarget = [0.3, -0.5, 0.1]
        R(X) = [X[1]^2 - Xtarget[1]^2, sin(X[2]) - sin(Xtarget[2]), X[3]^3 - Xtarget[3]^3]
        result = newton_solve(R, Xtarget .+ [0.05, -0.03, 0.02])
        @test result.status == Converged
        @test result.X ≈ Xtarget atol = 1e-5
    end

    @testset "near-singular Jacobian: LM damping still returns a finite, usable result" begin
        # A deliberately rank-deficient-near X0 system (mirrors the confirmed severe
        # ill-conditioning risk this solver was built to survive).
        R(X) = [X[1] + X[2] - 1.0, 1e-10 * (X[1] - X[2])]
        result = newton_solve(R, [0.0, 0.0])
        @test all(isfinite, result.X)
        @test result.status in (Converged, Stalled, MaxIterExceeded)
    end

    @testset "never throws on a pathological/non-finite-producing residual" begin
        R(X) = [1.0 / X[1], X[2]]  # blows up at X[1]=0
        result = newton_solve(R, [0.0, 1.0])
        @test result.status == Stalled
    end

    @testset "minnorm_gauss_newton: scale-invariant, minimum-norm, independent of the starting guess" begin
        # a tiny-scale affine system (residual ~1e-9, below newton_solve's absolute tolerance): the guess must
        # NOT be returned unchanged
        A = 1e-9 .* [2.0 0.3; -0.1 1.5]; b = 1e-9 .* [0.4, -0.2]
        R(X) = A * X - b
        X0 = [1e-3, 0.0]
        @test norm(R(X0)) < 1e-7                                   # newton_solve would accept X0 as is
        @test newton_solve(R, X0).X == X0
        res = minnorm_gauss_newton(R, X0)
        @test res.status == Converged
        @test res.X ≈ A \ b rtol = 1e-10
        # rank-deficient: x1 + x2 = 1 with an exactly null direction -> the minimum-norm solution [1/2, 1/2],
        # whatever the guess carries in the null direction
        R2(X) = [X[1] + X[2] - 1.0, 2 * (X[1] + X[2]) - 2.0]
        for X0b in ([0.0, 0.0], [5.0, -3.0], [-7.0, 1.0])
            r2 = minnorm_gauss_newton(R2, X0b)
            @test r2.status == Converged
            @test r2.X ≈ [0.5, 0.5] atol = 1e-12
        end
        # mildly nonlinear: same root as newton_solve
        Xt = [0.3, -0.5, 0.1]
        R3(X) = [X[1]^2 - Xt[1]^2, sin(X[2]) - sin(Xt[2]), X[3]^3 - Xt[3]^3]
        r3 = minnorm_gauss_newton(R3, Xt .+ [0.05, -0.03, 0.02])
        @test r3.status == Converged
        @test r3.X ≈ Xt atol = 1e-9
    end
end
