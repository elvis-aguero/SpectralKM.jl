# Convergence of the unilateral (Signorini) closure in the truncations (M, L), the node count n and the
# time step, against the default :crossing rule on the same settings. Observables only (t_c, alpha, delta,
# defined as in scripts/compare_fluid_experiment.jl); the pointwise pressure is not a converged output.
#
#   julia --project=. -t 3 scripts/signorini_convergence.jl
using SpectralKM, Printf, Base.Threads
include(joinpath(@__DIR__, "signorini_study_lib.jl"))

function main()
    cases = [("water", 0.725074), ("oil", 1.215765)]
    trunc = [(40, 40, 40), (60, 60, 40), (90, 90, 60), (60, 120, 200)]
    dts = [(1e-3, 0.02), (5e-4, 0.005)]
    jobs = Tuple[]
    for (w, We) in cases
        for (M, L, nq) in trunc, sel in (:crossing, :signorini)
            push!(jobs, (w, We, sel, M, L, nq, 1e-3, 0.02))
        end
        for (di, dm) in dts[2:2], sel in (:crossing, :signorini)
            push!(jobs, (w, We, sel, 60, 60, 40, di, dm))
        end
    end
    res = Vector{Any}(undef, length(jobs))
    @threads for j in eachindex(jobs)
        w, We, sel, M, L, nq, di, dm = jobs[j]
        res[j] = try run_case(w, We, sel; n=24, M=M, L=L, nq=nq, dt_init=di, dt_max=dm)
                 catch err; (tc=NaN, al=NaN, dl=NaN, wall=NaN, err=sprint(showerror, err)) end
    end
    last = ""
    for (j, job) in enumerate(jobs)
        w, We, sel, M, L, nq, di, dm = job
        tag = "$w We=$(round(We, digits=3))"
        tag != last && (@printf("\n%s\n  %-10s %-18s %-9s %8s %8s %8s %6s\n", tag, "closure", "(M,L,nq)", "dt_init/max", "tc", "alpha", "delta", "wall"); last = tag)
        r = res[j]
        haskey(r, :err) ? @printf("  %-10s (%d,%d,%d) dt %.0e/%.3f  FAILED: %s\n", sel, M, L, nq, di, dm, first(r.err, 90)) :
            @printf("  %-10s %-18s %.0e/%.3f %9.4f %8.4f %8.4f %5.0fs\n", sel, "($M,$L,$nq)", di, dm, r.tc, r.al, r.dl, r.wall)
    end
end
main()
