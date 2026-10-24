# Coupled-run study of the unilateral (Signorini) closure against the other contact-edge rules and the
# bath experiments (AlventosaEtAl2023). Metrics are exactly those of scripts/compare_fluid_experiment.jl:
# t_c and alpha from the centre-of-mass crossings of z = R (threshold_*), delta from max_penetration_depth.
#
#   julia --project=. -t auto scripts/signorini_study.jl [--M 60 --L 60 --nq 40] [--points water:0.725074,oil:1.215765]
#
# Output is a table; nothing is written to data/. Read the caveats in src/signorini.jl before quoting it:
# the pointwise pressure of this closure is not a converged output, only its moments are.
using SpectralKM, DelimitedFiles, Printf, Base.Threads

include(joinpath(@__DIR__, "signorini_study_lib.jl"))

function main()
    M, L, nq = argval("--M", 60), argval("--L", 60), argval("--nq", 40)
    pts = [(split(s, ":")[1], parse(Float64, split(s, ":")[2])) for s in split(argval("--points", "water:0.725074,water:1.938695,oil:1.215765"), ",")]
    configs = [("crossing", :crossing, 0), ("feasible", :feasible, 0), ("signorini n=16", :signorini, 16),
               ("signorini n=24", :signorini, 24), ("signorini n=32", :signorini, 32)]
    @printf("M=%d L=%d nq=%d, threads=%d\n", M, L, nq, nthreads())
    results = Dict()
    jobs = [(pt, cfg) for pt in pts for cfg in configs]
    @threads for j in eachindex(jobs)
        (which, We), (lbl, sel, n) = jobs[j]
        results[j] = try
            run_case(which, We, sel; n=max(n, 2), M=M, L=L, nq=nq)
        catch err
            (tc=NaN, al=NaN, dl=NaN, wall=NaN, diag=NamedTuple[], err=sprint(showerror, err))
        end
    end
    for (ip, (which, We)) in enumerate(pts)
        e = exp_point(which, We)
        @printf("\n%s We=%.4f   experiment: tc=%.3f±%.3f  alpha=%.3f±%.3f  delta=%.3f±%.3f\n", which, We, e.tc, e.tcsd, e.cor, e.corsd, e.dl, e.dlsd)
        @printf("  %-16s %8s %7s | %8s %7s | %8s %7s | %6s\n", "closure", "tc", "n_sd", "alpha", "n_sd", "delta", "n_sd", "wall")
        for (ic, (lbl, _, _)) in enumerate(configs)
            r = results[(ip - 1) * length(configs) + ic]
            haskey(r, :err) && (@printf("  %-16s FAILED: %s\n", lbl, first(r.err, 110)); continue)
            @printf("  %-16s %8.4f %+7.2f | %8.4f %+7.2f | %8.4f %+7.2f | %5.0fs\n", lbl, r.tc, (r.tc - e.tc) / e.tcsd,
                    r.al, (r.al - e.cor) / e.corsd, r.dl, (r.dl - e.dl) / e.dlsd, r.wall)
        end
    end
end
main()
