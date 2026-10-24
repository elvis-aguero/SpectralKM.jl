# Coupled-run study of the unilateral (Signorini) closure against the other contact-edge rules and the
# bath experiments (AlventosaEtAl2023). Metrics are exactly those of scripts/compare_fluid_experiment.jl:
# t_c and alpha from the centre-of-mass crossings of z = R (threshold_*), delta from max_penetration_depth.
#
#   julia --project=. -t auto scripts/signorini_study.jl [--M 60 --L 60 --nq 40] [--points water:0.725074,oil:1.215765]
#
# Output is a table; nothing is written to data/. Read the caveats in src/signorini.jl before quoting it:
# the pointwise pressure of this closure is not a converged output, only its moments are.
using SpectralKM, DelimitedFiles, Printf, Base.Threads

const DATA = joinpath(@__DIR__, "..", "data", "experiments")
const R_EXP = 3.5e-4
argval(flag, default) = (i = findfirst(==(flag), ARGS); i === nothing ? default :
                         default isa AbstractString ? ARGS[i+1] : parse(typeof(default), ARGS[i+1]))

function exp_point(which, We)
    raw = readdlm(joinpath(DATA, "bath_experiment_$(which).csv"), ','; header=true)
    d, h = raw[1], vec(raw[2]); col(n) = Float64.(d[:, findfirst(==(n), string.(h))])
    i = argmin(abs.(col("We") .- We))
    return (We=col("We")[i], tc=col("tc_over_tsigma")[i], tcsd=col("tc_sd")[i], cor=col("cor")[i], corsd=col("cor_sd")[i],
            dl=col("delta_over_R")[i], dlsd=col("delta_sd")[i])
end

function run_case(which, We, sel; n=48, M=60, L=60, nq=40, t_end=14.0, dt_init=1e-3, dt_max=0.02)
    f = which == "water" ? WATER : OIL_5CST
    c = conditions(drop=f, R=R_EXP, V0=sqrt(We * f.sigma / (f.rho * R_EXP)))
    p = Params(c; b=6.0, h0=3.0, M=M, L=L, N=3, nq=nq, selector=sel, signorini_nodes=n, check_budget=false)
    t0 = time()
    levels, diag, phases = run_simulation(p; t_end=t_end, dt_init=dt_init, dt_max=dt_max)
    ts = [l.t for l in levels]
    tc = threshold_contact_time(ts, levels); al = threshold_coefficient_of_restitution(ts, levels)
    return (tc=tc === nothing ? NaN : tc, al=al === nothing ? NaN : al, dl=max_penetration_depth(levels, p.L),
            wall=time() - t0, diag=diag)
end

function main()
    M, L, nq = argval("--M", 60), argval("--L", 60), argval("--nq", 40)
    pts = [(split(s, ":")[1], parse(Float64, split(s, ":")[2])) for s in split(argval("--points", "water:0.725074,water:1.938695,oil:1.215765"), ",")]
    configs = [("crossing", :crossing, 0), ("feasible", :feasible, 0), ("signorini n=24", :signorini, 24),
               ("signorini n=48", :signorini, 48), ("signorini n=96", :signorini, 96)]
    @printf("M=%d L=%d nq=%d, threads=%d\n", M, L, nq, nthreads())
    results = Dict()
    jobs = [(pt, cfg) for pt in pts for cfg in configs]
    @threads for j in eachindex(jobs)
        (which, We), (lbl, sel, n) = jobs[j]
        results[j] = run_case(which, We, sel; n=max(n, 2), M=M, L=L, nq=nq)
    end
    for (ip, (which, We)) in enumerate(pts)
        e = exp_point(which, We)
        @printf("\n%s We=%.4f   experiment: tc=%.3f±%.3f  alpha=%.3f±%.3f  delta=%.3f±%.3f\n", which, We, e.tc, e.tcsd, e.cor, e.corsd, e.dl, e.dlsd)
        @printf("  %-16s %8s %7s | %8s %7s | %8s %7s | %6s\n", "closure", "tc", "n_sd", "alpha", "n_sd", "delta", "n_sd", "wall")
        for (ic, (lbl, _, _)) in enumerate(configs)
            r = results[(ip - 1) * length(configs) + ic]
            @printf("  %-16s %8.4f %+7.2f | %8.4f %+7.2f | %8.4f %+7.2f | %5.0fs\n", lbl, r.tc, (r.tc - e.tc) / e.tcsd,
                    r.al, (r.al - e.cor) / e.corsd, r.dl, (r.dl - e.dl) / e.dlsd, r.wall)
        end
    end
end
main()
