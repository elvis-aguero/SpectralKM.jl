# Score one model configuration against all 17 bath experiments (5 water, 12 oil; Alventosa et al. 2023,
# data/experiments/bath_experiment_*.csv) on contact time, coefficient of restitution and maximum penetration.
#
#   julia --project=. scripts/score_vs_experiment.jl <selector> <solver> <contact_end> <out.csv> [water|oil|both]
#
# Writes one row per experimental point. Metrics as the experiment reads them (fitted_rebound_metrics: top of the drop
# at z = 2R, 40-frame parabolic fit at 10,000 fps) and max_penetration_depth; default numerical settings, b = 6, h0 = 3,
# t_end = 14. A run that stops (no admissible contact edge while the surfaces overlap) is recorded as FAILED.
#
# The comparison metric is fixed BEFORE any run (scripts/score_summary.jl): the pooled standardised mean squared
# error, mean over points and the three metrics of ((model - experiment)/sd)^2, with per-metric R^2 and raw MSE as
# secondary measures. A configuration with any failed point is not eligible as a default.
using SpectralKM, DelimitedFiles, Printf

sel, solver, cend, out = Symbol(ARGS[1]), Symbol(ARGS[2]), Symbol(ARGS[3]), ARGS[4]
which = length(ARGS) >= 5 ? ARGS[5] : "both"
DATA = joinpath(@__DIR__, "..", "data", "experiments")
fluids = which == "both" ? (:water, :oil) : (Symbol(which),)
open(out, "w") do io
    println(io, "fluid,We,selector,solver,contact_end,status,tc,alpha,delta,tc_exp,tc_sd,alpha_exp,alpha_sd,delta_exp,delta_sd,seconds")
    for fl in fluids
        raw = readdlm(joinpath(DATA, "bath_experiment_$(fl).csv"), ','; header=true)
        d, h = raw[1], vec(raw[2]); col(n) = Float64.(d[:, findfirst(==(n), h)])
        We, tce, tcs, ae, as, de, ds = col("We"), col("tc_over_tsigma"), col("tc_sd"), col("cor"), col("cor_sd"),
                                       col("delta_over_R"), col("delta_sd")
        f = fl == :water ? WATER : OIL_5CST
        for i in eachindex(We)
            R = 3.5e-4
            c = conditions(drop=f, R=R, V0=sqrt(We[i] * f.sigma / (f.rho * R)))
            p = Params(c; b=6.0, h0=3.0, selector=sel, solver=solver, contact_end=cend)
            t0 = time()
            status, tc, al, de_ = "ok", NaN, NaN, NaN
            try
                levels, _, phases = run_simulation(p; t_end=14.0)
                ts = [l.t for l in levels]
                # t_c and alpha as the experiment reads them: top of the drop crossing z = 2R, parabola over 40 frames at
                # 10,000 fps (Alventosa et al. 2023, section 2.2; FORMULATION.md section 3.7)
                fit = fitted_rebound_metrics(ts, levels, phases, p.L; frame_dt=1e-4 / t_sigma(f, R), nframes=40, point=:north)
                a = fit === nothing ? nothing : fit.tc; b_ = fit === nothing ? nothing : fit.alpha
                m = max_penetration_depth(levels, p.L)
                tc = a === nothing ? NaN : a; al = b_ === nothing ? NaN : b_; de_ = m === nothing ? NaN : m
                (isnan(tc) || isnan(al)) && (status = "no_contact_time")
            catch e
                status = "FAILED"
                @warn "run stopped" fl We[i] exception = e
            end
            @printf(io, "%s,%.6f,%s,%s,%s,%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.1f\n", fl, We[i], sel, solver, cend,
                    status, tc, al, de_, tce[i], tcs[i], ae[i], as[i], de[i], ds[i], time() - t0)
            flush(io)
        end
    end
end
