# The three ways to read t_c and alpha off a simulated trajectory, on the 17 bath experiments, default configuration:
#   com_instant  centre of mass crossing z = R, instantaneous velocities (Alventosa et al.'s model convention)
#   com_fitted   the experimental procedure (parabola over 40 frames at 10,000 fps) applied to the centre of mass
#   north_fitted the experimental procedure applied to the top of the drop crossing z = 2R (what the experiment tracks),
#                with 40 frames (the paper's minimum) and with 80
#
#   julia --project=. scripts/score_metric_definitions.jl <out.csv>
using SpectralKM, DelimitedFiles, Printf
out = ARGS[1]
DATA = joinpath(@__DIR__, "..", "data", "experiments")
open(out, "w") do io
    println(io, "fluid,We,tc_com_instant,alpha_com_instant,tc_com_fitted,alpha_com_fitted,tc_north_fitted,alpha_north_fitted,tc_north80,alpha_north80,delta,tc_exp,tc_sd,alpha_exp,alpha_sd,delta_exp,delta_sd")
    for fl in (:water, :oil)
        raw = readdlm(joinpath(DATA, "bath_experiment_$(fl).csv"), ','; header=true)
        d, h = raw[1], vec(raw[2]); col(n) = Float64.(d[:, findfirst(==(n), h)])
        f = fl == :water ? WATER : OIL_5CST
        R = 3.5e-4
        frame_dt = 1e-4 / sqrt(f.rho * R^3 / f.sigma)            # 10,000 fps in units of t_sigma
        for (i, We) in enumerate(col("We"))
            p = Params(conditions(drop=f, R=R, V0=sqrt(We * f.sigma / (f.rho * R))); b=6.0, h0=3.0)
            levels, _, phases = run_simulation(p; t_end=20.0)
            ts = [l.t for l in levels]
            a = threshold_contact_time(ts, levels); b = threshold_coefficient_of_restitution(ts, levels)
            cf = fitted_rebound_metrics(ts, levels, phases, p.L; frame_dt=frame_dt, point=:com)
            nf = fitted_rebound_metrics(ts, levels, phases, p.L; frame_dt=frame_dt, point=:north)
            n8 = fitted_rebound_metrics(ts, levels, phases, p.L; frame_dt=frame_dt, point=:north, nframes=80)
            v(x, s) = x === nothing ? NaN : getfield(x, s)
            @printf(io, "%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n", fl, We,
                    something(a, NaN), something(b, NaN), v(cf, :tc), v(cf, :alpha), v(nf, :tc), v(nf, :alpha),
                    v(n8, :tc), v(n8, :alpha),
                    max_penetration_depth(levels, p.L), col("tc_over_tsigma")[i], col("tc_sd")[i], col("cor")[i],
                    col("cor_sd")[i], col("delta_over_R")[i], col("delta_sd")[i])
            flush(io)
        end
    end
end
