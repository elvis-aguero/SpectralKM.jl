using SpectralKM, DelimitedFiles, Printf

# Shared by signorini_study.jl and signorini_convergence.jl.
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

