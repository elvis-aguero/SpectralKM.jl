# Summarise scripts/score_vs_experiment.jl outputs: one line per configuration.
#
#   julia --project=. scripts/score_summary.jl run1.csv run2.csv ...
#
# FIXED BEFORE THE RUNS WERE READ:
#   weights   w_i = We_i^(-log10 2), i.e. a point one decade lower in We weighs twice as much; normalised to
#             sum to one over the points of the comparison.
#   primary   weighted standardised MSE, sum_i w_i * mean over the three metrics of ((model - exp)/sd)^2.
#   secondary per-metric weighted R^2 = 1 - sum w (m - e)^2 / sum w (e - ebar_w)^2, and per-metric weighted
#             raw MSE.
#   A configuration with any failed point (or no contact time) is reported but is not eligible as a default.
using DelimitedFiles, Printf

const METRICS = (("tc", "tc_exp", "tc_sd"), ("alpha", "alpha_exp", "alpha_sd"), ("delta", "delta_exp", "delta_sd"))

function summarise(path)
    raw = readdlm(path, ','; header=true)
    d, h = raw[1], vec(raw[2])
    c(n) = d[:, findfirst(==(n), h)]
    We = Float64.(c("We")); status = String.(c("status"))
    cfg = join(string.(d[1, 3:5]), "/")
    nfail = count(!=("ok"), status)
    ok = status .== "ok"
    w = We[ok] .^ (-log10(2)); w ./= sum(w)
    zs = zeros(count(ok)); r2 = Float64[]; mse = Float64[]
    for (m, e, s) in METRICS
        mv, ev, sv = Float64.(c(m))[ok], Float64.(c(e))[ok], Float64.(c(s))[ok]
        zs .+= ((mv .- ev) ./ sv) .^ 2 ./ length(METRICS)
        eb = sum(w .* ev)
        push!(r2, 1 - sum(w .* (mv .- ev) .^ 2) / sum(w .* (ev .- eb) .^ 2))
        push!(mse, sum(w .* (mv .- ev) .^ 2))
    end
    primary = sum(w .* zs)
    @printf("%-30s n=%2d failed=%d | weighted std-MSE %.3f | R2 tc %+.3f alpha %+.3f delta %+.3f | MSE tc %.4f alpha %.2e delta %.4f%s\n",
            cfg, length(We), nfail, primary, r2..., mse..., nfail > 0 ? "  (NOT ELIGIBLE)" : "")
end

foreach(summarise, ARGS)
