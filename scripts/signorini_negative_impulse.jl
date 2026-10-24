# How much of the contact IMPULSE does negative pressure carry under the default :crossing rule?
#
#   julia --project=. scripts/signorini_negative_impulse.jl
#
# Measured: 0.63 %, 0.51 %, 4.40 %, 0.12 % of the positive impulse (water We 1.096 and 0.725; oil We 1.216 and
# 7.307, M = L = 60, nq = 40), although 40-55 % of contact steps have SOME negative region. This bounds how much
# the sign condition p >= 0 can move any observable; it is the reason the unilateral closure agrees with the
# edge-rule closures, and it is NOT evidence that p >= 0 is physically right.
using SpectralKM, Printf
using SpectralKM: geom_at_nodes, contact_quad, pressure_poly_raw
function impulse_split(which, We)
    f = which == "water" ? WATER : OIL_5CST; R = 3.5e-4
    c = conditions(drop=f, R=R, V0=sqrt(We*f.sigma/(f.rho*R)))
    p = Params(c; b=6.0, h0=3.0, M=60, L=60, N=3, nq=40, check_budget=false)
    levels, diag, phases = run_simulation(p; t_end=14.0)
    gn, gw = gauss_legendre_nodes(80)
    pos = 0.0; neg = 0.0; nneg = 0; nc = 0
    for (lv, ph) in zip(levels, phases)
        (ph == InContact && lv.X !== nothing) || continue
        chat = lv.X[1:p.N+1]; xc = cos(lv.X[p.N+2])
        x, wq, P, dP, _ = contact_quad(xc, p)               # same quadrature the closure uses
        xi, r, w = geom_at_nodes(lv.drop.beta, x, P, dP, p.L)
        pv = [pressure_poly_raw(chat, xc, xx) for xx in x]
        fp = 2 * sum(max(pv[i], 0) * w[i] * wq[i] for i in eachindex(x))
        fn = 2 * sum(min(pv[i], 0) * w[i] * wq[i] for i in eachindex(x))
        pos += fp * lv.dt; neg += fn * lv.dt; nc += 1; fn < -1e-9 && (nneg += 1)
    end
    @printf("%-5s We=%.3f: positive impulse %.4f, negative impulse %.4f (%.2f%% of positive); steps with any negative force part %d/%d\n",
            which, We, pos, neg, 100*abs(neg)/pos, nneg, nc)
end
for (w, We) in (("water",1.0958),("water",0.725074),("oil",1.215765),("oil",7.307001)); impulse_split(w, We); end
