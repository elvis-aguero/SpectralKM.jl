# A fully spectral kinematic match for droplet–bath impact

[![CI](https://github.com/elvis-aguero/SpectralKM.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/elvis-aguero/SpectralKM.jl/actions/workflows/ci.yml)

Numerical companion to a model of non-coalescing droplet impact on a liquid bath of finite size in which **both** interfaces are spectral and the contact pressure is solved for rather than assumed: a Fourier–Bessel bath, a Legendre droplet, and a shifted-Legendre pressure supported exactly on the contact patch. The model formulation is in [`FORMULATION.md`](FORMULATION.md).

![Impact simulation](assets/impact_demo.gif)

*Water, `We = 1.0958`, `Bo = 0.017`, `Oh = 0.006`, `M = L = 120`, `N = 6`. The inset is the solved pressure over the contact patch (not a converged profile; see below).*

It sits between two published models. From Alventosa, Cimpeanu & Harris (2023, *JFM* **958**, A24; [arXiv:2209.13276](https://arxiv.org/abs/2209.13276)), the "1PKM", it takes the spectral bath and droplet and their modal oscillator equations. From Agüero et al. (2026), the "full KM", it takes contact enforced over the whole pressed region with the pressure as an unknown. New here: the pressure is spectral rather than mesh-resolved, and the contact angle is an outer scalar parameter rather than an unknown in a joint system.

| Path | Contents |
|---|---|
| `src/`, `test/` | The `SpectralKM` package: model, time stepper, tests. |
| `FORMULATION.md` | Model formulation (JFM style, follows Alventosa et al. §3). |
| `scripts/` | Validation, rendering and sweep scripts (`sweep.jl` is parallel and resumable). |
| `data/experiments/` | Measured bath-rebound metrics (`bath_experiment_{water,oil}.csv`, the validation target) and the published 1PKM/DNS curves (`bath_model_curves_*`, model output, not measurement). See `data/experiments/PROVENANCE.md` for what is and is not verified. |
| `notebooks/tutorial.ipynb` | Guided tour, dependency-free (SVG figures). |
| `derivations/` | Full derivation (`paper-formulation.tex`), maintainer notes (`provenance.tex`), the two parent papers, and the CAS / numerical audit scripts behind every quantitative claim. Scripts marked `SUPERSEDED` verify withdrawn claims and are kept for the record. |

## Setup and a first run

Julia 1.12+; `ffmpeg` only for the renderer.

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

```julia
using SpectralKM
p = Params(We = 1.0958, Bo = 0.017, Oh = 0.006, b = 6.0, h0 = 3.0)   # M, L, N, nq default to 60, 120, 3, 200
levels, diag, phases = run_simulation(p; t_end = 8.0, dt_init = 1e-3)
```

`levels` holds, per accepted step, the bath modes `a_m`, droplet modes `β_l`, centre-of-mass height and velocity, and the solved `(ĉ_0…ĉ_N, θ_c)`. `diag` holds contact angle, net force, tangency residual and Newton status.

**From a fluid and a drop instead of dimensionless groups:**

```julia
c = conditions(drop = :oil, R = 3.5e-4, V0 = 0.2724)    # SI units; also bath_radius=, bath_depth=
p = Params(c; b = 6.0, h0 = 3.0)
```

Presets `:water` and `:oil` use Alventosa et al.'s table 1 verbatim (they reproduce the `(Bo, Oh)` printed on their figures); any other liquid is `Fluid(name, ρ, σ, ν)`. Drop and bath are the same liquid, as the model requires: a genuinely two-fluid specification throws an `ArgumentError`.

## Options that matter

| Option | Default | Notes |
|---|---|---|
| `M, L, N, nq` | 60, 120, 3, 200 | Judge convergence on **contact time**, not CoR (CoR settles first). At `We = 1.0958`, `M = L = 30` is within 0.55 % of the settled contact time and `60 → 120` changes it by < 1e-5 relative. **`N = 3` is cheap, not converged.** Raising `N` alone stops helping once it exceeds the operator's numerical rank, which grows with `M` and `L`; raised together (`N` to 32, `M` to 640) the pressure converges except at the edge of large patches, and at fixed `θ_c` the default is 0.9–2.8 % off in the force (`FORMULATION.md` §3.6). Do not quote pointwise pressure at the defaults as converged. |
| `viscous` | `:reid` | `:reid`: exact roots of Reid's (1960) equation, valid at any `Oh`. `:lamb`: the published small-`Oh` closure (overdamps by 4 % at `l = 2`, 23 % at `l = 120` even at `Oh = 0.006`). Effect on contact time is small (0.24 % in oil). |
| `selector` | `:crossing` | `θ_c` is the largest crossing of the two surfaces, previous value held if none (the 1PKM reference rule). `:feasible` (`θ_c = inf` of the admissible set) is degenerate and at high `We` loses contact (oil, `We = 7.307`: contact 0.19 against 5.26 measured); the run now stops with an error there instead of integrating the drop through the bath. Where both answer they agree to 0.03 % in contact time. See `derivations/tangency-selector.tex`. |
| `wall` | `:free` | `:free`: 90° line, no flux (both parent papers). `:clamped`: pinned line by a Lagrange multiplier on the same basis; conserves volume. `:pinned`: Dirichlet basis, **does not conserve volume** (17.5 R³ created from nothing); diagnostic only. The wall condition is negligible at `b = 6` (ΔCoR = −0.001) and matters for smaller baths. |

**Experimental: `selector = :signorini`** solves the unilateral contact problem (`p ≥ 0`, gap `≤ 0`, `p·gap = 0`) with the contact set as an output. Inside its verified envelope it matches `:feasible`/`:crossing` to ~1e-3 on contact time, `α` and `δ`, and it runs where `:feasible` fails (oil, `We = 7.307`). It adds **no accuracy** over `:crossing`: negative pressure carries only 0.1–4.4 % of the contact impulse (`scripts/signorini_negative_impulse.jl`), so the sign condition barely moves the observables. Outside the envelope it can return wrong numbers **without an error**, where `:crossing` degrades smoothly: use `M = L ≥ 50`, `signorini_nodes` 8–40 (default 24), `wall = :free`, `dt_init` in `[6e-4, 3e-3]`, and always cross-check against `:crossing` at the same settings. Its pointwise pressure and support are not converged outputs; only the moments are. See `FORMULATION.md` §3.7 (A3), `scripts/signorini_study.jl` and `scripts/signorini_convergence.jl`.

**Experimental: `selector = :tangency`** imposes the two edge conditions of Agüero et al. (2026), surfaces meeting (`C(θ_c) = 0`, replacing the top Galerkin moment) and meeting tangentially (`∂C/∂θ = 0`), with no overlap beyond the edge (`derivations/verify_edge_pinning.py` and `_model.jl` show it is exact for a static membrane). Mid-contact it is clean (positive force, tangency to 1e-9) and its observables are within 0.2 % (`t_c`), 2.2 % (`α`) and 0.1 % (`δ`) of `:crossing` on water `We = 1.0958`. It **fails** on patches below the modal resolution (onset, after separation) and stops early on oil `We = 7.307`. Not for production.

**Known solver defect.** The default inner solver accepts an absolute residual `< 1e-7`, so on the first and last contact steps it can return its starting guess or an unconverged iterate (wrong-sign force on 6 of 195 steps at water `We = 1.0958`; mid-contact it agrees with the true solution to < 1 %). `minnorm_gauss_newton` fixes it but changes `:crossing`'s end of contact; switching is pending a decision (`FORMULATION.md` §3.6).

**Which contact time?** A run can detach and re-attach a few times. `primary_contact_time` (the longest interval, normally the first) is the metric; the sum and the first-to-last span differ by up to 17 %. If `contact_intervals` shows the longest interval is not the first, the truncation is too coarse to trust.

## Validation

```bash
julia --project=. scripts/compare_fluid_experiment.jl water    # or oil
julia --project=. scripts/validate_experimental.jl             # trajectories vs measured top and bottom
```

Mean |relative error| over all 17 measured points (12 oil, 5 water), with the published models on the same points:

| | Oil: SpectralKM | 1PKM | DNS | Water: SpectralKM | 1PKM | DNS |
|---|---|---|---|---|---|---|
| contact time | **2.7 %** | 11.0 % | 3.3 % | 11.5 % | 7.9 % | 6.1 % |
| CoR `α` | 36.8 % | 5.3 % | 8.0 % | 4.9 % | 1.2 % | 3.7 % |
| `δ/R` (max depth) | 25.4 % | 12.3 % | 5.6 % | **8.9 %** | 11.9 % | 2.9 % |

- Oil contact time is the headline: closer to the measurement at every one of the 12 points than either published model. It also reproduces the flat dependence on `We`.
- Two open failures share a suspect: in oil (`Oh = 0.058`) `α` has the wrong trend (measured rises 0.237 → 0.269, model falls 0.177 → 0.139) and `δ` is under-predicted, while the contact timescale is right. This points at the viscous coupling, not the contact closure.
- Overall it does **not** supersede 1PKM: across 51 point comparisons it is closer at 16 and 1PKM at 35 (asserted in `test/test_fluid_sweep_results.jl`).
- **The sides measure different things.** The experiment times the north pole crossing `z = 2R`; `threshold_contact_time` uses the centre of mass crossing `z = R` (the authors quote a 5 % water / 2 % oil difference), so the water residual is of the order of the definitional offset.
- **Provenance limit.** The experimental CSVs come from `.fig` files whose published-PDF counterparts are unavailable here; conclusions resting on them (notably the oil `δ` deficit) inherit that. Conclusions about the solver do not.
- Against measured droplet top and bottom trajectories (99 points) the RMS error is 1.4 experimental error bars. No measured contact radius or width exists, so those comparisons are model-to-DNS only; the contact radius peaks early (`τ ≈ 0.89` against DNS `1.46`; 1PKM `1.22`) for reasons not yet resolved.
- `data/experiments/` also holds a droplet-on-**solid** dataset (`SOLID_SUBSTRATE_*`, contact times 1.96–3.36 against the bath's 4.44–8.09). It is a different problem and is guarded against accidental use.

## Sweeps, rendering, tests

```bash
julia --project=. -t auto scripts/sweep.jl --wall=free 0.2 0.4 1.0958 3.0   # worker count chosen by measurement; resumable
julia --project=. scripts/make_video.jl [We] [outfile.mp4] [M] [L] [N] [nq] # MP4 to output/media/, no plotting dependency
julia --project=. -e 'using Pkg; Pkg.test()'
```

A typical impact takes about 25 s at `M = L = 60`, roughly flat in `We`. CI runs the tests on every push to `main` and `dev`. The experiment comparison (about 3.5 min per point) is dispatched by hand: `gh workflow run CI -f fluid=oil [-f selector=crossing] [-f viscous=lamb]`. Tests assert physical facts (ballistic flight, Rayleigh–Lamb frequency, capillary–gravity dispersion, momentum balance to 1.3e-3), not just internal identities.

## Derivation audits

`derivations/` scripts back specific claims and are meant to be run, not trusted: `verify_formulation.py` (SymPy: reads the equations out of `FORMULATION.md` and checks they follow from first principles and from each other; a sign error in the text fails it), `verify_gauge_shift.py` (a constant pressure shift is a gauge change and cannot relax `p ≥ 0`), `verify_edge_pinning.py` with `verify_edge_pinning_model.jl` (track A4), `verify_legendre_pressure_basis.jl` (pressure basis degree, exact in `Rational{BigInt}`), `audit_compliance_operator.jl` (self-adjointness, rank), `audit_nested_closure.jl` (inner conditioning flat in `δ` against ~1e18 for the joint system), `audit_weak_determinacy.jl` (which functionals converge in `N`), `feasibility_pinned_contact_line.jl` (wall conditions), and `probe_global_basis_argmin.jl` (a formulation that **fails**, recorded for the record).
