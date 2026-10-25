module SpectralKM

using LinearAlgebra
using SpecialFunctions
using ForwardDiff

include("legendre.jl")
include("bessel.jl")
include("quadrature.jl")
include("reid.jl")
include("fluids.jl")
include("types.jl")
include("geometry.jl")
include("pressure.jl")
include("bessel_moments.jl")
include("affine.jl")
include("residual.jl")
include("newton.jl")
include("timestepper.jl")
include("signorini.jl")
include("tangency.jl")
include("postprocessing.jl")

export Params, Level, SimHistory, Phase, FreeFlight, InContact
export Fluid, WATER, OIL_5CST, fluid, conditions, ImpactConditions
export dynamic_viscosity, t_sigma, STANDARD_GRAVITY
export resolvable_rank_estimate, min_nq_for_exact_com
export reid_root, reid_root_tracked, reid_real_roots, reid_pole_pair
export sph_bessel_ratio, reid_first_singularity
export lamb_eigenvalue, drop_viscous_coeffs
export legendre_P, legendre_dP, bonnet_H
export bessel_zeros_J1, bessel_zeros_J0
export gauss_legendre_nodes
export pressure_poly, b_l_all, com_force_closed
export forward_map_r, forward_map_zd, xi_of_theta, w_of_x, r_of_x, zd_of_x, xi_of_x
export c_m_all
export fitted_rebound_metrics, north_pole_height
export bath_affine, drop_affine, com_affine, bdf2_coeffs
export residual, unpack_state, C_at_theta, tangency_residual
export check_nonintersect, check_monotone_r, check_positivity
export crossing_theta_c, select_theta_c_crossing
export select_theta_c_tangency, residual_pinned, inner_solve_pinned, max_overlap_beyond
export select_signorini, signorini_grid, signorini_assemble, signorini_theta_max, nnqp
export newton_solve, minnorm_gauss_newton, NewtonResult, NewtonStatus, Converged, MaxIterExceeded, Stalled
export run_simulation, contact_step, inner_solve, onset_theta_c, free_flight_step
export initial_level, gap_at_pole, select_theta_c, feasible_at, contact_quad
export coefficient_of_restitution, contact_time, max_penetration_depth
export contact_intervals, primary_contact_time, threshold_contact_time
export threshold_crossings, threshold_coefficient_of_restitution
export reconstruct_bath, reconstruct_drop

end # module SpectralKM
