# A fully spectral kinematic match for droplet–bath impact: model formulation

*Text and notation follow §3 of Alventosa, Cimpeanu & Harris (2023, *J. Fluid Mech.* **958**, A24; arXiv:2209.13276) with minimal edits. Their sentences are kept wherever the model is unchanged. The departures are: (i) the single-point match and the assumed pressure profile are replaced by a spectral pressure solved for over the whole contact region (§3.3–3.4); (ii) the bath radius $b$ is kept as a model parameter of a few droplet radii (Alventosa et al. use the same cylinder with $b=25R$), the finite-radius derivation is written out (§3.1), and the Fourier–Bessel normaliser is corrected; (iii) BDF2 replaces backward Euler, and an optional exact viscous closure replaces Lamb's for the droplet. Derivations and audits: `derivations/`.*

## 3. Linearized quasi-potential fluid model

In this section, we develop a model for droplet impact on a flat fluid interface from first principles. First, we use the linearized Navier–Stokes equations to model the flow within the bath and the shape of the bath interface. Then, we use an orthogonal function decomposition of the bath model to derive a single set of linear ordinary differential equations (ODEs) that govern the bath mode amplitudes. We then write a similar model for the droplet interface, which reduces to another set of linear ODEs governing droplet mode amplitudes. Finally we represent the pressure distribution and its extent acting on the bath and droplet during contact by a spectral series, couple the two sets of ODEs together using a kinematic match condition enforced over the whole contact region, and solve the system implicitly.

### 3.1. Bath interface model

The bath interface dynamics are modelled using a linearized, quasi-potential flow model following Galeano-Rios et al. (2017) and Dias et al. (2008). Assuming the flow to be incompressible, isothermal, and Newtonian, we can define the fluid velocity vector $\mathbf u=\nabla\phi+\nabla\times\boldsymbol\Psi$ and the bath interface shape $\eta=\eta(r,t)$. Here, $\phi$ is the scalar potential and $\boldsymbol\Psi$ is the vector stream function. We then linearize the governing equations and boundary conditions about the undisturbed free surface $z=0$, and recast them as

$$
\nabla^2\phi+\partial_z^2\phi=0,\quad z\le0;\qquad \partial_t\eta=\partial_z\phi+2\nu\nabla^2\eta,\quad z=0,\qquad\text{(3.1–3.2)}
$$

$$
\partial_t\phi=-g\eta-2\nu\,\partial_z^2\phi+\sigma\kappa/\rho-p_s/\rho,\quad z=0;\qquad \partial_z\phi=0\ \mathrm{at}\ z=-h_0.\qquad\text{(3.3–3.4)}
$$

Here, $p_s(r,t)$ is the contact pressure, $g$ is the gravitational acceleration, $\kappa=\nabla^2\eta$ is twice the linearized mean curvature of the interface, $h_0$ is the depth of the undisturbed bath, and $\rho,\sigma,\nu=\mu/\rho$ are the fluid density, surface tension, and kinematic viscosity, respectively. The tangential stress boundary conditions are automatically satisfied in these approximations. This leading-order theory is valid in the weakly viscous limit when $\mathrm{Re}\gg1$, where $\mathrm{Re}$ is $\mathrm{We}^{1/2}$ over $\mathrm{Oh}$.

We assume that the impact occurs in a bath of some viscous fluid contained in a cylinder of finite radius $b$, subject to two boundary conditions, $\partial_n\phi=0$ on the walls of the bath and $\partial_z\phi=0$ on the bottom of the bath, where $n$ is the outward facing normal of the walls (Benjamin & Ursell 1954). The former condition implies that $\partial_n\eta=0$ on the walls. These conditions correspond physically to a bath where the working fluid maintains a constant contact angle of $90^\circ$ at the walls, and has no-flux boundary conditions along the walls and bottom. For the axisymmetric problem the orthogonal basis functions are $S_{0,m}=J_0(k_mr)$, with $J_0$ the Bessel function of the first kind and $k_m$ the solutions of $J_0'(k_mb)=0$, i.e. $J_1(k_mb)=0$, together with the piston mode $k_0=0$. We then can express the free surface elevation as

$$
\eta(r,t)=\sum_{m=0}^{M}a_m(t)J_0(k_mr),\qquad\text{(3.6)}
$$

Each $J_0(k_mr)$ satisfies $(\nabla^2+k_m^2)J_0=0$ and, because $J_0'(k_mb)=0$, the wall condition at every depth. Writing $\phi=\sum_m\phi_m(t)J_0(k_mr)Z_m(z)$, Laplace's equation (3.1) gives $Z_m''=k_m^2Z_m$, and the bottom condition (3.4) selects $Z_m=\cosh k_m(z+h_0)/(k_m\sinh k_mh_0)$, normalised so that $Z_m'(0)=1$. The kinematic condition (3.2) then fixes $\phi_m=\dot a_m+2\nu k_m^2a_m$, so

$$
\phi=\sum_{m=0}^{M}\big(\dot a_m+2\nu k_m^2a_m\big)J_0(k_mr)\,Z_m(z),\qquad \kappa=-\sum_{m=0}^{M}k_m^2a_mJ_0(k_mr).\qquad\text{(3.9)}
$$

The functions $J_0(k_mr)$ are orthogonal under the weight $r$ on $[0,b]$, with $\int_0^bJ_0(k_mr)^2\,r\,dr=\tfrac12b^2\big[J_0(k_mb)^2+J_1(k_mb)^2\big]=\tfrac12b^2J_0(k_mb)^2$ since $J_1(k_mb)=0$. Substituting (3.9) into the dynamic condition (3.3), using $\partial_z^2\phi=k_m^2\phi$ mode by mode, projecting onto $J_0(k_mr)$, multiplying by $k_m\tanh(k_mh_0)$ and dropping the $O(\nu^2)$ term $4\nu^2k_m^4a_m$ (consistent with the weakly viscous limit), each wave mode in the bath is described by a forced, damped harmonic oscillator equation,

$$
\ddot a_m+4\nu k_m^2\dot a_m+\Big(\sigma k_m^2/\rho+g\Big)k_m\tanh(k_mh_0)\,a_m=-\,(d_m/\rho)\,k_m\tanh(k_mh_0),\qquad\text{(3.10)}
$$

where $d_m$ are the Fourier–Bessel coefficients of the pressure, $p_s=\sum_m d_mJ_0(k_mr)$,

$$
d_m=\frac{2}{(bJ_0(k_mb))^2}\int_0^b p_s(r)\,r\,J_0(k_mr)\,dr.\qquad\text{(3.21)}
$$

The preprint prints $J_1$ in this normaliser; $J_0(k_mb)$ is what the orthogonality relation above gives. The finite radius enters only through the discrete spectrum $J_1(k_mb)=0$ and this norm. Two consequences are used later. The piston mode $k_0=0$ is the only mode that carries volume ($\int_0^bJ_0(k_mr)\,r\,dr=bJ_1(k_mb)/k_m=0$ for $k_m\ne0$), and its forcing coefficient $k_0\tanh(k_0h_0)$ vanishes, so bath volume is conserved exactly. As $b\to\infty$ the zeros of $J_1$ become evenly spaced by $\pi/b$ and the sums tend to Hankel integrals, the deep-bath limit that $b=25R$ approximates; at $b$ of a few $R$ the discrete spectrum is a physical property of the container.

### 3.2. Droplet interface model

Additionally, we wish to recover a similar set of equations that describe the gravity-capillary waves present in the droplet, and then couple these equations to the motion of the bath. The full derivation of the droplet oscillation model can be found throughout prior works (Lamb 1924; Tsamopoulos & Brown 1983; Courty et al. 2006; Chevy et al. 2012; Balla et al. 2019) and is briefly summarized below. We begin by utilizing spherical harmonics to decompose the droplet radius in a spherical domain; due to the axisymmetry of the problem they reduce to Legendre polynomials, $\xi(\theta,t)=R+\sum_l\beta_l(t)P_l(\cos\theta)$. We assume that the velocity potential takes the same form as the decomposition of the interface (Lamb 1924; Balla et al. 2019), and turn to an energy conservation equation of the form $dT/dt=\dot W-\varepsilon$, with $T=K+G$, $\dot W$, and $\varepsilon$ the total energy of the drop, the rate of work done on the droplet interface, and the viscous dissipation, respectively. Utilizing the linearized kinematic boundary condition yields a set of forced, damped harmonic oscillators that describe the amplitude of each individual spherical mode,

$$
\ddot\beta_l+2\alpha_l\dot\beta_l+\omega_l^2\beta_l=-\frac{(2l+1)l}{2\rho R}\,\Pi_l+g\,\delta_{1l},\qquad\text{(3.14)}
$$

with $\alpha_l=(2l+1)(l-1)\mu/(\rho R^2)$ and $\omega_l^2=l(l-1)(l+2)\,\sigma/(\rho R^3)$. Here $\Pi_l$ is the projection of the pressure on mode $l$. In Alventosa et al. it is $\int_0^\pi p_s\sin\theta P_l\,d\theta$; we use $\Pi_l=\int_{x_c}^1p_s\,x\,P_l\,w\,dx$ with $x=\cos\theta$ and $w=-\tfrac12\,d(r^2)/dx$ (the area element, $r\,dr=w\,dx$). This differs from theirs at $O(\theta_c^2)$, reduces to it as $\theta_c\to0$, and makes the contact operator of §3.4 self-adjoint. The mode $l=1$ is the translational mode, i.e. the centre of mass, $\beta_1=-z_{cm}$. This model is valid in the weakly viscous limit, when $\mathrm{Oh}\ll1$; as an option (`viscous=:reid`) $\alpha_l,\omega_l$ are replaced by the exact roots of Reid's (1960) characteristic equation, valid at arbitrary $\mathrm{Oh}$.

### 3.3. Pressure forcing during impact

There is still an additional unknown in the bath mode (3.10) and drop mode (3.14) equations: the applied pressure distribution $p_s(r,t)$. However, in this work, we will assume that the viscosity of the ambient gas is small relative to the fluid bath such that the flow within the small air film is negligible, and the pressure acts solely to apply an upwards hydrodynamic force on the droplet. The radial extent of the pressure distribution is generally unknown for impact problems, and constitutes an additional problem that we must solve. The contact region, $A_c$, will be assumed to be a simply-connected disk, following Galeano-Rios et al. (2017) and Korobkin (1995), which allows us to write a single unknown $\theta_c(t)$, with $x_c=\cos\theta_c$, to fully describe the temporal evolution of this region of contact.

Unlike Alventosa et al., who assume a profile $p_s=F(t)H_r(r/r_c)/\pi R^2$ with $H_r=C[1-(r/r_c)^6]$, we assume no spatial profile. With $\psi=(x-x_c)/(1-x_c)\in[0,1]$, the pressure is a shifted-Legendre series in $x$,

$$
p_s=\frac{2\sigma}{R}\sum_{n=0}^{N}\hat c_n(t)\,\tilde P_n(\psi)\quad(x\ge x_c),\qquad p_s\equiv0\quad(x<x_c),\qquad\text{(3.17)}
$$

with $\tilde P_n(\psi)=P_n(2\psi-1)$. Each $\tilde P_n$ is exactly a degree-$n$ polynomial in $x$, the basis is orthogonal on the patch (Gram condition number $2N+1$, against Hilbert-matrix growth for monomials), and the vanishing of the pressure outside the patch is exact and costs no equations. The coefficients $\hat c_n$ are unknowns, determined by the kinematic match below. Substituting (3.17) into (3.10) and (3.14), the pressure enters only through its moments $c_m$ (the dimensionless $d_m$), $b_l$ (the dimensionless $\Pi_l$) and the net force $f$, computed by Gauss–Legendre quadrature. The droplet and centre-of-mass integrands are polynomial in $x$, so the quadrature is exact once the number of nodes $n_q$ exceeds half their degree; the bath integrals contain $J_0(k_mr(x))$ and are quadrature-approximate. The $l=1$ mode reduces to Newton's law for the centre of mass, $\ddot z_{cm}=-F/m+g$ (3.25), with $F$ the integral of the pressure over the contact area.

### 3.4. Modeling contact

The contact pressure is determined through a kinematic match condition. Alventosa et al. enforce contact only at a single point, the centre of the axisymmetric domain (the "1-Point" kinematic match, 1PKM), $\eta(0,t)=z_{cm}(t)-\xi(0,t)$, which determines the unknown contact force $F(t)$. We instead enforce contact throughout the contact region, the full kinematic match of Galeano-Rios et al. (2017, 2019, 2021) and Agüero et al. (2026), which circumvents the need for any assumption on the pressure profile shape. Writing $r(\theta)=\xi\sin\theta$, $z_d(\theta)=\xi\cos\theta$ for the droplet surface in the bath frame,

$$
C(\theta)\equiv\eta(r(\theta),t)-z_{cm}(t)+z_d(\theta)=0,\qquad0\le\theta\le\theta_c.\qquad\text{(3.26)}
$$

With $N+1$ coefficients this cannot hold pointwise, so it is imposed in the weighted-residual sense by testing against the pressure basis,

$$
\int_{x_c}^{1}C(\theta)\,\tilde P_n(\psi)\,w\,dx=0,\qquad n=0,\dots,N,\qquad\text{(3.26')}
$$

a square system for $\hat c_n$ at fixed $\theta_c$. Contact between the bath and drop ends when the magnitude of the contact force as predicted by the kinematic match becomes negative (here: non-positive for several consecutive steps, or the patch collapsing to its floor).

Lastly we turn to the unknown contact extent $\theta_c(t)$. The exact closure is a unilateral (Signorini) contact condition: the pressure is non-negative, the surfaces may not overlap, and pressure acts only where they touch,

$$
p\ge0,\qquad C\le0,\qquad p\,C=0\qquad(0\le\theta\le\pi).\qquad\text{(3.26a)}
$$

At frozen geometry $C=C_{\rm free}+\mathcal Ap$ (§3.6), with $C_{\rm free}$ the mismatch produced by the time-stepping history alone and $-\mathcal A$ symmetric positive semidefinite in $\langle u,v\rangle_w=\int uv\,w\,dx$. Then (3.26a) are precisely the optimality conditions of the convex problem

$$
\min_{p\ge0}\ \tfrac12\langle p,(-\mathcal A)p\rangle_w-\langle p,C_{\rm free}\rangle_w,\qquad\text{(3.26b)}
$$

so the gap response $\mathcal Ap$ is unique, and the contact set $\{p>0\}$ is part of the solution rather than an input. Tangency, $\partial_\theta C(\theta_c)=0$, is the usual consequence at the edge of the set, but it holds identically at $\theta=0$ by axisymmetry and so cannot select $\theta_c$ at onset.

What we solve is not (3.26b). We restrict the contact set to a single disk, carry the pressure as the polynomial (3.17), and choose its edge as Alventosa et al. do: the predictions of both interfacial models, advanced one step without pressure, are used to find the exact location where the two interfaces cross, and the largest crossing is the instantaneous edge, held at its previous value when there is none, with $\theta_c$ clamped to $[\arcsin0.01,\arcsin0.9998]$. Since $C$ is the vertical gap between the surfaces at radius $r(\theta)$, a root of $C$ is an intersection of the two curves. The rule is a predictor for the support of $p$: pressure pushes the surfaces apart, so the region it must press is bounded by where the unpressed surfaces overlap. This holds for the centre-of-mass and low-order bath contributions to $\mathcal A$ but is not guaranteed pointwise for its oscillatory high-wavenumber part. It is therefore a closure inherited from 1PKM, not a derived free-boundary condition, and it carries no error bound. An alternative, $\theta_c=\inf\{\theta:\text{no overlap beyond }\theta,\ \partial_xr<0\}$, is available as `selector=:feasible`; it is degenerate once non-overlap stops binding and is retained only to reproduce earlier results.

Consistency of the chosen $\theta_c$ with (3.26a) is therefore checked rather than assumed. The tangency residual $\partial_\theta C(\theta_c)$ is reported at every contact step, and $\min p$ on the patch and the overlap beyond the edge can be evaluated from the stored state. At $M=L=60$, $N=3$, negative pressure occurs on 2% of contact steps (90% at $N=12$, past the resolvable budget of §3.6). The contact radius peaks earlier than in the DNS ($\tau\approx0.89$ against $1.46$), so the closure is open to improvement. Solving (3.26b) directly, which would determine the contact set rather than assume it, is the natural next step and has not been done.

### 3.5. Summary

Choosing a timescale of $t_\sigma=\sqrt{\rho R^3/\sigma}$ (with $\tau=t/t_\sigma$), a length scale of $R$, a pressure scale of $2\sigma/R$, and a force scale of $2\pi\sigma R$ (with $f=F/2\pi\sigma R$), we recast the governing equations in non-dimensional form,

$$
\eta=\sum_{m=0}^{M}a_m(\tau)J_0(k_mr),\qquad \xi=1+\sum_{l=2}^{L}\beta_l(\tau)P_l(\cos\theta),\qquad\text{(3.27–3.28)}
$$

$$
\ddot a_m+2\gamma_m\dot a_m+\omega_m^2a_m=-2\,c_m\,k_m\tanh(k_mh_0),\qquad\text{(3.29)}
$$

with $\gamma_m=2\,\mathrm{Oh}\,k_m^2$ and $\omega_m^2=(k_m^2+\mathrm{Bo})\,k_m\tanh(k_mh_0)$,

$$
\ddot\beta_l+2\gamma_l\dot\beta_l+\omega_l^2\beta_l=-\,b_l\,(2l+1)\,l,\qquad\text{(3.30)}
$$

with $\gamma_l=\mathrm{Oh}\,(2l+1)(l-1)$ and $\omega_l^2=l(l-1)(l+2)$,

$$
\ddot z_{cm}=\tfrac32f-\mathrm{Bo},\qquad\text{(3.31)}
$$

$$
\int_{x_c}^{1}\big[\eta(r(\theta),\tau)-z_{cm}(\tau)+z_d(\theta)\big]\tilde P_n(\psi)\,w\,dx=0,\quad n=0,\dots,N,\qquad\text{(3.32)}
$$

with $c_m=\frac{2}{(bJ_0(k_mb))^2}\int p\,J_0(k_mr)\,r\,dr$, $b_l=\int_{x_c}^1p\,x\,P_l\,w\,dx$ and $f=2\int_{x_c}^1p\,w\,dx$, where $p$ is the pressure (3.17) scaled by $2\sigma/R$. Equations (3.29) and (3.30) describe the evolution of the bath and droplet oscillation modes, respectively. Equation (3.31) governs the vertical motion of the droplet's center of mass. Equation (3.32) couples these equations all together and, with the edge rule of §3.4, allows for determination of the unknowns $\hat c_n(\tau)$ and $\theta_c(\tau)$. The shape of the bath and droplet can be reconstructed at any time via the sums in (3.27) and (3.28). The complete model is valid when $\mathrm{Re}\gg1$ (with $\mathrm{Re}$ as above); $\mathrm{Oh}\ll1$ for the Lamb closure only. Also, since the model is linearized about the undeformed state, we anticipate it to hold when deformations remain small, further suggesting $\mathrm{Bo}\ll1$ and $\mathrm{We}\ll1$.

### 3.6. Numerical methods

We solve these equations using an implicit method, as the instantaneous size of the pressure distribution acting on both the droplet and bath at the next time step is unknown (Galeano-Rios et al. 2021); treating the pressure explicitly on either the droplet or the bath can lead to nonphysical behavior in the system. We use variable-step second-order backward differentiation (BDF2) rather than backward Euler. Each oscillator $\ddot u+2\gamma\dot u+\omega^2u=Fq$ of (3.29)–(3.31) then makes its state affine in its own pressure moment at the new level,

$$
u^{k+1}=\mathrm{history}+\kappa\,q^{k+1},\qquad \kappa=\delta^2F/(\mathfrak a^2+2\gamma\mathfrak a\,\delta+\delta^2\omega^2),\qquad\text{(3.33)}
$$

with step $\delta$ and BDF2 leading coefficient $\mathfrak a=(1+2s)/(1+s)$ for step ratio $s$. Hence $\eta$, $\xi$ and $z_{cm}$ are explicit affine functions of $\hat c$, and each time step reduces to the algebraic system (3.32). Writing $C=C_{\rm free}+\mathcal Ap$ at frozen geometry, the compliance operator $\mathcal A$ is self-adjoint and negative semidefinite in $\int uv\,w\,dx$, with a small numerical rank (7, 9, 14 at $M=L=20,40,80$). Hence refining $N$ at fixed $M,L$ does not converge the pointwise pressure, whereas the moments that enter (3.29)–(3.31) converge quickly, and $N$ is kept inside this budget (default $N=3$). Because $\theta_c$ is held fixed during each inner solve, the inner system has condition number of order one and flat in $\delta$, against order $10^{18}$ for a joint $(\hat c,\theta_c)$ system. It is solved by Levenberg–Marquardt-damped Newton with a forward-mode automatic-differentiation Jacobian. Defaults are $M=60$, $L=120$, $N=3$, $n_q=200$, with a bath of radius $b=6$ (Alventosa et al.: $M=150$, $L=55$, $b=25R$, no pressure unknowns). Convergence is judged on contact time.

*Wall condition.* The no-flux basis above (`wall=:free`) conserves bath volume exactly, since every mode with $k_m\ne0$ is volume-neutral and the piston mode cannot be driven. A pinned contact line, $\eta(b)=0$, is imposed (`wall=:clamped`) on the same basis by a scalar constraint $\sum_ma_mJ_0(k_mb)=0$ carried by a Lagrange multiplier (the rim's line force), a rank-one correction that preserves self-adjointness and volume conservation. A Dirichlet basis pins $\eta(b)$ but does not conserve volume (`wall=:pinned`, diagnostic only). At $b=6$ the wall condition changes the coefficient of restitution by $-0.001$; for smaller baths it matters.

*Performance against experiment.* Against the measured bath-rebound data of Alventosa et al., contact time is predicted to 2.7% mean error in oil (1PKM 11.0%, DNS 3.3%), while the coefficient of restitution and penetration depth are under-predicted in oil. Across all 51 point comparisons the model is closer than 1PKM at 16 (`README.md`; `test/test_fluid_sweep_results.jl`).

## References

- Agüero, E., Galeano-Rios, C. A., Ragazzo, C., Gabbard, C. T., Harris, D. M. & Milewski, P. A. 2026 Droplet rebounds off a fluid bath at low Weber numbers. *J. Fluid Mech.*
- Alventosa, L. F. L., Cimpeanu, R. & Harris, D. M. 2023 Inertio-capillary rebound of a droplet impacting a fluid bath. *J. Fluid Mech.* **958**, A24.
- Dias, F., Dyachenko, A. I. & Zakharov, V. E. 2008 Theory of weakly damped free-surface flows: a new formulation based on potential flow solutions. *Phys. Lett. A* **372**, 1297–1302.
- Galeano-Rios et al. 2017, 2019, 2021; Korobkin 1995; Benjamin & Ursell 1954; Lamb 1924; Reid 1960; and the droplet-oscillation references (Tsamopoulos & Brown 1983; Courty et al. 2006; Chevy et al. 2012; Balla et al. 2019), as cited in Alventosa et al. (2023).
