# Method notes, and where this reproduction differs from the paper

## What is reproduced exactly

* **Equations.** Ideal isothermal MHD in a periodic unit box, mean density and
  sound speed unity.  The paper uses Enzo with an adiabatic index of
  gamma - 1 = 1e-3 "to achieve a reasonable approximation to an isothermal
  equation of state"; we solve the isothermal system directly, which is the
  limit they are approximating and avoids the catastrophic cancellation that a
  gamma = 1.001 energy equation suffers.
* **Numerics.** Piecewise-linear reconstruction, the HLLD Riemann solver for
  isothermal MHD (Mignone 2007 - the solver reference the paper cites), and
  Dedner et al. (2002) mixed hyperbolic/parabolic divergence cleaning.  Time
  integration is SSP-RK2 at CFL 0.4 with an unsplit (method-of-lines) update.
* **Driving.** Stochastic Ornstein-Uhlenbeck acceleration on the 32 Fourier
  modes with 1 <= |k|/k_min <= 2, correlation time t_dyn = L_0/v_rms with
  L_0 = L/2, Helmholtz projection P_ij = zeta d_ij + (1-2 zeta) k_i k_j / k^2
  with zeta = 1/2, and an amplitude rescaled every step so that the injected
  kinetic energy is exactly Edot dt (Mac Low & Klessen 2004).  The driving
  field carries no net helicity and no net momentum, so the setup respects
  parity, as in the paper.
* **Parameters.** M_S = v_rms/c_s and M_A = v_rms/v_A with v_A = B_0/sqrt(rho),
  B_0 along x-hat; the grid is the paper's 7 x 3 = 21 combinations of
  M_S = 0.5, 1, 2, 3, 4, 5, 6 and M_A = 0.5, 1, 2.  Ten (here twelve) dynamical
  times are run and only the second half is analysed.
* **Observables.** T = int rho dz, Q and U from Eqs. (7)-(8), E/B from the
  flat-sky rotation (Eq. 9), shell-averaged auto- and cross-spectra, power-law
  fits C_k = A k^alpha over an inertial window, correlation coefficients
  r_k^XY, and the linear model alpha_q = a + b M_S + c M_A of Table 1.

## Deliberate differences

1. **Sign of Q, U.**  Eqs. (7)-(8) as printed use psi = the angle of the
   *magnetic field* in the plane of the sky.  Dust grains emit with their
   polarization perpendicular to the field, so the physical Stokes parameters
   carry the opposite sign, cos 2 psi_pol = -cos 2 psi_H.  Taken literally the
   printed equations flip the sign of E, B and r_TE.  We use the physical
   convention: it is the one that gives the positive T-E correlation that both
   the paper and Planck report, and it passes the analytic filament tests in
   `python/test_eb.py` (a filament threaded by a parallel field gives pure E
   with r_TE = +1; a field at 45 degrees gives pure B).

2. **Energy injection rate.**  The paper sets Edot proportional to M_S^3/L.
   That alone misses the target M_S badly at low M_A: a strong mean field makes
   the cascade much less dissipative, so the same Edot yields a ~60% larger
   v_rms at M_A = 0.5 than at M_A = 2.  We calibrated an extra M_A^1.5 factor
   from a first pass of the suite.  Edot is constant throughout each run, as in
   the paper, and every result is plotted against the *measured* Mach numbers.

3. **Resolution.**  The paper uses 512^3.  On a single 16-thread workstation
   the suite is run at 64^3 and 128^3, with 256^3 for convergence.  This is the
   one difference that matters quantitatively: with piecewise-linear
   reconstruction, numerical dissipation sets in near 20 cells per wavelength,
   i.e. k_diss ~ N/20 - so the paper's fit window k/k_min in [4,25] is only
   available at N >= 512.  At lower resolution the driving range (k <= 2) and
   the dissipative range nearly touch, the fit window is short, and *all*
   spectral slopes are biased steep.  Slope trends (the b and c coefficients),
   amplitude ratios and correlation coefficients are much less sensitive to
   this than the absolute slopes; the convergence table quantifies it.
   `cluster/run_suite.slurm` runs the identical code at 512^3.

4. **Grid-refinement restarts.**  Higher-resolution runs are started from the
   trilinearly interpolated final state of the next coarser run, and then
   evolved for several dynamical times before analysis.  The large scales are
   already in statistical equilibrium, so only the new small scales have to
   fill in; this costs ~1/3 of a run started from rest.

5. **Analysis window.**  We fit k/k_min in [3, N/20] (floored at 5), the same
   "20 cells per wavelength" dissipation criterion the paper's [4,25] at 512^3
   corresponds to.  Slope uncertainties are the scatter over snapshots and over
   variations of both window edges, as in the paper.
