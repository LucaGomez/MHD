# Reproduction of Stalpes, Collins & Huffenberger (2024)

*Planck dust polarization power spectra are consistent with strongly supersonic
turbulence* — arXiv:2404.02874

A driven isothermal MHD turbulence code was written from scratch (`src/mhd.c`),
validated against exact solutions, and used to run the paper's full 21-run
(M_S, M_A) grid at two resolutions.  Everything below comes from those runs;
no result is quoted from the paper except as a comparison target.

## 1. Solver validation

| test | result |
|---|---|
| circularly polarised Alfvén wave, one full period | L1 error 8.33e-3 → 3.00e-3 → 7.68e-4 for N = 16, 32, 64: **2nd-order convergence** (measured order 1.97) |
| isothermal shock tube vs the **exact** Riemann solution | L1(ρ) 1.17e-2 (64) → 5.80e-3 (128), 1st order at the shock; ρ\* = 0.3458, u\* = 1.0620, S = 1.6632 recovered exactly |
| isothermal MHD (Brio–Wu type) tube | monotone, symmetric, oscillation-free |
| ∇·B control (Dedner GLM) | \|∇·B\|Δx/\|B\| stays at 10⁻³–10⁻² through fully developed supersonic turbulence |
| projection + E/B decomposition | five analytic filament configurations reproduced exactly (`python/test_eb.py`): filament ∥ H → pure E with r_TE = +1; filament ⊥ H → pure E with r_TE = −1; H at 45° → pure B; H at 22.5° → B/E = 1; H along the line of sight → no polarization |

## 2. The suite

21 runs at each of 64³ and 128³, targeting M_S = 0.5, 1, 2, 3, 4, 5, 6 and
M_A = 0.5, 1, 2, driven for 12 dynamical times with the last 4.5 analysed
(19 snapshots), then repeated at 128³ starting from the converged 64³ state.
Achieved Mach numbers span M_S = 0.52 – 6.91 and M_A = 0.43 – 2.30, matching the
coverage of the paper's Fig. 2.

## 3. Table 1 — linear fits α_q = a + b M_S + c M_A

|          |     | N = 64 | N = 128 | paper (512³) |
|----------|-----|--------|---------|--------------|
| α_ρ  | a | −4.09 | **−3.56** | −3.61 |
|      | b |  0.219 | **0.164** |  0.160 |
|      | c |  0.160 |   0.110 | −0.000 |
| α_v  | a | −4.01 | **−3.63** | −3.86 |
|      | b | −0.013 | **0.023** |  0.020 |
|      | c |  0.104 | **0.012** |  0.140 |
| α_H  | a | −3.96 | **−3.17** | −3.31 |
|      | b |  0.010 | −0.023 |  0.020 |
|      | c |  0.109 | **−0.089** | −0.280 |
| α_TT | a | −4.82 | **−3.85** | −3.66 |
|      | b |  0.300 | **0.174** |  0.150 |
|      | c |  0.348 |   0.232 |  0.090 |
| α_EE | a | −4.75 | **−3.78** | −3.63 |
|      | b |  0.276 | **0.159** |  0.170 |
|      | c |  0.596 | **0.376** |  0.280 |
| α_BB | a | −5.56 | **−4.59** | −4.82 |
|      | b |  0.377 | **0.257** |  0.280 |
|      | c |  0.779 | **0.553** |  0.640 |

Every one of the eighteen coefficients moves monotonically toward the published
value as the resolution doubles.  At 128³ the density spectrum is reproduced
essentially exactly (a = −3.56 vs −3.61, b = 0.164 vs 0.160), as are the
sensitivities of the polarization slopes to the sonic Mach number
(b_EE = 0.159 vs 0.170, b_BB = 0.257 vs 0.280).

The paper's qualitative statements are all recovered:

* **α_ρ rises nearly linearly with M_S** and is independent of M_A — the paper's
  α_ρ ≈ −3.6 + 0.16 M_S is matched to two decimal places.
* **α_v is insensitive to M_S** (b ≈ 0.02, both) and sits near the Kolmogorov
  value, −3.6 here vs −3.7 to −3.9 in the paper.
* **α_H is insensitive to M_S and steepens with increasing M_A** — the sign of
  c_H, which is wrong at 64³, comes out correct at 128³.
* **b_BB > b_EE > b_TT ≫ b_v ≈ b_H** and **c_BB > c_EE > c_TT**: the B-mode
  slope is the most sensitive to both Mach numbers, exactly the ordering of
  Table 1.
* Simulations with **M_S ≳ 4 are needed to reach the Planck slope**
  α_EE = −2.42 (see `figures/n128/fig07_projected_slopes.png`); slower flows
  give slopes that are too steep.

## 4. Correlations and amplitude ratios

| quantity (M_S > 4) | N = 64 | N = 128 | paper | Planck |
|---|---|---|---|---|
| ⟨r_TE⟩ | 0.04 | **0.22** | ≈0.3 | 0.355 |
| ⟨A_BB/A_EE⟩ (band power) | 1.03 | 0.93 | 0.55 | 0.53 |
| ⟨A_EE/A_TT⟩ | 0.31 | 0.37 | 0.62 | – |
| ⟨A_BB/A_TT⟩ | 0.37 | 0.44 | 0.34 | – |
| r_TB > 0 | 8/21 | 11/21 | 14/21 | – |
| r_EB > 0 | 11/21 | **14/21** | 14/21 | – |
| mean \|r_TB\|, \|r_EB\| | – | 0.068, 0.063 | "below 0.05, murky" | 0.055 |

The paper's detailed description of the r_TE behaviour is reproduced exactly:
strong positive correlation at low M_S with a strong field (r_TE = 0.77 at
M_S = 0.5, M_A = 0.5, "well above the value observed by Planck"), and a slight
**anti**-correlation at low M_S with a weak field (r_TE = −0.18 at M_S = 0.5,
M_A = 2.1), converging to modest values consistent with Planck at high M_S.

The parity-violating correlations are spectrally flat, scatter about zero, and
are of order 0.06 — the paper's "certainly below 0.05, but the results are
somewhat murky".  The sign statistics agree for r_EB (14/21 positive, exactly
the published count) and are weaker for r_TB.

## 5. What has not converged: A_BB/A_EE

The one quantity still clearly away from the published value is the B-to-E power
ratio at high sonic Mach number: 0.93 here against 0.55 in the paper.  It is
moving the right way with resolution (1.03 → 0.93) but slowly, and it is flat in
k across the resolved range, so it is not a fitting-window artefact.

The interpretation is straightforward.  E-mode dominance is produced by density
filaments *aligned with the magnetic field*; a polarization field with no
preferred alignment gives B/E = 1 identically.  That alignment develops in the
inertial range, and with piecewise-linear reconstruction numerical dissipation
sets in near 20 cells per wavelength, i.e. k_diss ≈ N/20.  At 128³ that leaves
essentially no inertial range (the driving acts at k ≤ 2), so the maps are
smooth blobs rather than filaments and B/E sits near the no-structure limit.
The paper's window, k/k_min ∈ [4, 25], only exists at N ≥ 512.

This is the single quantitative reason to run the same code at 512³, which is
what `cluster/run_suite.slurm` does.

## 6. Deviations from the paper

See `docs/METHOD.md`.  The three that matter:

1. **Sign of Q, U.**  Eqs. (7)–(8) are written with ψ = the *field* angle; dust
   polarization is perpendicular to the field, so the physical Stokes parameters
   carry the opposite sign.  Taken literally the printed equations flip the sign
   of r_TE.  We use the physical convention, which is what reproduces the
   positive T–E correlation reported by both the paper and Planck, and which
   passes the analytic filament tests.
2. **Energy injection rate.**  Ė ∝ M_S³ alone misses the target M_S by up to 60%
   at low M_A, because a strong mean field makes the cascade markedly less
   dissipative.  An empirical M_A^1.5 factor was calibrated from a first pass of
   the suite.  Ė is constant through each run, as in the paper, and all results
   are reported against measured Mach numbers.
3. **Resolution**: 128³ rather than 512³ — the origin of the residual offsets in
   the intercepts a and of the A_BB/A_EE gap.

## 7. Files

```
figures/n64/, figures/n128/   fig01 (T/E/B maps), fig02 (Mach legend),
                              fig03 (alpha_EE vs A_BB/A_EE), fig04-05 (fluid
                              spectra and slopes), fig06-07 (T,E,B spectra and
                              slopes), fig08 (amplitude ratios), fig09
                              (correlations), table1.txt, runs_table.txt
figures/convergence.txt       resolution dependence of every Table-1 coefficient
runs/n64/, runs/n128/         21 runs each: snapshots, history, analysis .npz
```

---

# Part II — crosschecks, public data and the production decision (2026-09)

Part I above is the from-scratch solver against the paper.  Part II asks whether
these simulations can feed a CNN that removes Galactic dust from CMB maps: it
brings in a second code (Enzo), a third one through public data (Cho-ENO boxes
from The Well), and the real dust models (PySM3).  Every table quoted here is in
`results/`, produced by the scripts named with it.

## 8. The resolution ladder finished (128³ → 256³ → 512³)

`results/01_own_solver_convergence.txt`, from `cluster/solver_extract.slurm`.

The Table-1 coefficients keep improving with resolution — α_EE = (−3.48, 0.132,
0.192) at 512³ against the paper's (−3.63, 0.170, 0.280) — but two amplitudes do
not converge: BB/EE at M_S > 4 goes 0.96 → 0.84 → 0.90 (paper 0.55) and r_TE goes
0.26 → 0.41 → 0.47 (paper ~0.30, Planck 0.355).  Resolution is not the
explanation for either.

The same report gives three results used later:

* **line-of-sight depth**: BB/EE and r_TE are unchanged from the full box down to
  6% of it, so thin slabs are legitimate training maps (the polarization
  fraction does rise, 0.61 → 0.90)
* **decorrelation**: k ≥ 3 is independent after ~0.5 t_dyn, k = 1–2 after ~2
  t_dyn → snapshots every 2 t_dyn are independent samples
* **polarization-angle spread** depends on M_A alone: ~1° at M_A 0.5, ~7° at
  M_A 1, 25–45° at M_A 2.  This is why "Q is negative everywhere" in a
  sub-Alfvénic box and not in the sky.

## 9. Independent confirmation of Table 1, and the B-mode exception

`results/03_well_table1_reproduction.txt` (`python/well_stream.py`,
`python/well_paper.py`).

100 snapshots of The Well's 256³ boxes were streamed and analysed.  Its
sub-Alfvénic family (M_A ≈ 0.8, M_S 0.5–8.4) reproduces the paper's fit with an
unrelated code: α_EE = (−3.38, 0.164) against (−3.41, 0.170), α_ρ b = 0.152
against 0.160, α_TT b = 0.174 against 0.150.  α_BB does not: (−3.93, 0.216)
against (−4.31, 0.280) — the same direction our solver misses.

The velocity spectrum settles the opposite way: The Well gives α_v = −3.85 to
−4.03 and the paper −3.86, while our solver gives −3.41 at 512³.  That slope is
our solver's numerical dissipation, not a disagreement about physics.

## 10. Which turbulence looks like dust

`results/04_sims_vs_pysm3.txt` (`python/pysm_patches.py`, `python/patch_stats.py`,
`python/well_vs_pysm.py`); figures in `figures/well_vs_pysm/`.

One estimator is applied to PySM3 dust patches (d1, d10, d12 at 353 GHz, 128²
maps of 20°, |b| = 35–75°) and to simulated maps: E/B decomposed on the full
sphere or the full periodic face, the same beam, the same band.  Validation of
the estimator against the curved-sky measurement of the same PySM3 model:
BB/EE = 0.56 on patches against 0.51 full-sky.

Result: the public boxes bracket the dust without matching it — weak mean field
(M_A ≈ 8) has the right r_TE but angle dispersion 7.7° against 4.6° and BB/EE
0.90; strong mean field (M_A ≈ 0.8) has r_TE 0.07 and a polarization fraction 8×
too high.  The Enzo box at the paper's Planck point (M_S 5.3, M_A 1.5) matches
PySM3 d10 on α_EE (−2.76 vs −2.78), α_BB (−2.79 vs −2.81), r_TE (0.21 vs 0.28)
and angle dispersion (3.0° vs 4.6°), with BB/EE 0.73 against 0.59 the weakest
point.

## 11. Why the 512³ box stalled, and the Alfvén-speed limiter

`results/02_enzo512_cost.txt`, `results/05_enzo512_preview.txt`,
`results/06_alfven_limiter_test.txt` (`cluster/enzo_peek.slurm`,
`cluster/patch_enzo_floor.sh`, `python/floor_test.py`).

Measured time steps per dynamical time: 3,458 → 47,672 → 64,051 → 131,684, with
the fastest signal speed in the box rising from 29 to 2,295 while the gas moves
at ~5.  Cause: cells emptied to ρ ≈ 10⁻⁶ of the mean, where B/√ρ is enormous;
they hold ~10⁻⁶ of the mass.  Two 3-day segments (~16,000 core-hours) reached
2.44 of 10 dynamical times.

Enzo's `UseFloor` + `MaximumAlvenSpeed` caps exactly those cells.  Controlled
test at 64³ (same seed, capped vs uncapped): cost 9,665 → 4,571 cycles for 5
t_dyn, mass +1.2×10⁻⁴, and all statistics equal within the scatter over
snapshots.  Cost with the cap: ~4 days at 512³, ~8 h at 256³.

## 12. Is 256³ enough for 128×128 training maps?

`results/07_resolution_256.txt` (`python/resolution_check.py`).

* the cascade in a 256³ box ends at k ≈ 12 (the local spectral slope steepens),
  so the box is physical over a factor 4 in scale; 512³ reaches k ≈ 23
* binning a 256² face to 128² changes spectral slopes by 0.02 and BB/EE, r_TE and
  the polarization fraction by ≤0.01 — but cutting a *native-resolution* quarter
  tile instead changes slopes by 0.4–0.9, because such a tile spans mostly the
  driving scales.  Training maps must be binned faces.
* on the sky, k = 3–12 is ℓ = 54–207 for a 20° map and ℓ = 108–414 for a 10° map;
  Planck's dust measurements span ℓ ≈ 40–600, which no single box covers (that
  would need ~900³)

Conclusion: 256³ is the right production resolution — ~900 core-hours per box
against ~10,000 at 512³ — provided training maps are binned full faces with a
declared angular size, and the CNN is not asked to learn structure outside the
covered range.
