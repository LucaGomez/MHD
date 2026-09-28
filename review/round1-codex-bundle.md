# Round 1 — Codex findings (bundle run)

Codex could not execute commands (its sandbox cannot create user namespaces on
this machine), so the repository was supplied on stdin as an inlined bundle:
README, docs/METHOD.md, results/*.txt and the nine main analysis sources with
line numbers. Verbatim output follows.

## Findings

### [F01] One perpendicular projection uses the wrong coordinate order for E/B
- severity: blocking
- where: python/analysis.py:112–119,125–130; python/well_vs_pysm.py:59
- claim: For `axis=1`, Stokes components use sky coordinates `(z,x)`, but the projected arrays retain NumPy axis order `(x,z)`. The E/B transform consequently combines polarization components and Fourier directions expressed in different frames, contaminating every comparison pooling axes 1 and 2.
- evidence: Summing an `(x,y,z)` array over axis 1 leaves `(x,z)`. The code instead selects `Hh=Hz, Hv=Hx` and returns the maps without transposing them. Relative to the array frame, this reverses Q while leaving U unchanged; it is not an innocuous overall polarization sign. Both the PySM comparison and resolution study include this axis.
- suggested check: Transpose all three projected maps for axis 1, or consistently express their Stokes components in the retained array frame. Test an oblique pure-E pattern under all three axis permutations, then regenerate results/04 and results/07.

### [F02] The crop-to-box wavenumber conversion is reversed
- severity: blocking
- where: python/well_vs_pysm.py:11–14; python/resolution_check.py:13–16,92–94; README.md:175–177
- claim: A half-width tile samples **twice**, not half, the box wavenumber at a given patch k. This invalidates the explanation of cropping bias and makes the Enzo comparison sample a substantially different physical band.
- evidence: Equating wavelengths gives `Lpatch/kpatch = Lbox/kbox`, hence `kbox = kpatch × Lbox/Lpatch`. Well 128/256 tiles therefore probe box k=6–26 for patch k=3–13; Enzo 128/512 tiles probe k=12–52. The latter extends well beyond the project's asserted 512³ dissipation scale of approximately 23. Moreover, results/07's aggregate crop shifts are −0.54, −0.45 and −0.18 for TT, EE and BB—not uniformly 0.4–0.9.
- suggested check: Compare binned full faces at matched box wavenumbers and angular resolution. Repeat the Enzo comparison below an independently measured dissipation scale before interpreting its steep slopes as a preferred Mach number.

### [F03] The limiter experiment does not establish unchanged physics
- severity: major
- where: results/06_alfven_limiter_test.txt:4–27; python/floor_test.py:95–108; README.md:145–159
- claim: Overlapping snapshot scatter is not an equivalence test, and seven dumps from one seed at 64³ cannot establish negligible bias at 256³ or 512³. The measured changes are substantial compared with the claimed scientific precision.
- evidence: Capping changes α_EE by **+0.367** and BB/EE by **−0.219**, approximately a 20% decrease in the latter. As an illustrative calculation from the printed numbers, assuming seven independent observations per run and correcting the code's population standard deviations gives approximate 95% Welch intervals of **[−0.10,+0.84]** and **[−0.52,+0.08]** for those differences. These are not valid final intervals without temporal and cross-run covariance; they demonstrate that even a simple independence model cannot exclude large effects. The reported mass addition is a difference of *time-averaged* densities, not cumulative injected mass. No capped production-resolution timing measurement is supplied.
- suggested check: Specify acceptable changes before testing; use multiple paired seeds, account for temporal correlation, and measure at production resolution. Record cumulative injected mass and changes to momentum/energy, and label production runtime estimates as extrapolations.

### [F04] The actual slope estimator includes driving and nominally dissipative scales
- severity: major
- where: python/process_run.py:56–86,88–111; docs/METHOD.md:68–71; results/01_own_solver_convergence.txt:36; results/03_well_table1_reproduction.txt:2
- claim: The quoted fit window is not the window defining the reported mean slope. Spectral ratios and correlation statistics also have multiple definitions across the pipeline, preventing an unqualified “same measurement” claim.
- evidence: At 256³, the printed window is [3,13], but slopes average twelve fits with lower bounds **2,3,4** and upper bounds **13,14,15,16**. At 512³ the corresponding union is [2,29], versus the paper's quoted [4,25]. This includes forcing scales and scales beyond the claimed cutoff. `process_run` averages per-shell correlation coefficients, whereas `patch_stats` correlates summed shell powers; the printed `B/E` is a fitted intercept ratio while `bp_BB_EE` is a band ratio. Projection also sums density without the cell-depth factor, making absolute projected spectral amplitudes resolution dependent.
- suggested check: Publish a fixed-window estimate separately from window sensitivity, and rerun with matched physical windows. Explicitly identify which ratio each table uses and include cell-depth normalization wherever absolute amplitudes are compared.

### [F05] Table 1 is not reproduced except for a B-mode discrepancy
- severity: major
- where: README.md:49–79; results/01_own_solver_convergence.txt:10–39; results/03_well_table1_reproduction.txt:6–25
- claim: The supplied numbers show several discrepancies, and neither the Well regression nor the shared analysis supports isolating the problem to B-mode measurement. The statement that both codes miss α_BB “in exactly the same direction” is not valid throughout the parameter range.
- evidence: At 512³, the solver's velocity intercept is −3.410 versus −3.860, magnetic intercept −2.959 versus −3.310, and magnetic Mach slope −0.023 versus +0.020. BB/EE is 0.895 versus 0.550. The Well fit uses only five parameter settings near M_A≈0.8, including M_S=8.36 outside the stated calibration range, so it cannot reproduce the three-parameter Table 1 fit. Its BB residual changes from +0.52 at M_S=2.53 to −0.28 at M_S=8.36. Both simulation codes use the same analysis pipeline. Finally, a shallower velocity spectrum alone does not establish “too much small-scale damping.”
- suggested check: Compare predictions and uncertainties over the overlapping Mach range, with matched windows and driving. Exclude the extrapolated point as a sensitivity test and validate the spectra with an independent analysis implementation.

### [F06] “Only M_A≈1.5 matches dust” is not supported by the comparison
- severity: major
- where: README.md:100–109; results/04_sims_vs_pysm3.txt:21–164
- claim: One Enzo snapshot partitioned into sixteen tiles cannot establish a preferred parameter region or exclude all public alternatives. The claimed polarization-fraction mismatch also compares quantities with an unspecified relative normalization.
- evidence: Tile counts are not independent simulation realizations, and multiple Well tiles and sightlines share the same state. Enzo's BB/EE interval [0.498,1.067] overlaps those of both PySM and public simulations. Statistics omitted from the README comparison are less favorable: Enzo's S–p slope is −1.059 [−1.188,−1.007], versus d12's −0.372 [−0.546,−0.204]; its median log-intensity skewness is 0.016 versus d10's 0.545 and d12's 0.612. The “eight times too high” argument compares **0.746·p₀** with **0.088**; their ratio is **8.48·p₀**, not 8.48. No public-dataset inventory establishes universal absence of M_A≈1.5.
- suggested check: Treat this as a candidate parameter point. Compare multiple equilibrated seeds over a local Mach-number grid using a joint statistic, account for shared snapshots and overlapping sky patches, and specify or marginalize over p₀.

### [F07] The claimed polarization-projection validation bypasses the Q/U transformation
- severity: major
- where: python/pysm_patches.py:49–64,84–107; python/well_vs_pysm.py:78–82; python/patch_stats.py:61–89
- claim: Positive TE and plausible BB/EE do not validate the COSMO/IAU conversion or patch-frame rotation: those statistics use separately projected full-sky E/B and never use the transformed Q/U. The reference and beam treatment also leave systematic errors uncalibrated.
- evidence: `pysm_stats` supplies E/B explicitly, bypassing `qu_to_eb`; an erroneous U sign or angle rotation could therefore survive this comparison while affecting S and its p-dependence. The rotation derives an angle from `gradient(latitude)` on a non-conformal gnomonic grid, without a demonstrated spin-basis transport test. `fullsky_reference` masks Q/U before `anafast` and applies no mode-coupling correction, so it is a masked pseudo-spectrum with possible E/B mixing. Supplying global E/B **does avoid the patch-boundary spin-decomposition leakage**, but scalar tapering still convolves power between wavenumbers; division by a Gaussian at each shell center does not undo that convolution. Enzo tiles are additionally smoothed with periodic wrapping *after* cropping, unlike the Well faces. S uses a fixed pixel lag and ordinary pixelwise log–log regression, not a fixed sky separation across the distorted projection.
- suggested check: Pass known full-sky pure-E and mixed-E/B realizations through the entire projection pipeline. Check transformed Q/U against an explicit transported basis, quantify window/beam transfer and B-sign consistency, and test S at matched angular separations.

### [F08] Sound speeds are fitted to labels rather than established from simulation metadata
- severity: major
- where: python/well_stream.py:38–49,74–88,113–124
- claim: The sonic Mach numbers driving the regression depend on an unverified pressure lookup inferred from the desired labels. Cached states can subsequently acquire revised Mach numbers without consistent updates to their stored velocity and magnetic spectral amplitudes.
- evidence: `pick_cs` chooses whichever pressure makes instantaneous mean speed agree best with the filename label; the list adds 2.0 and 3.2 beyond the documented pressures. It runs independently for every snapshot, although sound speed should be fixed within an isothermal run. `load_state` revises only M_S when the list changes, leaving `Cv`, `CH` and metadata in their previous units. The mean-field M_A formula itself is consistent with the stated definition and invariant under joint v/B rescaling; this does not validate the inferred c_s. `divb_check` merely compares two axis assignments and logs a ratio—no supplied result establishes that the mapping passes.
- suggested check: Obtain each run's sound speed from original parameters or authoritative metadata, fix it per run, and regenerate cached states. Publish measured divergence diagnostics and the coordinate metadata supporting the axis assignment.

### [F09] The dissipation threshold and 512³ reach are not independently measured
- severity: major
- where: python/resolution_check.py:39–47,75–87,127–138; results/07_resolution_256.txt:4–11; README.md:165–207
- claim: The knee finder cannot detect dissipation at k≤10, and the claimed 512³ result is simply twice the 256³ value. Preservation of low-band statistics after binning does not demonstrate that a CNN's full input is physically resolved.
- evidence: `knee` searches only `k>10`; five of the twelve reported knees equal its first possible answer, 11. The median of the listed knees is **11.5**, displayed as 12; this explains why the sky table gives 207 rather than 216 at 20°. The actual spectral cutoff is then k≤11.5, so integer shell 12 is excluded despite the printed “[3,12]”. Line 86 prints `2*kd_typ` as if obtained by applying the criterion to a 512³ box, but no 512³ data are read. The 0.02 binning result is a difference of pooled medians, not a bound on individual-map changes.
- suggested check: Vary the reference interval, slope stencil and tolerance; inspect individual-state knees and matched-resolution spectra. Report the unrounded cutoff and paired map changes, and validate the CNN on explicitly band-limited inputs before claiming training sufficiency.

### [F10] The documented recipe cannot regenerate the advertised tables from a clean checkout
- severity: major
- where: README.md:227–237,273–274; python/well_vs_pysm.py:120–137; results/04_sims_vs_pysm3.txt:25,38
- claim: The supplied commands omit required datasets and intermediate products, while the comparison script silently skips missing inputs. Some README values also disagree with their cited result file.
- evidence: The streaming command requests only `Ma_0.7_Ms_7`; results/03 requires all ten pairs, and results/04 and results/07 additionally require both M_S≈2.5 runs and the weak-field M_S≈7 run. No command creates `enzo_preview/maps_4x4.npz`. Missing PySM files, Well pairs and Enzo files are silently omitted. README's d12 slopes **−2.78,−2.99** do not match results/04's **−2.735,−2.964**, which round to **−2.74,−2.96**. The argument syntax of the included Python commands is valid; missing inputs and provenance are the demonstrated failures.
- suggested check: Provide a complete command sequence with explicit pair lists and Enzo extraction steps, require expected sources to exist, and generate README tables directly from versioned result summaries.

## What I could not verify

- No commands were executed, as instructed. Numerical differences and illustrative intervals above were calculated from the supplied text; spectra, transfer functions and sensitivity tests could not be recomputed without arrays.
- `src/mhd.c` was absent. I could not inspect HLLD wave states and degeneracies, positivity handling, Dedner source terms, timestep stability, or actual driving implementation. METHOD's claims of exactly zero momentum and helicity need implementation evidence: absence of a zero Fourier mode does not alone ensure zero density-weighted momentum injection or zero helicity in each realization.
- The Enzo template, limiter implementation, patch script, cluster logs and dumps were absent. I could not independently establish the stall's cell-level cause, limiter update semantics, initial map identity, or projected production costs.
- `report/index.html`, `well_paper.py`, `figures.py`, `report_figures.py`, `extract_maps.py`, `enzo2snap.py` and cluster scripts were absent. Their numbers, parsers, references and provenance could not be audited.
- The paper and original dataset metadata were not included. Exact published conventions, pressure assignments and driving equivalence remain unverified; comparisons above use the coefficients and descriptions supplied in the bundle.
- `process_run.py` reads every matching snapshot without a time cut. Whether the solver or cluster extraction already restricts files to the claimed second-half analysis window cannot be determined.
- The dust projection implements an idealized density-weighted, optically thin model. Its omitted intrinsic polarization factor, emissivity variations and alignment variations cannot be calibrated from these data, nor can downstream CNN performance be established without training and independent validation.

## Overall

Claim 1 is overstated: several quantities disagree, and shared analysis defects prevent attributing the residual specifically to B modes; claim 2 is unsupported as an exclusive inference of M_A≈1.5.  
Claim 3's timestep diagnosis and low-resolution acceleration are supported by the reported numbers, but unchanged physics and production-resolution performance are unproven.  
Claim 4 supports only preservation of selected aggregate statistics under binning: the crop-scale interpretation is wrong, the dissipation reach is criterion dependent, and CNN sufficiency remains unsupported.
