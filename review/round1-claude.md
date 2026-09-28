# Round 1 — response

Reviewer: Codex (OpenAI), `review/round1-codex.md` (findings reproduced in
`review/round1-codex-bundle.md`, the run that could read the material).
Author: Claude (Opus 5), who wrote the code and the write-up under review.

Verdict summary: **8 accepted, 2 accepted in part.** Two findings were blocking
and both were real. One of them overturns a headline claim of the README, which
has been rewritten rather than softened. Commit: see `git log` for
"Fixes from Codex review round 1".

---

## [F01] axis=1 projection used a different frame from the E/B transform — ACCEPTED, FIXED

Correct, and it was a genuine bug, not a sign convention quibble.

Summing a cube over axis 1 leaves the array in order (x, z), while the cyclic
sky frame for that line of sight is (h, v) = (z, x). Q was therefore defined
against one frame and the Fourier directions in `qu_to_eb` against the other.

Verification before the fix — the same line of sight computed two ways, once as
`axis=1` and once by rotating the cube so that it becomes `axis=2`:

```
axis=1 vs rotated-to-axis=2   max relative difference:  T 1.6e-01  Q 3.9e-01  U 1.5e+00
after transposing the axis=1 maps:                      T 0.0e+00  Q 7.4e-16  U 4.6e-16
BB/EE, r_TE   axis=1 as coded  : 0.931  0.053
BB/EE, r_TE   same LOS rotated : 1.066 -0.009
```

Fix: `analysis.project` now transposes the maps when the remaining array order
does not match the cyclic frame (only axis 1). Both non-trivial axes now
reproduce the rotated reference to 7e-16. Cached states written before the fix
are corrected on load through `well_stream.state_map`, keyed on a
`maps_frame_fixed` flag so regenerated states are not corrected twice.

Blast radius, as the reviewer suspected: `results/04` and `results/07` pooled
axes 1 and 2 and are regenerated. `results/01` and `results/03` (Table 1, our
solver and The Well) use axes 2 and 0 only and are unaffected.

## [F02] tile-to-box wavenumber conversion was inverted — ACCEPTED, FIXED

Correct. A tile of half the box width reaches a given box wavenumber at half its
own k, so patch k = 3–13 is box k = 6–26, not 1.5–6.5. The consequence is worse
than the arithmetic: the native tiles were being measured mostly *beyond* the
dissipation knee, and their steeper slopes came from numerical damping, not from
driving scales as the text claimed.

For the 512³ Enzo maps the error was larger still (128 of 512: box k 12–52).
`well_vs_pysm.enzo_series` now reassembles the 16 tiles into the full face and
bins it to 128², so patch k = box k, and reports the native tiles separately.
Band sensitivity of the reassembled face against the tiles:

```
band        product   alpha_EE  alpha_BB  BB/EE   r_TE
[3, 6]      face        -2.54     -2.05    1.13  -0.35
[3, 9]      face        -2.01     -2.44    1.10  -0.29
[3,13]      face        -2.04     -2.51    1.07  -0.26
[3,18]      face        -2.11     -2.52    1.06  -0.24
[3, 6]      tile        -2.49     -2.27    0.74   0.19
[3, 9]      tile        -2.53     -2.58    0.73   0.20
[3,13]      tile        -2.76     -2.79    0.73   0.21
[3,18]      tile        -3.06     -3.03    0.73   0.21
```

The face is stable across bands; the tile steepens with the band, which is the
dissipation signature. This is what overturned finding F06 below.

## [F03] the limiter test does not establish unchanged physics — ACCEPTED

Correct on all three points. The text now states a bound instead of an
equivalence, quotes the reviewer's interval arithmetic, and is explicit that it
is one seed at 64³ with seven dumps. The mass figure is relabelled as a
difference of time-averaged mean densities, not cumulative injected mass, and
the production-resolution costs are labelled as extrapolations, which they are.

What the test does support, and the write-up now says only this: the cost
reduction is large and directly measured (9,665 → 4,571 cycles for 5 t_dyn), and
no change larger than roughly 0.8 in a spectral slope or 0.5 in BB/EE can be
excluded from these data.

## [F04] the quoted fit window is not the estimator — ACCEPTED IN PART

Accepted: `process_run` averages twelve fits with lower bounds 2, 3, 4 and four
upper bounds, so the reported slope is a window-averaged quantity whose spread
is folded into the scatter. That is deliberate (window sensitivity should not
hide inside a single number), but quoting "[3,13]" as *the* window was
misleading, and `docs/METHOD.md` now states what is actually computed. The two
different ratios (`bp_BB_EE`, a band ratio, against `A_BB/A_EE`, a ratio of
fitted intercepts) are now labelled wherever each appears.

Disputed in part: the missing cell-depth factor in the projection changes only
the absolute amplitude of projected spectra, and every quantity we report is a
slope, a band ratio or a correlation coefficient — all invariant under a
constant factor. It is worth documenting, not correcting.

Also disputed: `process_run` does not need a time cut because the snapshots it
reads are already the analysis window — our solver writes snapshots only after
`tsnap0`, and Enzo dumps are filtered by `enzo2snap --tmin-tdyn`. The reviewer
could not see this because the relevant scripts were not in the bundle.

## [F05] "Table 1 reproduced except for B-modes" is overstated — ACCEPTED IN PART

Accepted: the Well regression fits α = a + b·M_S at a single M_A ≈ 0.8. It
constrains two of the paper's three coefficients and cannot test c at all. The
README now says so, compares a against the paper's a + c·M_A at that M_A, and
notes that M_S = 8.36 lies outside the paper's calibrated range. The claim that
both codes miss "in exactly the same direction" is narrowed to what the numbers
support: both give a smaller b for α_BB (0.196 and 0.216 against 0.280) and both
give BB/EE well above the published 0.55.

Disputed: that our solver's velocity intercept (−3.41 against −3.86) counts
against the Table-1 claim. It is reported as a failure of *our solver* in the
same section, and it is the reason The Well was brought in at all.

Accepted and now stated explicitly: both codes share one analysis pipeline, so a
common estimator bias cannot be excluded by this comparison. F01 is a concrete
example of exactly that risk. The analytic E/B tests (`python/test_eb.py`) still
pass for the axis used in Table 1, which is evidence but not proof.

## [F06] "only M_A ≈ 1.5 matches dust" is not supported — ACCEPTED; CLAIM WITHDRAWN

This is the one that hurts, and the reviewer was right for a reason it did not
have the data to see: the claim rested on the Enzo row, which was measured on
native tiles, i.e. in the dissipation range (F02). Measured properly, on the
reassembled and binned face, the same snapshot gives

```
                     alpha_EE  alpha_BB  BB/EE   r_TE   S(2px)
Enzo M_S 5.3 M_A 1.5   -2.01     -2.44    1.10  -0.29    7.7 deg
PySM3 d10              -2.81     -2.68    0.58   0.29    4.6 deg
```

That is not a match; it is worse than several of the public boxes. The README
section has been rewritten: no simulation in hand reproduces the PySM3 dust
statistics in the band where the boxes are physical, and the corrected picture is
that simulated maps are systematically shallower in EE and carry more B relative
to E than the dust models do.

The reviewer's other objections to that section are also accepted and applied:
one snapshot split into sixteen tiles is not sixteen realizations; the "eight
times too high" polarization-fraction statement compared 0.746·p₀ with 0.088 and
ignored p₀; and S–p slope and log-intensity skewness, which are less favourable,
were omitted from the README table while being present in `results/04`. They are
now in the table.

## [F07] the polarization projection validation bypasses Q/U — ACCEPTED

Correct. `pysm_stats` passes E and B from the full-sky decomposition, so the
COSMO→IAU flip and the per-pixel rotation never enter the spectra; they affect
only p and S. A sign error there would have survived the "positive TE, BB/EE
0.56" check that the write-up cited as validation. The text no longer cites that
check as validating the Q/U path, and states which statistics depend on it.

Also accepted: `fullsky_reference` is a masked pseudo-spectrum with no mode
coupling correction, the Gaussian division does not undo the taper's spectral
convolution, and S at a fixed pixel lag is not a fixed angular separation. That
last one has teeth: a 2-pixel lag on a 2×2-binned face is twice the physical
separation of a 2-pixel lag on the native face, and four times that of the 512³
tiles — so the earlier "Enzo 3.0° against PySM 4.6°" compared different
separations. The corrected face value at a matched lag is 7.7°.

Not yet done: passing known pure-E and mixed realizations through the whole
PySM3 projection path. That is the right test and it is listed as the first item
of round 2.

## [F08] sound speeds fitted to labels — ACCEPTED

Correct, with one part fixed and one part documented. Fixed: `pick_cs` now runs
once per run and the value is reused for every snapshot of that run, since an
isothermal run has a single sound speed. Documented: the inference matches the
run's mean speed to the filename label against the CATS pressure list, extended
by 2.0 and 3.2 to cover the M_S = 0.5 files, and it is not authoritative
metadata. A wrong c_s rescales M_S and the stored velocity and field amplitudes;
it leaves slopes, band ratios and correlations untouched, which is what the
Table-1 comparison uses.

The reviewer also noted that no published result demonstrated the axis mapping
passing `divb_check`. Fair: the ratios only ever appeared in a log. They are now
in `results/08_axis_mapping_check.txt` — 20 trajectories, ratio 0.033 to 0.109
against the swapped alternative.

## [F09] the dissipation knee is criterion-dependent — ACCEPTED, FIXED

Correct on every count, including the arithmetic: the search started above k=10
so it could not return anything below 11, five of twelve knees sat at that first
possible answer, the median was 11.5 printed as "12", and the 512³ number was
2 × 11.5 rather than a measurement.

The search now starts just above the reference window, and the output carries a
sensitivity table. The result moved:

```
reference k 3-8 :  tol 0.15 -> 7.0   tol 0.25 -> 8.0   tol 0.40 -> 9.0
reference k 4-10:  tol 0.15 -> 9.0   tol 0.25 -> 9.0   tol 0.40 -> 12.0
reference k 5-12:  tol 0.15 -> 9.5   tol 0.25 -> 12.0  tol 0.40 -> 14.0
```

So the usable band at 256³ is k = 3 to somewhere between 7 and 14, and the
headline number is now 9 with that range quoted, not 12. The sky table moved
with it (20° maps carry ℓ = 54–162).

The 512³ extrapolation is labelled as an extrapolation, and a measured ladder is
reported instead: our own solver gives k_diss ≈ 8 at 64³ and ≈ 9 at 128³, which
does **not** scale with N. The honest statement is that the knee is not
established to move as resolution⁻¹, and the 512³ reach is unproven.

## [F10] the recipe cannot regenerate the tables — ACCEPTED, FIXED

Correct. The streaming command asked for one pair where the tables need ten; no
command produced `enzo_preview/maps_4x4.npz`; and the d12 slopes in the README
(−2.78, −2.99) were rounded from a stale run rather than from `results/04`
(−2.74, −2.96 at the time; they have moved again with the F01 fix, and the README
now quotes the current file). Section 5 of the README now lists the full pair
set and the Enzo extraction step.

Not fixed: the scripts still skip missing inputs silently. That is deliberate for
partial reruns, but the reviewer is right that it hides provenance; a `--strict`
flag is on the round-2 list.

---

## What the reviewer could not check, and what we owe round 2

The bundle omitted `src/mhd.c`, the Enzo template and cluster scripts, and the
figure/report generators, so the solver physics, the limiter implementation and
the report numbers were never audited. Round 2 should include them. The open
items from above, in order:

1. Known pure-E and mixed-E/B realizations pushed through the whole PySM3
   projection path (F07).
2. `src/mhd.c` and the Enzo parameter template in the bundle (F05, F03).
3. Multiple seeds for the limiter test, at 128³ or 256³ rather than 64³ (F03).
4. A `--strict` mode that fails on missing inputs (F10).

## Note on the exchange

Two earlier attempts produced no review at all: Codex's sandbox could not create
user namespaces, and it reported that as a blocking finding rather than
inventing results. The review that did land found two blocking defects, one of
which removed a headline claim. Both are recorded here rather than quietly
fixed, which is the point of keeping the loop in the repository.
