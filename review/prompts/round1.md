You are reviewing a computational astrophysics project written by another AI
agent working with a PhD student. Be a demanding referee: the student presents
this tomorrow and wants to know what is wrong before the audience does.

WHAT THE PROJECT CLAIMS (all in README.md, with numbers in results/):
1. A from-scratch isothermal MHD solver (src/mhd.c) plus public 256^3 boxes
   reproduce Table 1 of Stalpes, Collins & Huffenberger 2024 (arXiv:2404.02874),
   except the B-mode slope alpha_BB, which two independent codes miss the same way.
2. Simulated dust maps match PySM3 dust models only at Alfvenic Mach number
   M_A ~ 1.5, which no public dataset provides.
3. A 512^3 Enzo run stalled because a few near-empty cells drive the Alfven speed
   to ~2300; Enzo's UseFloor/MaximumAlvenSpeed limiter fixes it and was shown not
   to change the physics (results/06).
4. 256^3 is sufficient for 128x128 CNN training maps: the cascade runs to k~12,
   binning a 256^2 face to 128^2 costs 0.02 in spectral slope, cropping native
   tiles costs 0.4-0.9 (results/07).

WHAT TO CHECK, in this order of importance:

A. Does the analysis code do what the text says it does?
   - python/analysis.py: projection to Stokes T,Q,U; the flat-sky E/B
     decomposition; shell averaging; the power-law fits. Check conventions and
     signs, normalization of the spectra, and whether the fit windows quoted in
     the README match the code.
   - python/patch_stats.py: the estimator used on both simulations and PySM3.
     Check the taper, the beam deconvolution, the band limits, and whether the
     polarization-angle dispersion S and its p-dependence are computed as
     described. Does supplying E,B from a larger domain really avoid leakage?
   - python/pysm_patches.py: the gnomonic projection, the COSMO/IAU sign
     handling, the per-pixel rotation into the patch frame, and the full-sky
     reference used to validate it.
   - python/well_stream.py: the sound-speed inference (PRESSURES list and the
     rescaling of M_S in load_state), the axis mapping checked by divb_check,
     and whether M_A is measured against the mean field as the paper defines it.

B. Are the conclusions supported by the numbers actually produced?
   Read results/*.txt and check each headline claim in README.md sections 2-4
   against them. Flag any number in the README or report/index.html that does not
   appear in results/, any claim stronger than the data supports, and any place
   where scatter or sample size is ignored (for example: 7 snapshots at 64^3 in
   results/06, one snapshot for the Enzo 512^3 row of the comparison table).

C. Statistical validity of the comparisons.
   - E->B leakage, beam effects and the k-band choice in the PySM3 comparison
   - whether "agrees within the scatter" in results/06 is a fair statement given
     the sample size, and what difference the test could actually exclude
   - whether k_diss ~ 12 (results/07) is robust to the criterion used
     (python/resolution_check.py, function knee)

D. Reproducibility and correctness of the pipeline scripts.
   Would someone else get these numbers from this repository? Check the commands
   in README.md section 5 against the actual argument parsers, and look for stale
   or broken references.

E. The physics.
   Flag anything questionable in docs/METHOD.md and src/mhd.c: the driving, the
   isothermal HLLD solver, divergence cleaning, and the dust-emission model
   (optically thin, grains perpendicular to B, no p0 factor applied).

HOW TO WORK
- You are in a read-only sandbox. Read files, run python/grep/git as needed, but
  do not modify anything.
- Prefer checking a number over reasoning about it: recompute from results/ where
  you can, and say so when you cannot.
- Do not summarize the project back. Only findings.

OUTPUT FORMAT — plain markdown, no preamble:

## Findings
For each, in this shape:

### [ID] short title
- severity: blocking | major | minor | question
- where: file:line (or results/<file>, README section)
- claim: what is wrong or unverified, in one or two sentences
- evidence: the command you ran, the number you read, or the code path
- suggested check: the smallest thing that would settle it

Rank blocking first. Aim for the 10 findings that matter most; do not pad.

## What I could not verify
List anything you could not check and why (missing data, needs a cluster, needs
a long run).

## Overall
Three sentences: is the headline claim of each of the four numbered claims above
supported, overstated, or unsupported?
