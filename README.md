# MHD turbulence → dust polarization maps for CMB foreground removal

Can magnetohydrodynamic (MHD) turbulence simulations produce Galactic dust maps
realistic enough to train a CNN that removes the dust foreground from CMB
polarization data?

This repository answers the upstream half of that question. It contains a
from-scratch isothermal MHD solver, an Enzo-based production setup, and an
analysis pipeline that measures dust polarization statistics the way Planck and
the reference paper do — plus the measurements that decide which simulations are
worth generating.

**Reference paper:** Stalpes, Collins & Huffenberger (2024), *Planck dust
polarization power spectra are consistent with strongly supersonic turbulence*,
[arXiv:2404.02874](https://arxiv.org/abs/2404.02874).

**Status:** the simulation and analysis chain is complete and validated. The CNN
itself is **not trained yet** — the training set needs simulations we could not
finish within the available compute (see [The 512³ story](#the-5123-story-why-a-single-box-cost-more-than-the-whole-project)).

---

## 1. The scientific chain

```
MHD turbulence (rho, v, B on a periodic box)
        |  project along a line of sight, grains aligned perpendicular to B
        v
Stokes maps T, Q, U  ->  E and B modes  ->  power spectra, BB/EE, r_TE, ...
        |  compare with Planck / PySM3 dust
        v
128x128 training maps  ->  CNN foreground removal   [not reached: compute]
```

Dust emission is optically thin, so for a line of sight along `z`

```
T = ∫ rho dz ,   Q = -∫ rho (Hx² - Hy²)/H² dz ,   U = -∫ rho 2 Hx Hy / H² dz
```

with the sign convention that makes the polarization perpendicular to the
sky-projected field. `python/analysis.py` implements this, the flat-sky E/B
decomposition, the shell-averaged spectra and the power-law fits.

## 2. What was measured

| # | Question | Answer | Where |
|---|---|---|---|
| 1 | Does our own solver reproduce the paper's Table 1? | Slopes yes, B-mode amplitude no | [results/01](results/01_own_solver_convergence.txt) |
| 2 | Do public MHD data reproduce it? | Yes, with the same B-mode exception | [results/03](results/03_well_table1_reproduction.txt) |
| 3 | Do simulated maps look like real dust (PySM3)? | Only at M_A ≈ 1.5 | [results/04](results/04_sims_vs_pysm3.txt) |
| 4 | Why did the 512³ Enzo box stall? | Alfvén speed in empty cells | [results/02](results/02_enzo512_cost.txt) |
| 5 | Does the fix change the physics? | No, within the scatter | [results/06](results/06_alfven_limiter_test.txt) |
| 6 | Is 256³ enough for 128×128 CNN maps? | Yes for the statistics, with a stated scale limit | [results/07](results/07_resolution_256.txt) |

### 2.1 Table 1 of the paper, reproduced with an independent code

The paper fits each spectral slope as `α = a + b·M_S + c·M_A` over 21 runs. We
reproduced that fit from [The Well](https://polymathic-ai.org/the_well/datasets/MHD_256/)
(Burkhart's CATS boxes, Cho–Lazarian ENO code — unrelated to the paper's Enzo),
using its sub-Alfvénic family (M_A ≈ 0.8, M_S from 0.5 to 8.4, 10 snapshots each):

| slope | this work: a, b | paper: a, b |
|---|---|---|
| α_ρ | −3.47, 0.152 | −3.61, 0.160 |
| α_v | −3.85, 0.010 | −3.75, 0.020 |
| α_H | −3.57, 0.031 | −3.53, 0.020 |
| α_TT | −3.64, 0.174 | −3.59, 0.150 |
| α_EE | −3.38, 0.164 | −3.41, 0.170 |
| **α_BB** | **−3.93, 0.216** | **−4.31, 0.280** |

Everything matches except B-modes, which come out shallower with a weaker M_S
dependence. **Our own solver misses in exactly the same direction**, so two
independent codes disagree with the published B-mode fit in the same way. That
points at the B-mode measurement or definition rather than at either code.

The same comparison exonerates the paper on the velocity spectrum and indicts
our solver: The Well gives α_v = −3.85 to −4.03, the paper −3.86, our solver
−3.41 at 512³ — our solver damps small-scale velocity too strongly.

### 2.2 Which turbulence actually looks like dust

Identical estimator on simulated maps and on PySM3 dust (d1, d10, d12 at
353 GHz, 128×128 patches of 20°, 35° ≤ |b| ≤ 75°, E/B decomposed on the full
sphere, same beam, band k = 3–13):

| source | α_EE | α_BB | BB/EE | r_TE | angle spread S | polarization fraction |
|---|---|---|---|---|---|---|
| PySM3 d1 | −2.72 | −2.60 | 0.56 | 0.17 | 5.6° | 0.088 |
| PySM3 d10 | −2.78 | −2.81 | 0.59 | 0.28 | 4.6° | 0.035 |
| PySM3 d12 | −2.78 | −2.99 | 0.60 | 0.30 | 3.3° | 0.057 |
| **Enzo, M_S 5.3, M_A 1.5** | **−2.76** | **−2.79** | **0.73** | **0.21** | **3.0°** | 0.365·p₀ |
| The Well, M_S 7, M_A ≈ 8 | −2.39 | −2.53 | 0.90 | 0.25 | 7.7° | 0.185·p₀ |
| The Well, M_S 7, M_A ≈ 0.8 | −2.61 | −2.60 | 0.98 | 0.07 | 1.6° | 0.746·p₀ |
| The Well, M_S 2.5, M_A ≈ 8 | −2.85 | −2.97 | 0.81 | 0.26 | 6.4° | 0.162·p₀ |
| The Well, M_S 2.5, M_A ≈ 0.8 | −3.48 | −3.64 | 0.99 | 0.09 | 1.0° | 0.728·p₀ |

![statistics](figures/well_vs_pysm/stats.png)

**The public data bracket dust but never match it.** Weak mean field (M_A ≈ 8)
gives the right T–E correlation but polarization angles that wander twice as
much as dust and BB/EE near 1. Strong mean field (M_A ≈ 0.8) gives almost no
T–E correlation, angles far too ordered, and a polarization fraction about eight
times too high.

**M_A ≈ 1.5 — the paper's Planck point, and the value missing from every public
dataset — lands on top of PySM3 on nearly every statistic.** That is the single
most important result here: it says the training maps have to be simulated by us,
at that field strength, and it tells us exactly which parameters to use.

Side-by-side maps (PySM3 d10 vs a simulated face):
![maps](figures/well_vs_pysm/maps.png)

## 3. The 512³ story: why a single box cost more than the whole project

We ran Enzo at 512³ at the paper's Planck point (M_S 4.7, M_A 1.5) on one
112-core node. Parallel efficiency was fine (2.8×10⁵ cell-updates per core per
second, 62% of the time in the MHD solver, load imbalance 1.04). The run still
stalled:

| t / t_dyn | time steps per t_dyn | fastest signal speed in the box |
|---|---|---|
| 0 → 1 | 3,458 | ~29 |
| 1 → 2 | 47,672 | ~250 |
| 2 → 2.13 | 64,051 | ~900 |
| → 2.44 | 131,684 | ~2,295 |

The gas itself moves at about 5. The runaway comes from a handful of cells that
the turbulence empties to ρ ≈ 10⁻⁶ of the mean: their magnetic signal speed
B/√ρ sets the time step of the whole box. Those cells hold ~10⁻⁶ of the mass and
emit no measurable dust signal.

![cost](figures/report/fig_cost.png)

Two 3-day segments (about 16,000 core-hours) advanced the box to 2.44 of the 10
requested dynamical times. Finishing without a fix would have cost 40,000–70,000
core-hours and weeks of queue at FairShare 0.

**The fix:** Enzo's Alfvén-speed limiter (`UseFloor = 1`,
`MaximumAlvenSpeed = 50`) raises the density of exactly those cells until their
signal speed is at the cap. Stock Enzo prints one line per capped cell per step,
which at 256³ means gigabytes of log, so `cluster/patch_enzo_floor.sh` silences
that message and rebuilds.

**Validation** (two 64³ runs, identical parameters and forcing seed, one capped,
[results/06](results/06_alfven_limiter_test.txt)):

- cost: 9,665 → 4,571 cycles for 5 t_dyn, and the capped run's cost per t_dyn
  stays flat (~1,030) while the uncapped one keeps drifting
- mass added: 1.2×10⁻⁴ of the total
- maps identical to four decimals until the limiter first fires; afterwards the
  runs are different realizations (turbulence is chaotic), so the comparison is
  statistical
- every statistic agrees within the scatter over snapshots: M_S 4.66 ± 0.25 vs
  4.80 ± 0.26, σ(ln ρ) 1.31 ± 0.05 vs 1.34 ± 0.09, α_EE −2.38 ± 0.40 vs
  −2.02 ± 0.35, BB/EE 1.11 ± 0.23 vs 0.89 ± 0.24, α_v −3.90 ± 0.28 vs −4.01 ± 0.18

With the cap, 10 t_dyn costs about 4 days at 512³ (~10,000 core-hours) and about
8 hours at 256³ (~900 core-hours).

## 4. Is 256³ enough to train the CNN?

Three measurements, all on real 256³ data ([results/07](results/07_resolution_256.txt)):

**(a) Where the cascade ends.** The local slope of the 3D spectra is flat
through the inertial range and steepens where the scheme's numerical dissipation
takes over. Measured on four Well runs: **k_diss ≈ 11–12**, so a 256³ box is
physical over k = 3–12, a factor of 4 in scale. A 512³ box reaches k ≈ 23.

![resolution](figures/report/fig_resolution.png)

**(b) Does a 128×128 map keep the statistics?** Binning the 256² projected face
to 128² changes the spectral slopes by **0.02** and BB/EE, r_TE and the
polarization fraction by ≤0.01 — the CNN's input resolution costs nothing. What
does matter is *which* 128² map you cut: a native-resolution quarter tile spans
only box k 1.5–6, i.e. mostly the driving scales, and its slopes are steeper by
0.4–0.9. **Bin the whole face; don't crop tiles.**

| product | α_TT | α_EE | α_BB | BB/EE | r_TE |
|---|---|---|---|---|---|
| native 256² face | −2.59 | −2.25 | −2.61 | 1.10 | 0.10 |
| 128² binned face | −2.60 | −2.26 | −2.63 | 1.10 | 0.09 |
| 128² native tile | −3.13 | −2.70 | −2.79 | 0.92 | 0.14 |

(The one exception is the polarization-angle dispersion S, which is measured at a
fixed pixel lag: after binning, the same lag is a larger angle, so S rises from
4.3° to 7.0°. Comparisons of S must fix the lag in degrees, not pixels.)

**(c) What that means on the sky.** A 128² map of angular size Θ carries k = 3–12
at ℓ = 360k/Θ:

| map size | ℓ at k = 3 | ℓ at k = 12 | pixel |
|---|---|---|---|
| 20° | 54 | 207 | 9.4′ |
| 10° | 108 | 414 | 4.7′ |
| 5° | 216 | 828 | 2.3′ |

Planck measures dust over ℓ ≈ 40–600, a factor 15 in scale. **No single box
covers that**: 256³ gives a factor 4, 512³ a factor 8, and the full Planck range
would need ~900³. 

**Conclusion.** 256³ is enough to train the CNN, provided that
1. training maps are binned full faces, not native tiles;
2. each map is assigned an angular size (10° is a good compromise: ℓ ≈ 108–414,
   inside Planck's measured range);
3. the CNN is not asked to learn structure outside that range from these maps —
   a wider range needs either several box sizes or maps stitched across scales.

Given that 256³ costs ~900 core-hours against ~10,000 for 512³, the right
production plan is **several 256³ seeds rather than one 512³ box**.

## 5. Reproducing this

Requirements: `numpy`, `scipy`, `matplotlib`, `h5py`, `astropy`, `healpy`,
`pysm3`, `fsspec` (for streaming The Well), a C compiler with OpenMP, and Enzo
for the production runs.

```bash
# our own solver: build and run a small case
bash cluster/setup.sh              # gcc -O3 -march=native -ffast-math -fno-finite-math-only
                                   #     -fopenmp -DUSE_FLOAT src/mhd.c -lm -o src/mhd
./src/mhd n=128 ms=5 ma=1 ntdyn=10 out=runs/n128/ms5_ma1
python3 python/process_run.py runs/n128/ms5_ma1 --axis 2
python3 python/figures.py --root runs/n128 --outdir figures/n128      # Table 1 + Figs 1-9

# public 256^3 data, analysed without downloading the files (~4 min per snapshot)
python3 python/well_stream.py --out runs_well --pairs Ma_0.7_Ms_7
python3 python/well_paper.py --out figures/well_paper                 # Table 1 reproduction

# dust reality check against PySM3
python3 python/pysm_patches.py --models d1 d10 d12 --out pysm_patches
python3 python/well_vs_pysm.py --pairs Ma_2_Ms_7 Ma_0.7_Ms_7 --enzo enzo_preview/maps_4x4.npz \
        --out figures/well_vs_pysm

# resolution study and report figures
python3 python/resolution_check.py --out figures/resolution
python3 python/report_figures.py --out figures/report

# Enzo production box on a SLURM cluster (FASRC syntax)
bash cluster/install_enzo.sh                        # build Enzo
bash cluster/patch_enzo_floor.sh                    # silence the limiter's per-cell message
M=/n/netscratch/.../mhd
SBATCH_TIMELIMIT=0-12:00 sbatch --export=ALL,MHD_ROOT=$M,N=256,VA_CAP=50 cluster/enzo_box.slurm
sbatch -p test -t 0-03:00 --export=ALL,MHD_ROOT=$M,TMIN_TDYN=2 cluster/enzo_peek.slurm   # early look
```

## 6. Layout

```
src/mhd.c                 isothermal ideal MHD: HLLD, PLM, Dedner cleaning,
                          SSP-RK2, Ornstein-Uhlenbeck driving, checkpoints
python/analysis.py        spectra, T/Q/U projection, E/B, power-law fits
python/process_run.py     per-run averages over snapshots -> analysis_ax*.npz
python/figures.py         the paper's Table 1 and figures
python/well_stream.py     stream The Well's 256^3 boxes over HTTP and analyse
python/well_paper.py      Table-1 fit from those boxes
python/pysm_patches.py    PySM3 dust -> flat 128^2 T/Q/U/E/B patches at 353 GHz
python/patch_stats.py     one estimator used on every map, simulated or not
python/well_vs_pysm.py    the comparison table and figure
python/resolution_check.py  is 256^3 enough for 128^2 maps?
python/floor_test.py      does the Alfven limiter change the turbulence?
python/extract_maps.py    cut training maps (T,Q,U,E,B) from a snapshot
python/enzo_setup.py      write an Enzo parameter file matched to the paper
python/enzo2snap.py       Enzo dumps -> our snapshot format
cluster/                  SLURM jobs: solver suite, Enzo box, previews, scaling
results/                  every measurement quoted above, as produced
figures/                  figures produced by the scripts
docs/METHOD.md            numerics and conventions
docs/RESULTS.md           detailed results log
docs/AGENT_LOG.md         decision log: what was tried, what failed, why
```

Large outputs (`runs/`, `runs_well/`, `pysm_patches/`, Enzo dumps) are not in
git; every script regenerates them.

## 7. What is left

1. **Finish one 256³ Enzo box at M_S 4.7, M_A 1.5 with the limiter** (~8 h on one
   node). Submitted as job 47520681; its statistics go through
   `python/well_vs_pysm.py` for the same PySM3 comparison.
2. **Produce the training set:** several seeds, snapshots every 2 t_dyn (measured
   decorrelation time), binned full faces, random rotations. About 1,600 maps of
   128×128 need roughly 5 seeds × 900 core-hours.
3. **Train the CNN** and test it against PySM3 dust models it never saw — the
   robustness question of [arXiv:2603.12364](https://arxiv.org/abs/2603.12364).
4. **Open question worth chasing:** the B-mode discrepancy in §2.1. Two
   independent codes disagree with the published α_BB fit in the same direction.

## 8. Data sources and credits

- **The Well**, MHD_256 (CC BY 4.0) — Ohana et al. 2024; boxes by B. Burkhart,
  [CATS](https://www.mhdturbulence.com), Burkhart et al. 2020, ApJ 905, 14.
- **PySM3** dust models d1, d10, d12 — Thorne et al. 2017, Zonca et al. 2021.
- **Enzo** — Bryan et al. 2014; stochastic forcing by Ph. Grete (ProblemType 59).
- **Reference paper** — Stalpes, Collins & Huffenberger 2024, arXiv:2404.02874.
- Planck dust values quoted from Planck 2018 XI and XII.

Cluster work ran on the FASRC Cannon cluster (Harvard).
