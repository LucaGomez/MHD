# Decision log

What was tried, what failed, and what each failure changed. Written for readers
who want to judge the research process, not only the result. Every number here
is measured in this repository; nothing is quoted from the paper except as a
comparison target.

---

## 1. Write the solver first, then distrust it

A from-scratch isothermal MHD solver (`src/mhd.c`) was validated before any
science: second-order convergence on a circularly polarized Alfvén wave (measured
order 1.97), an isothermal shock tube against the **exact** Riemann solution
(`python/exact_iso_tube.py`), ∇·B control through fully developed supersonic
turbulence, and five analytic filament configurations for the E/B pipeline
(`python/test_eb.py`).

That last test mattered later: it fixed the sign convention (polarization
perpendicular to the projected field) that decides the sign of r_TE. Taking the
paper's equations literally flips Q, U, E, B and r_TE.

**Why it mattered:** when the solver later disagreed with the paper, the
validation suite made it possible to say *which* parts of the disagreement could
not be a coding error.

## 2. The resolution ladder, and a result that refused to converge

21 runs at 128³, 256³ and 512³ (`cluster/run_suite.slurm`). The slopes converged
toward the paper's Table 1, but two amplitudes did not:

| | 128³ | 256³ | 512³ | paper |
|---|---|---|---|---|
| BB/EE at M_S > 4 | 0.96 | 0.84 | 0.90 | 0.55 |
| r_TE at M_S > 4 | 0.26 | 0.41 | 0.47 | 0.30 |

More resolution did not help. At that point there were three candidate
explanations — our solver, our analysis, or the paper — and no way to separate
them from inside the project.

**Decision:** stop adding resolution; get an independent code.

## 3. A bookkeeping failure that looked like three physics failures

The health report showed run `ms5_ma1` broken at 128³, 256³ *and* 512³. It looked
like an instability that survived resolution. It was not: the 128³ run blew up,
and the higher resolutions were started by upsampling its last snapshot, so NaNs
were copied up the ladder.

**Changes made:** analysis tools now skip non-finite snapshots
(`analysis.finite_snapshot`), the batch report filters runs whose Mach number is
not a finite number, and the fix was tested by deliberately corrupting a copy of
a good run.

**Lesson recorded here on purpose:** an agent that reports "three independent
failures" when the cause is one chained dependency is wrong in a way that costs
days of compute.

## 4. Getting Enzo to run at all

The paper uses Enzo, so the crosscheck had to. Four failures worth naming, each
found from the job logs:

| symptom | cause | fix |
|---|---|---|
| "MPI processes: 1" on 8 ranks | OpenMPI 5 needs PMIx | `srun --mpi=pmix` in `cluster/mpi_launch.sh`, plus a guard that refuses to resubmit if Enzo reports the wrong rank count |
| `GLIBCXX_3.4.32 not found` | system libraries ahead of the compiler's | `enzo_env.sh` keeps only compiler/MPI library paths |
| implicit-declaration errors | GCC ≥14 makes them errors | build C as `gnu11 -fpermissive` |
| post-processing died silently | `set -u` + conda hooks exit 127 | save shell flags from `$-`, restore after sourcing |

A first version of the MPI claim was wrong: the run had been serial, and the
"MPI reproduces serial" statement was retracted once the rank count was checked
explicitly. The guard exists so that cannot repeat silently.

## 5. Public data instead of more simulation — and why the labels lie

Before burning more compute, the obvious question: does the data already exist?
A search found [The Well](https://polymathic-ai.org/the_well/datasets/MHD_256/)
(10 parameter pairs × 10 trajectories × 100 snapshots at 256³, CC BY 4.0) and the
CATS catalogue behind it.

Rather than trusting the dataset card, the parameters were **measured** from the
data — and they differ from the labels:

- the files labelled `Ma = 2` have a mean-field Alfvén Mach number of **≈ 8**;
  the label uses the *total local* field, the paper uses the mean field
- the files labelled `Ms = 7` measure M_S ≈ 7.7–8.5 rms
- the sound speed is not stored at all; it was recovered by matching each run's
  mean |v| to the CATS pressure list, and the two `Ms = 0.5` files needed a
  pressure outside the documented list
- decorrelation was measured too: about 20 snapshots, so each trajectory holds
  ~5 independent states, not 100

**Consequence:** the public data bracket the dust but never sit at M_A ≈ 1.5,
which is the value the dust actually picks (§2.2 of the README). Had the labels
been trusted, the project would have concluded that public data covered the
Planck point. They do not.

**Also:** 376 GB per training file made downloading unattractive, so
`python/well_stream.py` reads single snapshots over HTTP range requests
(~4 min per snapshot) and analyses them in place. 100 snapshots were analysed
on a laptop with 13 GB of RAM.

## 6. Building an estimator that can be trusted on both sides

Comparing simulations with PySM3 dust required one estimator that is fair to
both. Three traps, each caught by a check rather than by reasoning:

1. **E→B leakage.** Cutting a small tapered patch leaks E into B: PySM3 d1
   measured BB/EE = 0.71 on patches against 0.51 on the full sky. Fix: decompose
   E/B on the full sphere (or the full periodic box face) and *then* cut.
   Patch value became 0.56 against 0.51.
2. **The beam.** Smoothing steepens fitted slopes; it is now divided out of every
   spectrum (`patch_stats.beam_sigma`).
3. **Conventions.** HEALPix (COSMO) Q/U, the IAU convention and our flat-sky
   convention differ by signs that change r_TE and can swap E and B. The
   orientation is derived per pixel from the true longitude/latitude grid, and
   the result is checked against the curved-sky measurement of the same model.

The validated estimator (`python/patch_stats.py`) is then applied identically to
PySM3 patches, Well boxes and Enzo boxes. This is what makes the table in §2.2 of
the README a like-for-like comparison.

## 7. The 512³ box: diagnosing a stall instead of buying more compute

The production box slowed to 131,684 time steps per dynamical time. The
temptation is to ask for a bigger allocation. Instead the time-step history was
reconstructed from the log and traced to its cause: a handful of cells emptied to
ρ ≈ 10⁻⁶, where the magnetic signal speed B/√ρ reaches ~2,300 while the gas moves
at ~5. Those cells carry ~10⁻⁶ of the mass.

Enzo's own Alfvén-speed limiter fixes exactly this, and reading the source showed
two things the documentation does not: `SmallRho` is *not* a general floor for
this solver (it only acts when a density goes negative), and the limiter prints
one line per capped cell per step — gigabytes of log at 256³. Hence the one-line
patch in `cluster/patch_enzo_floor.sh`.

**The limiter was then tested, not assumed** (`python/floor_test.py`): two 64³
runs, identical seed, capped and uncapped. Cost halved; mass changed by 1.2×10⁻⁴;
every statistic agreed within the snapshot-to-snapshot scatter. The runs
decorrelate after the limiter first fires, so the honest statement is
"consistent within the scatter", not "identical" — and the scatter at 64³ cannot
exclude differences smaller than ~0.3 in a spectral slope.

## 8. Choosing the production resolution by measurement

The final decision — 256³ or 512³ — was made by measuring where each box stops
being turbulent (`python/resolution_check.py`): the local spectral slope departs
from its inertial value at k ≈ 12 at 256³ and k ≈ 23 at 512³. Combined with the
map test (binning 256² → 128² changes slopes by 0.02) this gives a quantitative
answer: 256³ at ~900 core-hours per box, not 512³ at ~10,000.

## 9. Working notes on the process itself

- **Numbers over adjectives.** Every claim in the README has a file in
  `results/` produced by a script in `python/`.
- **Measure metadata.** Dataset labels, parameter names and documentation were
  wrong or ambiguous three times (The Well's M_A, its sound speeds, Enzo's
  `SmallRho`). Each was caught by measuring instead of reading.
- **Cheap experiments before expensive ones.** The limiter was validated at 64³
  on a laptop before touching a 112-core node; the estimator was validated
  against a curved-sky calculation before being applied to 100 snapshots.
- **Retractions are part of the record.** The "MPI reproduces serial" claim (§4)
  and the "three independent failures" report (§3) were both wrong and are kept
  here rather than quietly deleted.
- **Long-running work is checkpointed and self-resubmitting**, because a 3-day
  queue limit against a multi-day run is otherwise a silent data loss.
