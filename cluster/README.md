# Running this on a cluster

## Dependencies

**Solver — none.**  `src/mhd.c` is C99 + OpenMP + libm.  No MPI, no HDF5, no
FFT library, no Fortran.  The only headers it uses are `stdio.h`, `stdlib.h`,
`math.h`, `string.h`, `sys/stat.h` and `omp.h`.

```bash
gcc -O3 -march=native -ffast-math -fno-finite-math-only -fopenmp -DUSE_FLOAT -o src/mhd src/mhd.c -lm
# icc/icx: icx -O3 -xHost -fp-model=precise -qopenmp -DUSE_FLOAT -o src/mhd src/mhd.c -lm
```

`-DUSE_FLOAT` stores the state in single precision (all arithmetic stays in
double); it is ~1.3x faster and halves the memory, and reproduces the
double-precision results to 5 significant figures on every validation test.
Drop it for a double-precision build.

Check the build before queueing anything:

```bash
OMP_NUM_THREADS=8 ./src/mhd test=cpaw      # expect 2nd-order convergence, ~1.97
OMP_NUM_THREADS=8 ./src/mhd test=divb      # |div B| dx/|B| stays ~1e-2
```

**Analysis — numpy, matplotlib, scipy.**  numpy does the FFTs; scipy is used
only by `python/exact_iso_tube.py` (the shock-tube validation).  Nothing else.

```bash
python3 -m pip install numpy scipy matplotlib
```

The analysis is cheap and can be run afterwards on a login node or locally.

## Resource requirements

Memory is 3 state arrays x 8 variables x (N+4)^3, plus 3 driving arrays of N^3
in double:

| N | memory (`-DUSE_FLOAT`) | memory (double) |
|---|---|---|
| 256 | 2.1 GB | 3.8 GB |
| 512 | 16.5 GB | 30 GB |
| 1024 | 130 GB | 240 GB |

One node per (M_S, M_A) pair; the code is shared-memory only, so use one task
with all the node's cores.  Thread scaling is sublinear beyond ~8 threads
(the scheme is memory-bandwidth bound), so on a large node it is usually
better to run 2-4 runs concurrently with a fraction of the cores each than one
run on all of them.

## Two ways to run the 512^3 suite

### (a) From scratch — simplest, most expensive

```bash
sbatch --export=N=512,NTDYN=12,NSNAP=19 cluster/run_suite.slurm
```

~30k steps per run.  Budget 24-48 h per run depending on the node.

### (b) Grid-refinement ladder — ~3x cheaper, recommended

Each resolution starts from the converged state of the previous one, so only
the new small scales have to fill in.  Nothing needs to be transferred: run the
cheap resolution on the cluster first.

```bash
# 128^3 from rest, 12 dynamical times   (~1 h/run on a 64-core node)
sbatch --export=N=128,NTDYN=12,TSNAP0=7.5,NSNAP=19 cluster/run_suite.slurm

# 256^3, 3 dynamical times, from the 128^3 states
sbatch --export=N=256,NTDYN=3,TSNAP0=1.0,NSNAP=11,UPSAMPLE=$SCRATCH/mhd/n128 \
       cluster/run_suite.slurm

# 512^3, 3 dynamical times, from the 256^3 states
sbatch --export=N=512,NTDYN=3,TSNAP0=1.0,NSNAP=11,UPSAMPLE=$SCRATCH/mhd/n256 \
       cluster/run_suite.slurm
```

If you would rather skip the 128^3 rung, the states computed here are in
`ic_n128/` (21 x 59 MB = 1.2 GB); copy that directory over and pass
`UPSAMPLE=ic_n128`.  The 128^3 results in `runs/n128/` are already analysed, so
this only saves cluster time, not information.

Runs checkpoint every `CHKINT` steps (default 500) and restart automatically
from `checkpoint.bin` if the job is requeued, so a 48 h wall limit is not a
problem — just resubmit.

## Analysis

```bash
python3 python/process_run.py $SCRATCH/mhd/n512/*/ --axis 2   # perpendicular
python3 python/process_run.py $SCRATCH/mhd/n512/*/ --axis 0   # along the field
python3 python/figures.py --root $SCRATCH/mhd/n512 --outdir figures/n512
python3 python/compare_res.py runs/n64 runs/n128 $SCRATCH/mhd/n256 $SCRATCH/mhd/n512
```

`process_run.py` writes one `analysis_ax2.npz` per run (a few MB).  If you want
to bring results back rather than the raw snapshots, those `.npz` files plus
`history.txt` and `params.txt` are all that `figures.py` and `compare_res.py`
need — about 40 MB for the whole suite instead of ~250 GB of snapshots.

The fit window scales with resolution as k/k_min in [3, N/20]; at N = 512 that
is [3, 26], matching the paper's [4, 25].

## What to expect

The residual differences from the published values at 128^3 are all consistent
with too little inertial range.  At 512^3 the quantities to watch are

* `A_BB/A_EE` for M_S > 4 — 0.93 at 128^3, paper 0.55, Planck 0.53
* `<r_TE>` for M_S > 4 — 0.22 at 128^3, paper ~0.3, Planck 0.355
* the intercepts a of the Table-1 fits, which are 0.1-0.2 too steep at 128^3

The b and c coefficients already agree at 128^3 and should stay put.
