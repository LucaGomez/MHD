# Send to the cluster and run

Replace `USER@CLUSTER` and the remote path with yours.

## 1. Copy the project (8.6 MB — code, scripts, docs, figures)

```bash
rsync -avz --progress \
  --exclude 'runs/' --exclude 'ic_n128/' --exclude '__pycache__/' \
  --exclude 'src/mhd' --exclude 'src/mhd_*' \
  /home/lgomez/MHD/ USER@CLUSTER:~/MHD/
```

## 2. Build and validate (one command, ~1 minute)

```bash
ssh USER@CLUSTER 'cd ~/MHD && bash cluster/setup.sh'
```

Expect: Alfvén-wave convergence order ≈1.97, |∇·B|Δx/|B| ~1e-2, and five E/B
unit tests passing.  If `gcc` is not the right compiler there, use e.g.
`CC=icx bash cluster/setup.sh`.

## 3. Submit the jobs

FASRC partitions (docs.rc.fas.harvard.edu/kb/running-jobs):

| partition | cores/node | RAM/node | max time | job limit |
|---|---|---|---|---|
| **sapphire** | 112 | 990 GB | 3 days | none |
| shared | 48 | 172 GB | 3 days | none |
| test | 112 | 990 GB | 12 h | **5 jobs / 112 cores** |

Each array index counts as one job, so a 21-task array does not fit in `test` —
that is the `QOSMaxSubmitJobPerUserLimit` error.  Use **sapphire** for the real
rungs.  Defaults in `run_suite.slurm` are now `-p sapphire -c 32 -t 1-00:00
--mem=32000`.

### Pick `-c` from measurement, not guesswork

The scheme is memory-bandwidth bound, so scaling flattens well before the core
count and `-c 112` mostly burns fairshare.  Measure it once:

```bash
salloc -p test -c 112 -t 0-01:00 --mem=64000
bash cluster/bench.sh 256          # ~2 min; prints steps/s vs thread count
exit
```

Take the smallest thread count still near the best rate.  (On my 8-core laptop
the curve is flat past 8 threads; a sapphire node with 8-channel DDR5 should
keep scaling considerably further, which is the number worth knowing.)

### The three rungs

Wait for each rung to finish before submitting the next — each starts from the
previous one's final state.

```bash
ssh USER@CLUSTER
cd ~/MHD

# smoke test, 3 runs (the M_S=4 triplet), fits inside the -p test limits
sbatch -p test --array=12,13,14 -c 16 -t 0-08:00 --mem=8000 \
       --export=ALL,N=128,NTDYN=12,TSNAP0=7.5,NSNAP=19 cluster/run_suite.slurm

# rung 1: 128^3 from rest, 12 dynamical times
sbatch -c 16 -t 0-08:00 --mem=8000 \
       --export=ALL,N=128,NTDYN=12,TSNAP0=7.5,NSNAP=19 cluster/run_suite.slurm

# rung 2: 256^3, 3 dynamical times, from the 128^3 states
sbatch -c 32 -t 1-00:00 --mem=16000 \
       --export=ALL,N=256,NTDYN=3,TSNAP0=1.0,NSNAP=11,UPSAMPLE=<root>/n128 \
       cluster/run_suite.slurm

# rung 3: 512^3, 3 dynamical times, from the 256^3 states  <-- paper's resolution
sbatch -c 64 -t 2-12:00 --mem=48000 \
       --export=ALL,N=512,NTDYN=3,TSNAP0=1.0,NSNAP=11,UPSAMPLE=<root>/n256 \
       cluster/run_suite.slurm
```

`<root>` is whatever `run_suite.slurm` printed as `output root:` (see below).

### If you are stuck on a partition with a low job cap

`NCHUNKS=n` makes each array task run its share of the 21 pairs sequentially,
so the whole suite fits in `n` jobs.  Striding by `n` gives each task one M_A
value across all M_S, so the tasks cost about the same:

```bash
# whole 128^3 suite in 3 jobs, 7 runs each
sbatch -p test --array=0-2 -c 16 -t 0-12:00 --mem=8000 \
       --export=ALL,NCHUNKS=3,N=128,NTDYN=12,TSNAP0=7.5,NSNAP=19 \
       cluster/run_suite.slurm
```

If a chunk hits the wall limit, resubmit the identical command: runs that are
already finished resume from their checkpoint and exit in seconds, so the task
picks up where it stopped.

### Where the output goes

On FASRC `$SCRATCH` is the bare mount point `/n/netscratch`, which is **not**
writable — the writable part is the lab subdirectory underneath it.  The script
resolves this itself, trying in order:

```
/n/netscratch/<lab>/Everyone/$USER/mhd
/n/netscratch/<lab>/Lab/$USER/mhd
/n/netscratch/<lab>/$USER/mhd
$SLURM_SUBMIT_DIR/runs                (last resort, next to the code)
```

with `<lab>` = `protopapas_lab` by default.  It prints `output root: ...` and
the free space at the top of every `.out` file.  If none of those work it exits
immediately with instructions rather than half-starting.  To choose explicitly:

```bash
ls -ld /n/netscratch/*/Everyone/$USER /n/netscratch/*/Lab/$USER 2>/dev/null

sbatch --export=ALL,N=128,NTDYN=12,TSNAP0=7.5,NSNAP=19,\
MHD_ROOT=/n/netscratch/protopapas_lab/Everyone/$USER/mhd cluster/run_suite.slurm
```

`analyse.slurm` resolves the root the same way; pass it the same `MHD_ROOT`.

### Storage per rung (21 runs, float32 snapshots, 7 fields)

| rung | per snapshot | per rung |
|---|---|---|
| 128^3, 19 snapshots | 59 MB | **23 GB** |
| 256^3, 11 snapshots | 470 MB | **108 GB** |
| 512^3, 11 snapshots | 3.8 GB | **870 GB** |

The 512^3 rung must go on netscratch, not a lab home directory.  If space is
tight, analyse each rung and delete its snapshots before starting the next —
the `.npz` files are ~2 MB per run and are all the plotting needs.  Lower
`NSNAP` to shrink it proportionally (11 snapshots over the analysis window is
already modest; the paper averages over more).

Memory per node with the `-DUSE_FLOAT` build: 0.4 GB at 128^3, 2.5 GB at
256^3, 17 GB at 512^3.

Runs checkpoint every 500 steps and resume by themselves, so hitting the wall
limit is harmless — resubmit the identical command: finished runs are skipped
and partial ones continue from `checkpoint.bin`.  That is what makes a 12 h
limit workable for 512^3 if you have to stay on a short partition.

`run_suite.slurm` needs no python, no module load and no conda env: `src/mhd`
is a pure C/OpenMP binary.

## 4. Analyse, then bring back ~40 MB instead of ~250 GB

This one *does* need python (numpy, scipy, matplotlib).  Submit it from a shell
with `venv_gpu` active, so the inherited `CONDA_PREFIX` fallback in the script
can rescue it if `source activate` silently does nothing:

```bash
sbatch --export=ALL,N=512 cluster/analyse.slurm
```

It prints which interpreter and library versions it ended up with, and fails
immediately with a clear message if numpy/scipy/matplotlib are missing, rather
than dying halfway through.

That writes `analysis_ax2.npz` / `analysis_ax0.npz` next to each run.  Those,
plus `history.txt` and `params.txt`, are everything the plotting needs:

```bash
rsync -avz --progress \
  --include '*/' --include '*.npz' --include 'history.txt' --include 'params.txt' \
  --exclude '*' \
  USER@CLUSTER:'$SCRATCH/mhd/n512/' /home/lgomez/MHD/runs/n512/
```

## 5. Figures and the comparison with the paper

```bash
cd /home/lgomez/MHD
python3 python/figures.py --root runs/n512 --outdir figures/n512
python3 python/compare_res.py runs/n64 runs/n128 runs/n256 runs/n512
```

`figures/n512/table1.txt` prints our Table 1 next to the published one, and
`compare_res.py` shows every coefficient as a function of resolution.
