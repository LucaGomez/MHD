# One Enzo 512^3 box at the Planck point -> 16 maps of 128x128 + measured cost

Physics: driven isothermal MHD turbulence matched to Stalpes, Collins &
Huffenberger (2024) at their Planck point, M_S = 4.7, M_A = 1.5.  Enzo with
Dedner/PLM/HLLD, gamma = 1.001, stochastic forcing with zeta = 1/2 at L/2 and
correlation time t_dyn, 10 t_dyn.  Calibration measured on the laptop puts the
achieved point at M_S 4.7, M_A 1.5 (the job prints the achieved values).

Maps: the final snapshot (10 t_dyn) is integrated through the FULL box along z
(perpendicular to the mean field) and the 512x512 sky is cut into 4x4 tiles ->
16 maps of 128x128 each of T, Q, U, E, B.

## 0. Copy the project (laptop)

```bash
rsync -avz --progress \
  --exclude 'runs/' --exclude 'runs_enzo/' --exclude 'runs_xcheck/' --exclude 'ic_n128/' \
  --exclude '__pycache__/' --exclude 'src/mhd' --exclude 'src/mhd_*' \
  /home/lgomez/MHD/ lgomez@CLUSTER:~/MHD/
```

## 1. Install Enzo (cluster login node, once, ~15-20 min)

```bash
cd ~/MHD
source activate venv_gpu          # its serial HDF5 is used for Enzo
mkdir -p slurm_logs
bash cluster/install_enzo.sh
```

It must end with `enzo.exe starts OK` and `python: numpy/scipy/h5py/matplotlib present`.

- `module load` fails -> `module spider openmpi`, then
  `ENZO_MODULES="gcc/<ver> openmpi/<ver>" bash cluster/install_enzo.sh`
- `enzo.exe does not start` with a GLIBCXX error -> build against the module
  HDF5 instead: `module load hdf5; HDF5_ROOT=$HDF5_HOME bash cluster/install_enzo.sh`
- python libs missing -> `pip install h5py scipy matplotlib` inside venv_gpu

## 2. Smoke test on the test partition (~10 min) -- do this first

The laptop can only run Enzo serially, so this is the first time Enzo runs on
MPI ranks.  It exercises everything the big job does, at 64^3.

```bash
ROOT=/n/netscratch/protopapas_lab/Everyone/lgomez/mhd
sbatch -p test --ntasks=8 -t 0-00:30 --mem=16000 \
  --export=ALL,N=64,NTDYN=1,DUMP_TDYN=0.5,TMIN_TDYN=0.5,MARGIN_H=0.1,MHD_ROOT=$ROOT \
  cluster/enzo_box.slurm
```

Check `slurm_logs/enzo512_<jobid>.out` for: `Enzo exit 0`, `axis check`,
`16 maps of 16x16`, `status : FINISHED`.

## 3. The 512^3 box

```bash
ROOT=/n/netscratch/protopapas_lab/Everyone/lgomez/mhd
sbatch --export=ALL,MHD_ROOT=$ROOT cluster/enzo_box.slurm
```

Defaults: sapphire, 1 node, 112 MPI ranks, 3 days, 200 GB.  When a segment
approaches the wall limit Enzo stops cleanly, writes a dump, and the job
resubmits itself (only if it made progress; at most 6 segments).  Restart is
exact: tested, the resumed run reproduces the uninterrupted one dump by dump.

## 4. Watching it

```bash
squeue -u $USER
R=$ROOT/enzo_box/ms4.7_ma1.5_n512_s100001

# live progress (StopTime is 1.0638)
grep -oE "time = [0-9.e+-]+ +cycle = [0-9]+" $R/enzo.log | tail -1

# cost so far + projected total (updated at the end of each segment)
python3 python/enzo_cost.py $R
```

## 5. Outputs, in $R

| file | contents |
|---|---|
| `maps/maps_4x4.npz` | `T Q U E B`, each (16, 128, 128) float32; `meta`, per-map `table` (JSON) |
| `maps/maps.png` | quick look at 8 of the maps |
| `snapshots/analysis_ax2.npz` | full-box statistics at 512^3, 6 snapshots over 5-10 t_dyn: slopes, BB/EE, r_TE |
| `snapshots/analysis_ax0.npz` | same, projected along the mean field |
| `cost_summary.txt` | wall time, core-hours, cycles, throughput, projection for N boxes |
| `DD*/` | Enzo dumps (~14 GB each, ~200 GB total). `DELETE_DUMPS=1` removes them after conversion |

```python
import numpy as np, json
d = np.load("maps/maps_4x4.npz")
Q, U = d["Q"], d["U"]                     # (16, 128, 128)
meta = json.loads(str(d["meta"]))         # M_S, M_A, c_s, time, axis, ...
table = json.loads(str(d["table"]))       # per map: <T>, median p, BB/EE
```

## What to expect

Measured on the laptop: Enzo does 3.6e5 zone-updates/s per core, and a 32^3
box needs ~6,300 cycles for 10 t_dyn.  Cycles grow at least in proportion to N,
so at 512^3:

| | estimate |
|---|---|
| cycles | ~100,000-130,000 |
| wall time on 112 ranks | ~4.5-6 days (2 segments) |
| core-hours | ~12,000-20,000 |
| disk | ~200 GB dumps + 22 GB snapshots |

These are extrapolations (per-core speed on sapphire, MPI efficiency at 112
ranks, and how much the time step shrinks at 512^3 are all unmeasured).  After
the first 3-day segment `python3 python/enzo_cost.py $R` prints the measured
projection, which is the number to plan the ensemble with.
