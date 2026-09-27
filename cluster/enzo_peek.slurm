#!/bin/bash
# ---------------------------------------------------------------------------
#SBATCH --job-name=enzopeek
#SBATCH -p sapphire
#SBATCH -c 32
#SBATCH -t 0-03:00
#SBATCH --mem=128000
#SBATCH -o slurm_logs/enzopeek_%j.out
#SBATCH -e slurm_logs/enzopeek_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gomezlucajavier@gmail.com
# ---------------------------------------------------------------------------
# Early look at an Enzo box that is still in progress, from the dumps already
# on disk.  Does not touch the run: everything goes to <run>/preview/.
#
#   sbatch --export=ALL,MHD_ROOT=/n/netscratch/protopapas_lab/Everyone/$USER/mhd \
#          cluster/enzo_peek.slurm
#
#   1. dump table: simulation time, cycle, cycles per t_dyn
#   2. time step history from enzo.log -> where the cost comes from
#   3. dust statistics (M_S, M_A, slopes, BB/EE, r_TE) per dump with
#      t >= TMIN_TDYN, then the dust maps of the newest one
#
# Dumps before ~2 t_dyn are not yet in the statistically steady state: the
# numbers are a preview, not the result (the paper averages 5-10 t_dyn).
# Run it while the next Enzo segment is PENDING: a dump being written at the
# same moment would make the conversion fail (harmless, just resubmit).
# ---------------------------------------------------------------------------
set -u
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
mkdir -p slurm_logs

N=${N:-512}; MS=${MS:-4.7}; MA=${MA:-1.5}; SEED=${SEED:-100001}
TMIN_TDYN=${TMIN_TDYN:-1}
if [ -z "${RUN:-}" ]; then
    ROOT=${MHD_ROOT:?set MHD_ROOT or RUN}
    RUN=$ROOT/enzo_box/ms${MS}_ma${MA}_n${N}_s${SEED}
fi
[ -d "$RUN" ] || { echo "no run directory $RUN"; exit 1; }
TDYN=$(python3 -c "print(0.5/$MS)")
PV=$RUN/preview

INHERITED_CONDA_PREFIX="${CONDA_PREFIX:-}"
# Lmod/conda hooks reference unset variables; under set -u sourcing them kills
# the job silently.  Switch the checks off around them, then restore.
_had_e=0; _had_u=0; _had_pf=0
case $- in *e*) _had_e=1 ;; esac
case $- in *u*) _had_u=1 ;; esac
if [ -o pipefail ]; then _had_pf=1; fi
set +eu
if command -v module >/dev/null 2>&1; then
    module load python/3.10.9-fasrc01
    source activate venv_gpu || echo "aviso: 'source activate venv_gpu' devolvio error"
fi
if [ "$_had_e" = 1 ]; then set -e; fi
if [ "$_had_u" = 1 ]; then set -u; fi
if [ "$_had_pf" = 1 ]; then set -o pipefail; fi
if ! python3 -c "import numpy, scipy, h5py" 2>/dev/null; then
    if [ -n "$INHERITED_CONDA_PREFIX" ] && [ -x "$INHERITED_CONDA_PREFIX/bin/python3" ]; then
        export PATH="$INHERITED_CONDA_PREFIX/bin:$PATH"
    fi
fi
python3 -c "import numpy, scipy, h5py, matplotlib" || { echo "python env lacks numpy/scipy/h5py/matplotlib"; exit 1; }
export MHD_FFT_WORKERS=${SLURM_CPUS_PER_TASK:-8}

if command -v squeue >/dev/null && squeue -u "$USER" -h -t R -n enzo512 | grep -q .; then
    echo "note: an enzo512 segment is RUNNING; a dump written right now could break the conversion"
fi

rm -rf "$PV"; mkdir -p "$PV"
{
echo "=== preview of $RUN   ($(date))"
echo "    t_dyn = $TDYN code units, N = $N"

echo; echo "################ 1. dumps"
python3 - "$RUN" "$TDYN" <<'PY'
import sys, glob, os
run, tdyn = sys.argv[1], float(sys.argv[2])
rows = []
for d in sorted(glob.glob(os.path.join(run, "DD[0-9]*"))):
    p = os.path.join(d, "data" + d.rsplit("DD", 1)[1])
    kv = {}
    if os.path.exists(p):
        for l in open(p, errors="replace"):
            s = l.split("=", 1)
            if len(s) == 2:
                kv[s[0].strip()] = s[1].strip()
    if "InitialTime" not in kv:
        print(f"  {os.path.basename(d)}: incomplete, skipped"); continue
    rows.append((os.path.basename(d), float(kv["InitialTime"]) / tdyn, int(kv["InitialCycleNumber"])))
print(f"  {'dump':<8}{'t/t_dyn':>9}{'cycle':>9}{'cycles per t_dyn since previous':>34}")
prev = None
for name, t, c in rows:
    rate = ""
    if prev and t - prev[0] > 0.05:
        rate = f"{(c - prev[1]) / (t - prev[0]):,.0f}"
    print(f"  {name:<8}{t:9.3f}{c:9d}{rate:>34}")
    if not prev or t - prev[0] > 0.05:
        prev = (t, c)
print("  (dumps a few cycles apart at the end are restart points written when a segment stops)")
PY

echo; echo "################ 2. time step history (enzo.log)"
python3 - "$RUN/enzo.log" "$TDYN" "$N" <<'PY'
import sys, re
import numpy as np
log, tdyn, n = sys.argv[1], float(sys.argv[2]), int(sys.argv[3])
pat = re.compile(r"TopGrid dt = ([0-9.eE+-]+)\s+time = ([0-9.eE+-]+)\s+cycle = (\d+)")
dt, t = [], []
try:
    with open(log, errors="replace") as fh:
        for l in fh:
            m = pat.search(l)
            if m:
                dt.append(float(m.group(1))); t.append(float(m.group(2)))
except FileNotFoundError:
    sys.exit(f"  no {log}")
if not dt:
    sys.exit("  no 'TopGrid dt' lines found")
dt, t = np.array(dt), np.array(t) / tdyn
cfl, dx = 0.4, 1.0 / n
# the step is CFL * dx / (fastest signal speed): turn it back into that speed
print(f"  {'t/t_dyn':>12}{'cycles':>9}{'median dt':>12}{'min dt':>11}{'max signal speed':>18}")
edges = np.arange(0, t.max() + 0.25, 0.25)
for lo, hi in zip(edges[:-1], edges[1:]):
    m = (t >= lo) & (t < hi)
    if m.sum():
        print(f"  {lo:5.2f}-{hi:5.2f}{m.sum():9d}{np.median(dt[m]):12.3e}{dt[m].min():11.3e}"
              f"{cfl * dx / np.median(dt[m]):18.1f}")
print(f"\n  for scale: v_rms ~ M_S c_s ~ 5, fast speed of the mean field ~ 4.")
print("  A signal speed far above that comes from a few very-low-density cells where")
print("  the Alfven speed B/sqrt(rho) is huge: they set the step for the whole box.")
PY

echo; echo "################ 3. dust statistics per dump (t >= $TMIN_TDYN t_dyn)"
python3 python/enzo2snap.py "$RUN" "$PV/all" --tdyn "$TDYN" --tmin-tdyn "$TMIN_TDYN" --ms "$MS" --ma "$MA" \
    || { echo "!!! conversion failed"; exit 1; }
for S in "$PV"/all/snap_*.bin; do
    b=$(basename "$S" .bin)
    mkdir -p "$PV/$b"
    ln -sf "$S" "$PV/$b/snap_0000.bin"; ln -sf "${S%.bin}.txt" "$PV/$b/snap_0000.txt"
    cp "$PV/all/params.txt" "$PV/$b/"
    t=$(awk '$1=="time"{print $2}' "${S%.bin}.txt")
    echo "--- $b  (t = $(python3 -c "print(f'{$t/$TDYN:.2f}')") t_dyn)"
    python3 python/process_run.py "$PV/$b" --axis 2 | sed "s|$PV/||"
done
if [ "$(ls "$PV"/all/snap_*.bin | wc -l)" -gt 1 ]; then
    echo "--- average of all of them"
    python3 python/process_run.py "$PV/all" --axis 2 | sed "s|$PV/||"
fi
echo "    (paper at M_S>4: BB/EE 0.55, r_TE ~0.3; Planck: 0.53, 0.355;"
echo "     our solver at 512^3: BB/EE 0.90, r_TE 0.47)"

echo; echo "--- dust maps of every converted dump (16 tiles each; maps.png = newest)"
LAST=$(ls "$PV"/all/snap_*.bin | tail -1)
for S in "$PV"/all/snap_*.bin; do
    b=$(basename "$S" .bin)
    png=""; [ "$S" = "$LAST" ] && png="--png $PV/maps.png"
    echo "  $b:"
    python3 python/extract_maps.py "$S" --out "$PV/maps_${b}.npz" --tiles 4 --axis 2 $png | tail -2
done
echo; echo "=== done $(date)"
} 2>&1 | grep -v -i "runtimewarning" | tee "$PV/peek.txt"
echo "report: $PV/peek.txt    maps: $PV/maps.png"
