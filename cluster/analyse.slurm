#!/bin/bash
# ---------------------------------------------------------------------------
#SBATCH --job-name=mhdana
#SBATCH -p sapphire
#SBATCH -c 16
#SBATCH -t 0-08:00
#SBATCH --mem=48000
#SBATCH -o slurm_logs/mhdana_%A_%a.out
#SBATCH -e slurm_logs/mhdana_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gomezlucajavier@gmail.com
# ---------------------------------------------------------------------------
# Spectra, projections and E/B decomposition.  THIS is the job that needs
# python (numpy, scipy, matplotlib).
#
# NO --array by default: one job loops over all 21 runs.  A 21-task array
# would trip QOSMaxSubmitJobPerUserLimit once the simulation arrays are also
# queued.  Split it only if you need to:
#     sbatch --array=0-2 --export=ALL,NCHUNKS=3,N=512,MHD_ROOT=... cluster/analyse.slurm
#
#   sbatch --export=ALL,N=512,MHD_ROOT=/n/netscratch/<lab>/Everyone/$USER/mhd \
#          cluster/analyse.slurm
#
# Cost: the FFTs are threaded (scipy.fft, workers = -c) and single precision.
# ~1 min per run at 128^3, ~10 min at 256^3, ~1 h at 512^3, so the default
# 8 h covers a whole 512^3 suite in one job.  Runs already analysed are
# skipped, so resubmitting after a timeout is cheap.
# ---------------------------------------------------------------------------
set -u
mkdir -p slurm_logs
failed=""

# CONDA_PREFIX del shell que hizo el submit.  sbatch exporta el entorno por
# defecto, asi que si lanzaste con venv_gpu activo esto ya apunta al env bueno.
# Es el fallback fiable: en un shell no interactivo no se sourcea ~/.bashrc y
# `source activate` puede no hacer nada, dejando el python3 del modulo (sin
# numpy/scipy).
INHERITED_CONDA_PREFIX="${CONDA_PREFIX:-}"

# Lmod/conda hooks reference unset variables; under set -u/-e sourcing them
# kills the job silently.  Switch the checks off around them, then restore.
# Read the flags from $- directly: $(set +o) runs in a subshell, where bash
# drops errexit, so restoring from it would silently leave -e off.
_had_e=0; _had_u=0; _had_pf=0
case $- in *e*) _had_e=1 ;; esac
case $- in *u*) _had_u=1 ;; esac
if [ -o pipefail ]; then _had_pf=1; fi
set +eu
module load python/3.10.9-fasrc01
source activate venv_gpu || echo "aviso: 'source activate venv_gpu' devolvio error"
if [ "$_had_e" = 1 ]; then set -e; fi
if [ "$_had_u" = 1 ]; then set -u; fi
if [ "$_had_pf" = 1 ]; then set -o pipefail; fi

if ! python3 -c "import numpy" 2>/dev/null; then
    if [ -n "$INHERITED_CONDA_PREFIX" ] && [ -x "$INHERITED_CONDA_PREFIX/bin/python3" ]; then
        echo "usando python heredado: $INHERITED_CONDA_PREFIX/bin/python3"
        export PATH="$INHERITED_CONDA_PREFIX/bin:$PATH"
    fi
fi

# fallar temprano y con un mensaje claro, en vez de a mitad del analisis
python3 - <<'PY'
import sys, numpy, scipy, matplotlib
print("python :", sys.executable)
print("numpy  :", numpy.__version__, " scipy:", scipy.__version__,
      " mpl:", matplotlib.__version__)
PY

export MHD_FFT_WORKERS=${SLURM_CPUS_PER_TASK:-8}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

N=${N:-512}
if [ -n "${MHD_ROOT:-}" ]; then
    ROOT="$MHD_ROOT/n$N"
elif [ -z "${ROOT:-}" ]; then
    LAB=${MHD_LAB:-protopapas_lab}
    for cand in "/n/netscratch/$LAB/Everyone/$USER/mhd" \
                "/n/netscratch/$LAB/Lab/$USER/mhd" \
                "/n/netscratch/$LAB/$USER/mhd" \
                "${SLURM_SUBMIT_DIR:-$PWD}/runs"; do
        [ -d "$cand/n$N" ] && { ROOT="$cand/n$N"; break; }
    done
fi
echo "reading from: ${ROOT:-<unresolved>}"
DIRS=($(ls -d $ROOT/ms*_ma* 2>/dev/null))
if [ ${#DIRS[@]} -eq 0 ]; then
    echo "no runs found in '${ROOT:-<unresolved>}' -- pass MHD_ROOT=... if the"
    echo "output went elsewhere (run_suite.slurm prints 'output root:')"
    exit 1
fi

k=${SLURM_ARRAY_TASK_ID:-0}
NCHUNKS=${NCHUNKS:-1}
echo "${#DIRS[@]} runs found; this task handles every ${NCHUNKS}th starting at $k"

for ((j=k; j<${#DIRS[@]}; j+=NCHUNKS)); do
    D=${DIRS[$j]}
    newest=$(ls -t "$D"/snap_*.bin 2>/dev/null | head -1)
    if [ -z "$newest" ]; then
        echo "--- $(basename $D): no snapshots, skipping"; continue
    fi
    if [ -f "$D/analysis_ax2.npz" ] && [ -f "$D/analysis_ax0.npz" ] \
       && [ "$D/analysis_ax0.npz" -nt "$newest" ]; then
        echo "--- $(basename $D): already analysed, skipping"; continue
    fi
    echo "--- $(basename $D)  ($(date))"
    # Do NOT let one bad run kill the whole job: `set -e` plus 21 runs in a
    # loop means a single failure loses the other 20.  Record and continue.
    ok=1
    python3 python/process_run.py "$D" --axis 2 || ok=0
    python3 python/process_run.py "$D" --axis 0 || ok=0
    if [ $ok = 0 ]; then
        echo "!!! $(basename $D) FAILED, continuing with the rest"
        failed="$failed $(basename $D)"
    fi
done
if [ -n "${failed:-}" ]; then
    echo "analysis task $k finished WITH FAILURES:$failed  ($(date))"
    exit 1
fi
echo "analysis task $k done, all runs OK  ($(date))"
