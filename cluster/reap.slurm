#!/bin/bash
# ---------------------------------------------------------------------------
#SBATCH --job-name=mhdreap
#SBATCH -p sapphire
#SBATCH -c 16
#SBATCH -t 2-00:00
#SBATCH --mem=48000
#SBATCH -o slurm_logs/mhdreap_%j.out
#SBATCH -e slurm_logs/mhdreap_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gomezlucajavier@gmail.com
# ---------------------------------------------------------------------------
# Unattended housekeeping: every WATCH seconds, analyse any run that has
# reached its end time and delete its snapshots.  Keeps the suite's disk
# footprint bounded while the rest of it is still computing.
#
#   sbatch --export=ALL,N=512,MHD_ROOT=/n/netscratch/<lab>/Everyone/$USER/mhd \
#          cluster/reap.slurm
#
# Stops on its own when nothing is left running, or after MAXH hours.
# Keeps analysis_ax*.npz, history.txt, params.txt and checkpoint.bin.
# Set KEEP_LAST=1 for a resolution a later rung still upsamples from.
# ---------------------------------------------------------------------------
set -e
mkdir -p slurm_logs

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
        export PATH="$INHERITED_CONDA_PREFIX/bin:$PATH"
    fi
fi
python3 -c "import numpy, scipy, matplotlib, sys; print('python:', sys.executable)"

export MHD_FFT_WORKERS=${SLURM_CPUS_PER_TASK:-8}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

N=${N:-512}
WATCH=${WATCH:-1200}
MAXH=${MAXH:-44}
OPTS="--yes --watch $WATCH --max-hours $MAXH"
[ "${KEEP_LAST:-0}" = "1" ] && OPTS="$OPTS --keep-last"

bash cluster/reap.sh "$N" $OPTS
