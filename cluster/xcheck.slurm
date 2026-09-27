#!/bin/bash
# ---------------------------------------------------------------------------
#SBATCH --job-name=xcheck
#SBATCH -p sapphire
#SBATCH -t 2-00:00
#SBATCH --mem=64000
#SBATCH -o slurm_logs/xcheck_%j.out
#SBATCH -e slurm_logs/xcheck_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gomezlucajavier@gmail.com
# ---------------------------------------------------------------------------
# Enzo <-> our solver crosscheck at a Planck-consistent point.  The same
# physical setup in both codes, and the SAME analysis (python/process_run.py)
# on both outputs.
#
#   Enzo is MPI-only, so it takes --ntasks; our solver is OpenMP, so -c:
#
#   sbatch --ntasks=64 --export=ALL,CODE=enzo,N=256,MHD_ROOT=$ROOT cluster/xcheck.slurm
#   sbatch -c 64      --export=ALL,CODE=mhd,N=256,MHD_ROOT=$ROOT  cluster/xcheck.slurm
#
#   then:  python3 python/compare_codes.py $ROOT/xcheck/mhd/<tag> $ROOT/xcheck/enzo/<tag>
#
# Defaults: M_S 4.7, M_A 1.5 (Stalpes et al.'s Planck point), 10 t_dyn, snapshots
# every 0.5 t_dyn over t in [5,10] t_dyn, calibration measured on the laptop.
# Resubmitting after a wall-limit kill resumes: Enzo from its last dump, our
# solver from checkpoint.bin.
# ---------------------------------------------------------------------------
set -eu
mkdir -p slurm_logs

CODE=${CODE:?set CODE=enzo or CODE=mhd}
N=${N:-256}
MS=${MS:-4.7}
MA=${MA:-1.5}
SEED=${SEED:-100001}
NTDYN=${NTDYN:-10}
CALIB=${CALIB:-1.132}        # Enzo: achieved/input M_S, measured
CS=${CS:-1.12}               # Enzo: mean c_s over the analysis window, measured
EDOTFAC=${EDOTFAC:-0.597}    # our solver: calibrated at the same point
ROOT=${MHD_ROOT:?set MHD_ROOT=/n/netscratch/<lab>/Everyone/$USER/mhd}
TAG=ms${MS}_ma${MA}_n${N}_s${SEED}
OUT=$ROOT/xcheck/$CODE/$TAG
TDYN=$(python3 -c "print(0.5/$MS)")

# --- python env for conversion / analysis (your usual block) ---------------
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
if ! python3 -c "import numpy, h5py" 2>/dev/null; then
    if [ -n "$INHERITED_CONDA_PREFIX" ] && [ -x "$INHERITED_CONDA_PREFIX/bin/python3" ]; then
        export PATH="$INHERITED_CONDA_PREFIX/bin:$PATH"
    fi
fi
python3 -c "import numpy, scipy, h5py, matplotlib; print('python ok')"

echo "== $CODE  M_S=$MS M_A=$MA N=$N seed=$SEED  -> $OUT   ($(date))"

if [ "$CODE" = enzo ]; then
    source cluster/enzo_env.sh || { echo "run bash cluster/install_enzo.sh first"; exit 1; }
    RUN=$ROOT/enzo_runs/$TAG
    if [ ! -f "$RUN/turbulence.enzo" ]; then
        python3 python/enzo_setup.py --ms $MS --ma $MA --n $N --seed $SEED \
            --ntdyn $NTDYN --calib $CALIB --cs $CS --out $RUN
    fi
    cd "$RUN"
    NT=${SLURM_NTASKS:-1}
    if [ -f RunFinished ]; then
        echo "Enzo already finished"
    else
        last=$(ls -d DD[0-9]* 2>/dev/null | sort | tail -1)
        if [ -n "$last" ]; then
            echo "restarting Enzo from $last"
            srun -n "$NT" "$ENZO" -r "$last/data${last#DD}" >> enzo.log 2>&1
        else
            srun -n "$NT" "$ENZO" turbulence.enzo > enzo.log 2>&1
        fi
    fi
    cd - >/dev/null
    python3 python/enzo2snap.py "$RUN" "$OUT" --tdyn "$TDYN" --tmin-tdyn 5 --ms $MS --ma $MA
else
    [ -x ./src/mhd ] || { echo "build first: bash cluster/setup.sh"; exit 1; }
    export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16} OMP_PROC_BIND=close OMP_PLACES=cores
    ./src/mhd n=$N ms=$MS ma=$MA ntdyn=$NTDYN tsnap0=5 nsnap=11 \
              edotfac=$EDOTFAC seed=$SEED chkint=500 out="$OUT"
fi

export MHD_FFT_WORKERS=${SLURM_CPUS_PER_TASK:-${SLURM_NTASKS:-8}}
python3 python/process_run.py "$OUT" --axis 2
python3 python/process_run.py "$OUT" --axis 0
echo "== done  ($(date))"
