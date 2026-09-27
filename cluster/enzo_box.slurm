#!/bin/bash
# ---------------------------------------------------------------------------
#SBATCH --job-name=enzo512
#SBATCH -p sapphire
#SBATCH -N 1
#SBATCH --ntasks=112
#SBATCH -t 3-00:00
#SBATCH --mem=200000
#SBATCH -o slurm_logs/enzo512_%j.out
#SBATCH -e slurm_logs/enzo512_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gomezlucajavier@gmail.com
# ---------------------------------------------------------------------------
# ONE Enzo 512^3 box of driven isothermal MHD turbulence at the Planck point
# (M_S 4.7, M_A 1.5; Stalpes, Collins & Huffenberger 2024), then 16 dust maps
# of 128x128 and the measured cost.
#
#   sbatch --export=ALL,MHD_ROOT=/n/netscratch/protopapas_lab/Everyone/$USER/mhd \
#          cluster/enzo_box.slurm
#
# What the job does, unattended:
#   1. writes the paper-matched Enzo parameter file (python/enzo_setup.py)
#   2. runs Enzo on MPI ranks, or restarts it from the newest COMPLETE dump
#   3. stops Enzo cleanly before the wall limit (StopCPUTime = time left minus
#      MARGIN_H), so a limit kill never leaves a half-written dump
#   4. if the run isn't done, resubmits itself (only if this segment made
#      progress, at most MAX_SEGMENTS times)
#   5. when done: converts the dumps from 5 t_dyn on, runs the same analysis
#      as for our solver, cuts 16 x 128^2 maps from the last snapshot, and
#      prints the cost (python/enzo_cost.py)
#
# Output under $MHD_ROOT/enzo_box/<tag>/:
#   DD*/ (Enzo dumps), enzo.log, cost.log, cost_summary.txt,
#   snapshots/ (analysis_ax2.npz, analysis_ax0.npz), maps/maps_16x128.npz + maps.png
#
# Smoke test on the test partition first (64^3, 1 t_dyn, a few minutes):
#   sbatch -p test --ntasks=8 -t 0-00:30 --mem=16000 \
#          --export=ALL,N=64,NTDYN=1,TMIN_TDYN=0.5,MHD_ROOT=... cluster/enzo_box.slurm
# ---------------------------------------------------------------------------
set -u
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
mkdir -p slurm_logs

N=${N:-512}
MS=${MS:-4.7}
MA=${MA:-1.5}
SEED=${SEED:-100001}
NTDYN=${NTDYN:-10}           # paper: driven 10 t_dyn
DUMP_TDYN=${DUMP_TDYN:-1.0}  # one dump per t_dyn (13.9 GB each at 512^3)
TMIN_TDYN=${TMIN_TDYN:-5}    # paper: analyse the second half
CALIB=${CALIB:-1.132}        # measured: achieved/input M_S for this template
CS=${CS:-1.12}               # measured: mean c_s over the analysis window
MARGIN_H=${MARGIN_H:-1.5}    # stop Enzo this long before the wall limit
MAX_SEGMENTS=${MAX_SEGMENTS:-6}
TILES=${TILES:-4}
DELETE_DUMPS=${DELETE_DUMPS:-0}
VA_CAP=${VA_CAP:-}           # e.g. 50: Alfven-speed limiter (needs the patched Enzo)
ROOT=${MHD_ROOT:?set MHD_ROOT, e.g. /n/netscratch/protopapas_lab/Everyone/$USER/mhd}

TAG=ms${MS}_ma${MA}_n${N}_s${SEED}
RUN=$ROOT/enzo_box/$TAG
NT=${SLURM_NTASKS:-1}
TDYN=$(python3 -c "print(0.5/$MS)")
STOPTIME=$(python3 -c "print($NTDYN*0.5/$MS)")

if [ -z "${ENZO:-}" ]; then
    # shellcheck disable=SC1091
    source cluster/enzo_env.sh 2>/dev/null || { echo "no cluster/enzo_env.sh: run bash cluster/install_enzo.sh"; exit 1; }
fi
[ -x "$ENZO" ] || { echo "ENZO=$ENZO is not executable"; exit 1; }
# shellcheck disable=SC1091
source cluster/mpi_launch.sh

# Refuse to start if enzo.exe can't load its libraries here (e.g. an old system
# libstdc++ shadowing the compiler's): fail in seconds with the reason.
PRE=$( { "$ENZO" -h 2>&1; ldd "$ENZO" 2>&1; } | grep -iE "not found|error while loading|cannot open shared" | sort -u )
if [ -n "$PRE" ]; then
    echo "!! enzo.exe cannot load its libraries on this node:"
    echo "$PRE" | head -8 | sed 's/^/   /'
    echo "   LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-}"
    echo "   rerun:  bash cluster/install_enzo.sh --env-only"
    exit 1
fi

echo "== $TAG   ranks=$NT   $(date)"
echo "   run dir: $RUN"

# --- 1. parameter file ------------------------------------------------------
if [ ! -f "$RUN/turbulence.enzo" ]; then
    python3 python/enzo_setup.py --ms "$MS" --ma "$MA" --n "$N" --seed "$SEED" \
        --ntdyn "$NTDYN" --dump-tdyn "$DUMP_TDYN" --calib "$CALIB" --cs "$CS" \
        ${VA_CAP:+--va-cap "$VA_CAP"} --out "$RUN" || exit 1
fi
touch "$RUN/cost.log"

# newest dump whose parameter file, hierarchy and every data file exist
latest_dump() {
    local d num p f ok
    for d in $(ls -d "$RUN"/DD[0-9]* 2>/dev/null | sort -r); do
        num=${d##*/DD}; p="$d/data$num"
        [ -s "$p" ] && [ -s "$p.hierarchy" ] && grep -q "^InitialTime" "$p" || continue
        ok=1
        for f in $(awk -F'= ' '/^BaryonFileName/{print $2}' "$p.hierarchy" | sort -u); do
            [ -s "$RUN/${f#./}" ] || { ok=0; break; }
        done
        [ $ok = 1 ] && { echo "$p"; return 0; }
    done
    return 1
}
dump_field() { awk -v k="$2" '$1==k {print $3; exit}' "$1"; }
is_done() {   # finished = newest complete dump is at StopTime
    local p t
    p=$(latest_dump) || return 1
    t=$(dump_field "$p" InitialTime)
    awk -v t="$t" -v s="$STOPTIME" 'BEGIN{exit !(t >= s*(1-1e-4))}'
}

# --- 2/3. run or restart Enzo ----------------------------------------------
if is_done; then
    echo "   Enzo already reached StopTime = $STOPTIME"
else
    # seconds Slurm will still give us (falls back to 72 h outside Slurm)
    LEFT=$(python3 - "${SLURM_JOB_ID:-}" <<'PY'
import subprocess, sys
jid = sys.argv[1]
s = ""
if jid:
    try:
        # Python 3.6 compatible (system python3 on Rocky 8): no capture_output/text
        s = subprocess.run(["squeue", "-h", "-j", jid, "-o", "%L"],
                           stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                           universal_newlines=True).stdout.strip()
    except Exception:
        s = ""
def parse(x):
    if not x or x in ("UNLIMITED", "NOT_SET", "INVALID"):
        return 72 * 3600
    d = 0
    if "-" in x:
        d, x = x.split("-"); d = int(d)
    p = [int(v) for v in x.split(":")]
    while len(p) < 3:
        p.insert(0, 0)
    return d * 86400 + p[0] * 3600 + p[1] * 60 + p[2]
print(parse(s))
PY
)
    STOPCPU=$(python3 -c "print(max(300, int($LEFT - $MARGIN_H*3600)))")
    echo "   wall time left: $((LEFT/3600)) h;  Enzo StopCPUTime = $STOPCPU s"

    if p0=$(latest_dump); then
        c0=$(dump_field "$p0" InitialCycleNumber); t0=$(dump_field "$p0" InitialTime)
        # CPUTime is stored in the dump and would count the previous segments
        # against the new limit: reset it, and set this segment's limit.
        sed -i -E "s/^(CPUTime[[:space:]]*=).*/\1 0/; s/^(StopCPUTime[[:space:]]*=).*/\1 $STOPCPU/" "$p0"
        grep -q "^StopCPUTime" "$p0" || echo "StopCPUTime = $STOPCPU" >> "$p0"
        REL=${p0#$RUN/}
        echo "   restarting from $REL  (t=$t0, cycle $c0)"
        ARGS=(-r "$REL")
    else
        c0=0
        sed -i '/^StopCPUTime/d' "$RUN/turbulence.enzo"
        echo "StopCPUTime = $STOPCPU" >> "$RUN/turbulence.enzo"
        echo "   starting from t=0"
        ARGS=(turbulence.enzo)
    fi

    [ "$NT" -gt 1 ] && echo "   MPI launcher: $(enzo_launcher)"
    tA=$(date +%s)
    ( cd "$RUN" || exit 1; enzo_launch "$NT" "${ARGS[@]}" ) >> "$RUN/enzo.log" 2>&1
    rc=$?
    tB=$(date +%s)
    # Did Enzo really run as ONE job on $NT ranks?  If not, stop here: going
    # on (or resubmitting) would burn a node at one-core speed.
    if [ "$NT" -gt 1 ]; then
        GOT=$(enzo_nproc "$RUN/performance.out")
        if [ "$GOT" != "$NT" ]; then
            echo "!! Enzo reports $GOT MPI process(es), not $NT: the MPI launch is broken."
            echo "   Each rank ran as a separate 1-process Enzo.  Try MPI_LAUNCH=mpirun or"
            echo "   MPI_LAUNCH=srun-pmix.  NOT resubmitting."
            echo "$tA $tB $NT - - -" >> "$RUN/cost.log"
            exit 1
        fi
        echo "   Enzo confirms $GOT MPI processes"
    fi
    if p1=$(latest_dump); then
        c1=$(dump_field "$p1" InitialCycleNumber); t1=$(dump_field "$p1" InitialTime)
    else
        c1=$c0; t1=0
    fi
    echo "$tA $tB $NT $c0 $c1 $t1" >> "$RUN/cost.log"
    echo "   Enzo exit $rc after $(( (tB-tA)/60 )) min: cycles $c0 -> $c1, t -> $t1 of $STOPTIME"
    # Enzo's positivity fallbacks (a cell whose energy or density would go
    # negative in a strong shock is repaired locally).  A handful during the
    # start-up transient is normal: 6 in the whole first t_dyn at 64^3, the
    # same serial and on 8 MPI ranks.  A count growing steadily through the
    # run is not.
    ntau=$(grep -c "tau < 0" "$RUN/enzo.log" 2>/dev/null); nrho=$(grep -c "rho <0" "$RUN/enzo.log" 2>/dev/null)
    # "Fallback failed, give up" = a cell whose update failed even at first
    # order; Enzo keeps the state and moves on.  Seen once in the 128^3
    # start-up (time step then shrank 3.4x and recovered).  If this count keeps
    # growing after the first t_dyn, stop and investigate before trusting the run.
    nfo=$(grep -c "Falling back to zero order" "$RUN/enzo.log" 2>/dev/null); ngu=$(grep -c "Fallback failed" "$RUN/enzo.log" 2>/dev/null)
    echo "   numerical fallbacks in enzo.log so far: tau<0 ${ntau:-0}, rho<0 ${nrho:-0}, first-order ${nfo:-0}, failed ${ngu:-0}"
    if [ "${ngu:-0}" -gt 0 ]; then
        grep -B1 "Fallback failed" "$RUN/enzo.log" | grep -oE "cycle = [0-9]+" | tail -1 | sed 's/^/   last failed fallback before /'
    fi

    if ! is_done; then
        nseg=$(grep -c . "$RUN/cost.log")
        if [ "$c1" -le "$c0" ]; then
            echo "!! no progress in this segment -- NOT resubmitting.  Look at $RUN/enzo.log"
            tail -20 "$RUN/enzo.log"
            exit 1
        fi
        if [ "$nseg" -ge "$MAX_SEGMENTS" ]; then
            echo "!! $nseg segments already -- NOT resubmitting (raise MAX_SEGMENTS)"
            exit 1
        fi
        if [ -n "${SLURM_JOB_ID:-}" ] && command -v sbatch >/dev/null; then
            nj=$(sbatch --parsable --dependency=afterany:"$SLURM_JOB_ID" --export=ALL \
                        "${SLURM_SUBMIT_DIR}/cluster/enzo_box.slurm")
            echo "   not finished: resubmitted as job $nj (segment $((nseg+1)))"
        else
            echo "   not finished: rerun this script to continue"
        fi
        python3 python/enzo_cost.py "$RUN" 2>/dev/null | sed 's/^/   /'
        exit 0
    fi
fi

# --- 5. convert, analyse, cut maps, cost -----------------------------------
echo "== post-processing  $(date)"
INHERITED_CONDA_PREFIX="${CONDA_PREFIX:-}"
# Lmod and conda's activate/deactivate hooks reference unset variables: under
# 'set -u' (or -e) sourcing them kills this shell SILENTLY -- which is how the
# first smoke test died right after "== post-processing".  Switch the checks
# off around them and restore the caller's options afterwards.
# Read the flags from $- directly: $(set +o) runs in a subshell, where bash
# drops errexit, so restoring from it would silently leave -e off.
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
export MHD_FFT_WORKERS=$NT

# -o pipefail: the cost report is piped through tee, and without pipefail a
# failing python step there would still let the job report success.
set -e -o pipefail
python3 python/enzo2snap.py "$RUN" "$RUN/snapshots" --tdyn "$TDYN" --tmin-tdyn "$TMIN_TDYN" --ms "$MS" --ma "$MA"
python3 python/process_run.py "$RUN/snapshots" --axis 2
python3 python/process_run.py "$RUN/snapshots" --axis 0
LAST=$(ls "$RUN"/snapshots/snap_*.bin | tail -1)
mkdir -p "$RUN/maps"
python3 python/extract_maps.py "$LAST" --out "$RUN/maps/maps_${TILES}x${TILES}.npz" \
        --tiles "$TILES" --axis 2 --png "$RUN/maps/maps.png"
python3 python/enzo_cost.py "$RUN" | tee "$RUN/cost_summary.txt"
if [ "$DELETE_DUMPS" = 1 ]; then rm -rf "$RUN"/DD[0-9]*; echo "   Enzo dumps deleted"; fi
echo "== done  $(date)"
