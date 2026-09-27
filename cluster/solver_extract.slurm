#!/bin/bash
# ---------------------------------------------------------------------------
#SBATCH --job-name=solverx
#SBATCH -p sapphire
#SBATCH -c 32
#SBATCH -t 0-08:00
#SBATCH --mem=128000
#SBATCH -o slurm_logs/solverx_%j.out
#SBATCH -e slurm_logs/solverx_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gomezlucajavier@gmail.com
# ---------------------------------------------------------------------------
# Everything useful from our solver's ladder (128^3 -> 256^3 -> 512^3), in one
# report.  No new simulations: only analysis of what is already on disk.
#
#   sbatch --export=ALL,MHD_ROOT=/n/netscratch/protopapas_lab/Everyone/$USER/mhd \
#          cluster/solver_extract.slurm
#
#   1. health of every run at each resolution
#   2. analysis of any run not yet analysed (skips the ones that are)
#   3. convergence of Table 1, <BB/EE> and <r_TE> across resolutions
#   4. Table 1 and figures at the highest resolution (both projections)
#   5. at the highest resolution, on runs that still have snapshots:
#        polarization-angle order vs M_A, line-of-sight depth dependence
#        (M_S > 4 runs) and snapshot decorrelation
#
# Output: figures/solver_extract/ (report.txt + figures).  Small: rsync it back.
# ---------------------------------------------------------------------------
set -u
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
mkdir -p slurm_logs

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
if ! python3 -c "import numpy, scipy, matplotlib" 2>/dev/null; then
    if [ -n "$INHERITED_CONDA_PREFIX" ] && [ -x "$INHERITED_CONDA_PREFIX/bin/python3" ]; then
        export PATH="$INHERITED_CONDA_PREFIX/bin:$PATH"
    fi
fi
python3 -c "import numpy, scipy, matplotlib" || { echo "python env lacks numpy/scipy/matplotlib"; exit 1; }

ROOT=${MHD_ROOT:?set MHD_ROOT}
RES_REQ=${RES:-"128 256 512"}
SECTIONS=${SECTIONS:-"1 2 3 4 5"}     # e.g. SECTIONS=5 to redo only the snapshot checks
want() { case " $SECTIONS " in *" $1 "*) return 0 ;; esac; return 1; }
OUT=${OUT:-figures/solver_extract}
export MHD_FFT_WORKERS=${SLURM_CPUS_PER_TASK:-8}
export MHD_ROOT
mkdir -p "$OUT"
REPORT=$OUT/report${SECTIONS// /}.txt

RES=""; for N in $RES_REQ; do [ -d "$ROOT/n$N" ] && RES="$RES $N"; done
RES=${RES# }
[ -n "$RES" ] || { echo "no suites under $ROOT for N in: $RES_REQ"; exit 1; }

last_ms() { awk '$1=="ms"{print $2}' "$(ls "$1"/snap_*.txt 2>/dev/null | tail -1)" 2>/dev/null; }

{
echo "=== our solver, resolution ladder: $RES   ($(date))"
echo "    root: $ROOT"

if want 1; then
echo; echo "################ 1. health"
bash cluster/check_health.sh $RES
fi

if want 2; then
echo; echo "################ 2. analysis of runs not yet analysed"
for N in $RES; do
    for D in "$ROOT"/n$N/ms*_ma*; do
        [ -d "$D" ] || continue
        newest=$(ls -t "$D"/snap_*.bin 2>/dev/null | head -1)
        [ -n "$newest" ] || continue                           # reaped: nothing to redo
        # runs have finished, so an existing analysis is final (mtimes are not
        # trustworthy after rsync)
        if [ -f "$D/analysis_ax2.npz" ] && [ -f "$D/analysis_ax0.npz" ]; then continue; fi
        echo "--- analysing $(basename "$D") at $N^3"
        python3 python/process_run.py "$D" --axis 2 || echo "!!! $D axis 2 failed"
        python3 python/process_run.py "$D" --axis 0 || echo "!!! $D axis 0 failed"
    done
done

fi

# only resolutions that actually have analysed runs enter the comparison
ARES=""; for N in $RES; do
    ls "$ROOT"/n$N/*/analysis_ax2.npz >/dev/null 2>&1 && ARES="$ARES $N"
done
ARES=${ARES# }
[ -n "$ARES" ] || { echo "nothing analysed anywhere"; exit 1; }
NMAX=${ARES##* }
echo; echo "analysed resolutions: $ARES   (highest: $NMAX)"

if want 3; then
echo; echo "################ 3. convergence across resolutions"
python3 python/compare_res.py $(for N in $ARES; do echo "$ROOT/n$N"; done) --out "$OUT/convergence.txt"

fi

if want 4; then
echo; echo "################ 4. Table 1 at ${NMAX}^3"
python3 python/figures.py --root "$ROOT/n$NMAX" --outdir "$OUT/n$NMAX" --label "N=$NMAX" \
    || echo "!!! figures (perpendicular) failed"
python3 python/figures.py --root "$ROOT/n$NMAX" --outdir "$OUT/n${NMAX}_parallel" --axis 0 \
    --label "N=$NMAX along B0" >/dev/null || echo "!!! figures (parallel) failed"

fi

if want 5; then
echo; echo "################ 5. ${NMAX}^3 runs that still have snapshots"
HAS=(); PLANCK=()
for D in "$ROOT"/n$NMAX/ms*_ma*; do
    [ "$(ls "$D"/snap_*.bin 2>/dev/null | wc -l)" -ge 4 ] || continue
    ms=$(last_ms "$D")
    # "nan" must not pass as a string comparison: force numeric (nan+0 = 0)
    if ! awk -v m="${ms:-0}" 'BEGIN{exit !(m+0 > 0)}'; then
        echo "skipping $(basename "$D"): non-finite M_S (blown-up run)"; continue
    fi
    HAS+=("$D")
    awk -v m="$ms" 'BEGIN{exit !(m+0 > 4)}' && PLANCK+=("$D")
done
echo "with snapshots: ${#HAS[@]}   (M_S > 4: ${#PLANCK[@]})"
if [ ${#HAS[@]} -gt 0 ]; then
    echo; echo "--- 5a. polarization-angle order vs M_A"
    python3 python/pol_angle_stats.py "${HAS[@]}"
fi
if [ ${#PLANCK[@]} -gt 0 ]; then
    echo; echo "--- 5b. line-of-sight depth dependence (M_S > 4)"
    python3 python/los_depth.py "${PLANCK[@]}" --max-snaps 3
fi
if [ ${#HAS[@]} -gt 1 ]; then
    echo; echo "--- 5c. snapshot decorrelation"
    python3 python/snapshot_correlation.py "${HAS[@]}"
fi
fi
echo; echo "=== done $(date)"
} 2>&1 | grep -v -i "runtimewarning\|arrmean\|ret = ret" | tee "$REPORT"
echo "report: $REPORT"
