#!/bin/bash
# ---------------------------------------------------------------------------
#SBATCH --job-name=enzoscale
#SBATCH -p test
#SBATCH -N 1
#SBATCH --ntasks=112
#SBATCH -t 0-01:00
#SBATCH --mem=64000
#SBATCH -o slurm_logs/enzoscale_%j.out
#SBATCH -e slurm_logs/enzoscale_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gomezlucajavier@gmail.com
# ---------------------------------------------------------------------------
# How fast Enzo really is on a sapphire node, before spending days on 512^3.
#
#   sbatch --export=ALL,MHD_ROOT=/n/netscratch/protopapas_lab/Everyone/$USER/mhd \
#          cluster/enzo_scaling.slurm
#
# The same Planck-point box at 128^3 and 256^3, a fixed number of cycles
# (CYCLES) at several rank counts, no data dumps.  Per case it reports the
# cell-updates per core per second, the share of time in the MHD solve vs the
# ghost-zone exchange, and the CPU binding srun gave each rank -- ranks that
# all share one core show up as a single distinct mask.
#
# The 512^3 cost follows from the 256^3 rows: 8x the cells, ~2x the cycles.
# ---------------------------------------------------------------------------
set -u
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
mkdir -p slurm_logs

if [ -z "${ENZO:-}" ]; then
    # shellcheck disable=SC1091
    source cluster/enzo_env.sh || { echo "run bash cluster/install_enzo.sh first"; exit 1; }
fi
# shellcheck disable=SC1091
source cluster/mpi_launch.sh
ROOT=${MHD_ROOT:?set MHD_ROOT}
CASES=${CASES:-"128:1 128:8 128:28 128:56 128:112 256:28 256:56 256:112"}
CYCLES=${CYCLES:-30}
SKIP=${SKIP:-5}
AVAIL=${SLURM_NTASKS:-1}
OUT=$ROOT/enzo_scaling/$(date +%Y%m%d_%H%M%S)
mkdir -p "$OUT"
# one core per rank; drop CPU settings inherited from an interactive salloc
unset SLURM_CPUS_PER_TASK SLURM_TRES_PER_TASK

echo "== Enzo scaling  $(date)   node $(hostname)   $CYCLES cycles per case"
echo "   MPI launcher: $(enzo_launcher)"
echo "   results: $OUT/scaling.txt"
HDR="    N  ranks cycles cu/core/s  cu/s total  MHDsolve  ghosts  imbal  masks  efficiency"
echo "$HDR" | tee "$OUT/scaling.txt"

declare -A REF          # per-core throughput with the fewest ranks, per N
for c in $CASES; do
    N=${c%%:*}; R=${c##*:}
    if [ "$R" -gt "$AVAIL" ]; then echo "   skip $N^3 on $R ranks (only $AVAIL)"; continue; fi
    D=$OUT/n${N}_r${R}
    python3 python/enzo_setup.py --ms 4.7 --ma 1.5 --n "$N" --seed 100001 \
        --ntdyn 10 --dump-tdyn 1000 --out "$D" >/dev/null || exit 1
    echo "StopCycle = $CYCLES" >> "$D/turbulence.enzo"
    ( cd "$D" && ENZO_SRUN_EXTRA="--cpu-bind=verbose,cores" enzo_launch "$R" turbulence.enzo ) > "$D/enzo.log" 2>&1
    if [ ! -s "$D/performance.out" ]; then
        echo "   $N^3 on $R ranks: no performance.out -- see $D/enzo.log"; tail -5 "$D/enzo.log"; continue
    fi
    MASKS=$(grep -oE "mask 0x[0-9a-fA-F]+" "$D/enzo.log" | sort -u | wc -l)
    GOT=$(enzo_nproc "$D/performance.out")
    if [ "$R" -gt 1 ] && [ "$GOT" != "$R" ]; then
        echo "   $N^3 on $R ranks: Enzo reports $GOT MPI process(es) -- launch is broken, row skipped"
        rm -rf "$D"/DD[0-9]*; continue
    fi
    ROW=$(python3 python/enzo_perf.py "$D/performance.out" --skip "$SKIP" --oneline) || continue
    CU=$(echo "$ROW" | awk '{print $4}')
    [ -z "${REF[$N]:-}" ] && REF[$N]=$CU
    EFF=$(awk -v a="$CU" -v b="${REF[$N]}" 'BEGIN{printf "%.2f", a/b}')
    printf "%s %6s  %s\n" "$ROW" "$MASKS" "$EFF" | tee -a "$OUT/scaling.txt"
    rm -rf "$D"/DD[0-9]*
done

cat <<'TXT' | tee -a "$OUT/scaling.txt"

cu/core/s   cell updates per core per second (the cost driver); laptop reference
            for this Enzo build on one core: ~3.6e5 at 32^3-64^3
MHDsolve    share of cycle time in the MHD update;  ghosts: in ghost-zone exchange
imbal       slowest / fastest rank;  masks: distinct CPU masks srun bound ranks to
            (must equal the rank count -- 1 with many ranks = all on one core)
efficiency  cu/core/s relative to the fewest ranks run at the same N
512^3 cost  ~ 512^3 x cycles(512) / (cu/core/s x ranks);  cycles(512) ~ 2 x cycles(256)
TXT
echo "== done  $(date)"
