#!/bin/bash
# Measure OpenMP thread scaling on the node you will actually run on, so -c is
# chosen from data instead of guesswork.  The scheme is memory-bandwidth bound,
# so scaling flattens well before the core count -- asking for -c 112 when the
# curve flattens at 32 just burns fairshare.
#
# Run it inside an interactive allocation:
#     salloc -p test -c 112 -t 0-01:00 --mem=64000
#     bash cluster/bench.sh            # defaults to N=256
#     bash cluster/bench.sh 512        # if you have the memory for it
#
# Takes a couple of minutes.  Pick the smallest -c whose "speedup" is still
# close to the best one; that is the efficient point.
set -e
export LC_ALL=C          # keep printf/bc decimal points as "."
cd "$(dirname "$0")/.."
N=${1:-256}
STEPS_TDYN=${2:-0.02}     # starting slice; auto-enlarged if runs are too short

[ -x ./src/mhd ] || { echo "build first: bash cluster/setup.sh"; exit 1; }
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

MAXC=${SLURM_CPUS_PER_TASK:-$(nproc)}
echo "node has $MAXC cores available;  benchmarking N=$N"

# Calibrate the slice so that even the fastest configuration runs long enough
# for the timing to mean something (start-up and the first snapshot dominate
# anything under ~10 s).  Probe once at the largest thread count.
probe=$MAXC
rm -rf "$TMP/p"
ps=$(date +%s.%N)
pout=$(OMP_NUM_THREADS=$probe ./src/mhd n=$N ms=4 ma=1 ntdyn=$STEPS_TDYN \
       nsnap=0 out="$TMP/p" 2>&1 | tail -1)
pe=$(date +%s.%N)
pw=$(echo "$pe - $ps" | bc)
if [ "$(echo "$pw < 15" | bc)" = "1" ] && [ "$(echo "$pw > 0.2" | bc)" = "1" ]; then
    STEPS_TDYN=$(echo "scale=4; $STEPS_TDYN * 20 / $pw" | bc)
    echo "calibrating: probe took ${pw}s, raising the slice to ntdyn=$STEPS_TDYN"
fi
echo
printf "%8s %12s %12s %10s\n" threads "wall (s)" "steps/s" speedup
base=""
for t in 4 8 16 32 64 112; do
    [ "$t" -gt "$MAXC" ] && continue
    rm -rf "$TMP/b"
    st=$(date +%s.%N)
    out=$(OMP_NUM_THREADS=$t OMP_PROC_BIND=close OMP_PLACES=cores \
          ./src/mhd n=$N ms=4 ma=1 ntdyn=$STEPS_TDYN nsnap=0 out="$TMP/b" 2>&1 | tail -1)
    en=$(date +%s.%N)
    nst=$(echo "$out" | grep -o '[0-9]\+ steps' | grep -o '[0-9]\+' || echo 0)
    wall=$(echo "scale=1; ($en - $st)/1" | bc)
    rate=$(echo "scale=3; $nst / ($en - $st)" | bc)
    [ -z "$base" ] && base=$rate
    short=$(echo "$wall < 5" | bc); [ "$short" = "1" ] && tooshort=1
    sp=$(echo "scale=2; $rate / $base" | bc)
    printf "%8d %12s %12s %9sx\n" "$t" "$wall" "$rate" "$sp"
done
echo
if [ "${tooshort:-0}" = "1" ]; then
    echo "NOTE: some timings were still under 5 s -- treat those rows as noisy."
    echo "      Re-run with a bigger slice: bash cluster/bench.sh $N <ntdyn>"
    echo
fi
echo "Two ways to read this:"
echo "  * fastest wall clock  -> the thread count with the highest steps/s"
echo "  * best use of fairshare -> the largest thread count whose speedup is"
echo "    still close to proportional (speedup / (threads/4) near 1)"
echo "All 21 array tasks run at once, so per-run time IS the rung wall time."
echo "Memory needed per task (-DUSE_FLOAT):  128^3 0.4 GB, 256^3 2.5 GB,"
echo "512^3 17 GB -- so several tasks fit on one 990 GB sapphire node."
