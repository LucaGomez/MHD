#!/bin/bash
# Analyse finished runs and then delete their snapshots, so the suite's disk
# footprint stays bounded while the rest of it is still computing.
#
# A full ladder is ~1.1 TB at peak (512^3 snapshots alone are 869 GB), which
# is right at a 1.1 TB quota.  Reaping finished runs as they land keeps the
# peak far below that.
#
#   bash cluster/reap.sh 512                 # DRY RUN: says what it would do
#   bash cluster/reap.sh 512 --yes           # actually delete
#   bash cluster/reap.sh 256 --yes --keep-last   # keep 1 snapshot per run
#   bash cluster/reap.sh 512 --yes --watch 1200 --max-hours 36
#
# What is kept, always: analysis_ax*.npz, history.txt, params.txt, snap_*.txt,
# and checkpoint.bin (so a run can still resume).  Only snap_*.bin go.
#
# --keep-last is REQUIRED for a resolution that a later rung still upsamples
# from: run_suite.slurm resolves its upsample source by globbing the coarser
# suite's snapshots, and exits with an error if none are left.
set -u
cd "$(dirname "$0")/.."

N=${1:?usage: reap.sh <N> [--yes] [--keep-last] [--watch SECONDS] [--max-hours H]}
shift || true
DO_IT=0; KEEP_LAST=0; WATCH=0; MAXH=0
while [ $# -gt 0 ]; do
    case "$1" in
        --yes) DO_IT=1 ;;
        --keep-last) KEEP_LAST=1 ;;
        --watch) WATCH=${2:?}; shift ;;
        --max-hours) MAXH=${2:?}; shift ;;
        *) echo "unknown option: $1"; exit 1 ;;
    esac
    shift
done

if [ -n "${MHD_ROOT:-}" ]; then
    ROOT="$MHD_ROOT/n$N"
else
    LAB=${MHD_LAB:-protopapas_lab}
    for cand in "/n/netscratch/$LAB/Everyone/$USER/mhd" \
                "/n/netscratch/$LAB/Lab/$USER/mhd" \
                "/n/netscratch/$LAB/$USER/mhd" \
                "$PWD/runs"; do
        [ -d "$cand/n$N" ] && { ROOT="$cand/n$N"; break; }
    done
fi
[ -d "${ROOT:-}" ] || { echo "no such suite directory for N=$N (set MHD_ROOT)"; exit 1; }
echo "suite: $ROOT      $([ $DO_IT = 1 ] && echo 'DELETING' || echo 'dry run (pass --yes to delete)')"

start=$(date +%s)
while : ; do
    freed=0; done_n=0; busy_n=0
    for D in "$ROOT"/ms*_ma*; do
        [ -d "$D" ] || continue
        tag=$(basename "$D")
        tmax=$(awk '/^tmax/{print $2}' "$D/params.txt" 2>/dev/null)
        # Completion is judged from the LAST SNAPSHOT's time, not from
        # history.txt: history is only written every 10 steps, so its final
        # line can sit short of tmax even for a finished run.  The final
        # snapshot, by construction, is written exactly at tmax.
        lastmeta=$(ls "$D"/snap_*.txt 2>/dev/null | tail -1)
        tnow=$(awk '/^time/{print $2}' "$lastmeta" 2>/dev/null)
        if [ -z "$tmax" ] || [ -z "$tnow" ]; then
            busy_n=$((busy_n+1)); continue
        fi
        # awk, not bc: history.txt writes times as %.8e and bc cannot parse
        # scientific notation (it errors out and the test silently passes)
        if awk -v a="$tnow" -v b="$tmax" 'BEGIN{exit !(a < 0.999*b)}'; then
            busy_n=$((busy_n+1)); continue          # still running
        fi
        snaps=($(ls "$D"/snap_*.bin 2>/dev/null))
        [ ${#snaps[@]} -eq 0 ] && continue           # already reaped
        newest=$(ls -t "$D"/snap_*.bin | head -1)
        if [ ! -f "$D/analysis_ax2.npz" ] || [ ! -f "$D/analysis_ax0.npz" ] \
           || [ ! "$D/analysis_ax0.npz" -nt "$newest" ]; then
            echo "  analysing $tag"
            if [ $DO_IT = 1 ]; then
                # on failure: continue WITHOUT deleting anything for this run
                python3 python/process_run.py "$D" --axis 2 >/dev/null \
                    || { echo "    analysis failed, keeping snapshots"; continue; }
                python3 python/process_run.py "$D" --axis 0 >/dev/null \
                    || { echo "    analysis failed, keeping snapshots"; continue; }
            else
                echo "    (dry run: would analyse, then delete ${#snaps[@]} snapshots)"
                continue
            fi
        fi
        keep=""
        [ $KEEP_LAST = 1 ] && keep=$(ls "$D"/snap_*.bin | tail -1)
        mb=0
        for f in "${snaps[@]}"; do
            [ "$f" = "$keep" ] && continue
            sz=$(stat -c%s "$f"); mb=$((mb + sz/1048576))
            [ $DO_IT = 1 ] && rm -f "$f"
        done
        freed=$((freed + mb)); done_n=$((done_n+1))
        echo "  $tag: $([ $DO_IT = 1 ] && echo freed || echo would free) ${mb} MB"
    done
    echo "== $done_n reaped, $busy_n still running, ${freed} MB $([ $DO_IT = 1 ] && echo freed || echo reclaimable)  ($(date))"
    [ "$WATCH" = "0" ] && break
    if [ "$MAXH" != "0" ]; then
        el=$(( ($(date +%s) - start) / 3600 ))
        [ "$el" -ge "$MAXH" ] && { echo "reached --max-hours $MAXH, stopping"; break; }
    fi
    [ "$busy_n" = "0" ] && { echo "nothing left running, stopping"; break; }
    sleep "$WATCH"
done
