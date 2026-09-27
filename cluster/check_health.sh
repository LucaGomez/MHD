#!/bin/bash
# Report card for a suite: progress, completion, analysis and numerical health
# of every run, robust to runs whose snapshots the reaper has already deleted.
#
#   bash cluster/check_health.sh 512            (or: 128 256 512)
#
# Completion and health are read from snap_*.txt metadata, which the reaper
# keeps, not from snap_*.bin, which it deletes.  A run is:
#   done      last snapshot written at tmax
#   running   not yet at tmax (shows % of tmax reached)
#   BAD       some snapshot has a non-finite or absurd Mach number -> recompute
set -u
cd "$(dirname "$0")/.."
for N in "$@"; do
    if [ -n "${MHD_ROOT:-}" ]; then ROOT="$MHD_ROOT/n$N"; else
        LAB=${MHD_LAB:-protopapas_lab}
        for c in "/n/netscratch/$LAB/Everyone/$USER/mhd" "$PWD/runs"; do
            [ -d "$c/n$N" ] && { ROOT="$c/n$N"; break; }
        done
    fi
    [ -d "${ROOT:-}" ] || { echo "no suite for N=$N"; continue; }
    echo
    echo "=== N=$N   $ROOT"
    printf "  %-12s %9s %6s %6s %9s %7s  %s\n" run progress snaps .bin analysed nfix health
    n_done=0; n_run=0; n_bad=0; n_ana=0; bad=""; unana=""
    for D in "$ROOT"/ms*_ma*; do
        [ -d "$D" ] || continue
        tag=$(basename "$D")
        tmax=$(awk '/^tmax/{print $2}' "$D/params.txt" 2>/dev/null)
        nmeta=$(ls "$D"/snap_*.txt 2>/dev/null | wc -l)
        nbin=$(ls "$D"/snap_*.bin 2>/dev/null | wc -l)
        nfix=$(tail -1 "$D/history.txt" 2>/dev/null | awk '{print (NF>=14)?$14:"-"}')
        # progress: prefer the last snapshot's time, else the latest history line
        if [ "$nmeta" -gt 0 ]; then
            tlast=$(awk '/^time/{print $2}' "$(ls "$D"/snap_*.txt | tail -1)")
        else
            tlast=""
        fi
        thist=$(tail -1 "$D/history.txt" 2>/dev/null | awk '{print $2}')
        pct=$(awk -v a="${tlast:-0}" -v h="${thist:-0}" -v t="${tmax:-0}" \
              'BEGIN{ x=(a>h?a:h); if (t>0) printf "%.0f%%", 100*x/t; else print "?" }')
        done_flag=$(awk -v a="${tlast:-0}" -v t="${tmax:-1}" 'BEGIN{print (a>=0.999*t)?1:0}')
        nbadsnap=0
        [ "$nmeta" -gt 0 ] && nbadsnap=$(awk '/^ms /{ if ($2 != $2 || $2 ~ /nan/ || $2+0 > 50) c++ }
                                               END{print c+0}' "$D"/snap_*.txt)
        ana="-"; [ -f "$D/analysis_ax2.npz" ] && ana="yes"
        # history NaN in the last lines means it is blowing up right now
        hnan=$(tail -20 "$D/history.txt" 2>/dev/null | awk '$4 != $4 || $4 ~ /nan/ {c++} END{print c+0}')
        if [ "$nbadsnap" -gt 0 ] || [ "$hnan" -gt 0 ]; then
            health="BAD ($nbadsnap bad snapshots)"; n_bad=$((n_bad+1)); bad="$bad $tag"
        elif [ "$done_flag" = 1 ]; then
            health="done"; n_done=$((n_done+1))
            [ "$ana" = "yes" ] && n_ana=$((n_ana+1)) || unana="$unana $tag"
        else
            health="running"; n_run=$((n_run+1))
        fi
        printf "  %-12s %9s %6s %6s %9s %7s  %s\n" "$tag" "$pct" "$nmeta" "$nbin" "$ana" "$nfix" "$health"
    done
    echo "  -----"
    echo "  $n_done done ($n_ana analysed), $n_run still running, $n_bad BAD"
    if [ -n "$bad" ];   then echo "  RECOMPUTE:$bad"; fi
    if [ -n "$unana" ]; then echo "  done but not analysed yet:$unana"; fi
done
exit 0
