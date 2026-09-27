#!/bin/bash
# ---------------------------------------------------------------------------
#SBATCH --job-name=mhdturb
#SBATCH --array=0-20
#SBATCH -p sapphire
#SBATCH -c 32
#SBATCH -t 1-00:00
#SBATCH --mem=32000
#SBATCH -o slurm_logs/mhdturb_%A_%a.out
#SBATCH -e slurm_logs/mhdturb_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gomezlucajavier@gmail.com
# ---------------------------------------------------------------------------
# One array task per (M_S, M_A) pair: 21 runs of driven isothermal MHD.
#
# NO python here.  ./src/mhd is a pure C99 + OpenMP binary (libm only), so no
# module load and no conda env are needed -- that whole class of failure does
# not apply to this job.  The analysis job (cluster/analyse.slurm) is the one
# that needs python.
#
#   sbatch --export=ALL,N=128,NTDYN=12,TSNAP0=7.5,NSNAP=19 cluster/run_suite.slurm
#
# FASRC partitions (docs.rc.fas.harvard.edu/kb/running-jobs):
#   sapphire  112 cores/node, 990 GB, 3 days, NO job limit   <- use this
#   shared     48 cores/node, 172 GB, 3 days
#   test      112 cores/node, 990 GB, 12 h, max 5 jobs / 112 cores
# Array indices each count as one job, so a 21-task array does not fit in
# -p test: use --array with 3-5 indices plus NCHUNKS (below), or sapphire.
#
# The SBATCH lines above are sized for the 256^3 rung.  Override on the command
# line for the others (command line beats the file):
#
#   smoke test on -p test (3 jobs is the QOS cap there): the M_S=4 triplet
#           sbatch -p test --array=12,13,14 \
#              --export=ALL,N=128,NTDYN=12,TSNAP0=7.5,NSNAP=19 \
#              cluster/run_suite.slurm
#
#   full 21 runs where only 3 jobs may be submitted: 3 tasks x 7 runs each
#           sbatch -p test --array=0-2 \
#              --export=ALL,NCHUNKS=3,N=128,NTDYN=12,TSNAP0=7.5,NSNAP=19 \
#              cluster/run_suite.slurm
#
#   128^3:  sbatch -c 16 -t 0-08:00 --mem=8000 \
#              --export=ALL,N=128,NTDYN=12,TSNAP0=7.5,NSNAP=19 \
#              cluster/run_suite.slurm
#
#   256^3:  sbatch -c 32 -t 1-00:00 --mem=16000 \
#              --export=ALL,N=256,NTDYN=3,TSNAP0=1.0,NSNAP=11,UPSAMPLE=<root>/n128 \
#              cluster/run_suite.slurm
#
#   512^3:  sbatch -c 64 -t 2-12:00 --mem=48000 \
#              --export=ALL,N=512,NTDYN=3,TSNAP0=1.0,NSNAP=11,UPSAMPLE=<root>/n256 \
#              cluster/run_suite.slurm
#
# Why sapphire: the scheme is memory-bandwidth bound, so the 8-channel DDR5 on
# Sapphire Rapids buys as much as the extra cores do.
#
# On choosing -c: thread scaling is sublinear (bandwidth, not flops), so do NOT
# reflexively ask for -c 112 -- it wastes fairshare for little gain, and with
# --mem=48000 several 512^3 tasks share a 990 GB node happily anyway.  Measure
# it once instead:  bash cluster/bench.sh   (see that script).
#
# Memory actually used (-DUSE_FLOAT build): 128^3 0.4 GB, 256^3 2.5 GB,
# 512^3 17 GB, 1024^3 135 GB.
#
# Every run checkpoints every CHKINT steps and resumes by itself, so hitting
# the wall limit is harmless: resubmit the identical command.  Finished runs
# are skipped, partial ones continue from checkpoint.bin.
# ---------------------------------------------------------------------------
set -e
mkdir -p slurm_logs

N=${N:-128}
NTDYN=${NTDYN:-12}
TSNAP0=${TSNAP0:-7.5}
NSNAP=${NSNAP:-19}
EDOTFAC=${EDOTFAC:-0.42}
CHKINT=${CHKINT:-500}
UPSAMPLE=${UPSAMPLE:-}          # coarser suite to start from, or empty

# --- where to write ---------------------------------------------------------
# On FASRC $SCRATCH is the bare mount point /n/netscratch, which is NOT
# writable: the writable part is the lab subdirectory underneath it.  Try the
# sensible candidates in order and use the first one we can actually create.
# Override any of this with  --export=ALL,...,MHD_ROOT=/path/you/want
if [ -n "${MHD_ROOT:-}" ]; then
    ROOT="$MHD_ROOT/n$N"
elif [ -n "${ROOT:-}" ]; then
    :
else
    LAB=${MHD_LAB:-protopapas_lab}
    for cand in "/n/netscratch/$LAB/Everyone/$USER/mhd" \
                "/n/netscratch/$LAB/Lab/$USER/mhd" \
                "/n/netscratch/$LAB/$USER/mhd" \
                "${SCRATCH:-}/$LAB/$USER/mhd" \
                "${SLURM_SUBMIT_DIR:-$PWD}/runs"; do
        [ -z "$cand" ] && continue
        parent=$(dirname "$cand")
        if [ -d "$parent" ] && [ -w "$parent" ]; then ROOT="$cand/n$N"; break; fi
        if mkdir -p "$cand" 2>/dev/null; then ROOT="$cand/n$N"; break; fi
    done
fi

if [ -z "${ROOT:-}" ] || ! mkdir -p "$ROOT" 2>/dev/null; then
    echo "ERROR: could not create an output directory."
    echo "  tried lab '${MHD_LAB:-protopapas_lab}' under \$SCRATCH='${SCRATCH:-unset}'"
    echo "  Find your writable scratch space, e.g.:"
    echo "      ls -ld /n/netscratch/*/Everyone/\$USER /n/netscratch/*/Lab/\$USER 2>/dev/null"
    echo "  then resubmit with, for example:"
    echo "      sbatch --export=ALL,N=$N,...,MHD_ROOT=/n/netscratch/<lab>/Everyone/$USER/mhd \\"
    echo "             cluster/run_suite.slurm"
    exit 1
fi
echo "output root: $ROOT"

# Storage needed for a full 21-run rung (snapshots are float32, 7 fields):
#   128^3, 19 snaps :   23 GB        256^3, 11 snaps : 108 GB
#   512^3, 11 snaps :  870 GB        <-- netscratch, not a lab home directory
avail=$(df -Pk "$ROOT" 2>/dev/null | awk 'NR==2{print int($4/1048576)}')
[ -n "$avail" ] && echo "space available at output root: ${avail} GB"

MS_LIST=(0.5 0.5 0.5 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6)
MA_LIST=(0.5 1   2   0.5 1 2 0.5 1 2 0.5 1 2 0.5 1 2 0.5 1 2 0.5 1 2)

# Which of the 21 (M_S, M_A) pairs does this array task own?
#
#   default            : one pair per task, needs --array=0-20
#   NCHUNKS=<n> set    : task k runs pairs k, k+n, k+2n, ... sequentially,
#                        so the whole suite fits in n array tasks.  Use this
#                        where the QOS caps how many jobs you may submit
#                        (-p test allows 3):
#                            sbatch --array=0-2 --export=ALL,NCHUNKS=3,...
#                        Striding by n means each task gets one M_A value
#                        across all M_S, so the tasks cost about the same.
k=${SLURM_ARRAY_TASK_ID:-0}
IDX=()
if [ -n "${NCHUNKS:-}" ]; then
    for ((j=k; j<${#MS_LIST[@]}; j+=NCHUNKS)); do IDX+=($j); done
else
    IDX=($k)
fi
echo "array task $k handles pair index/indices: ${IDX[*]}"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export OMP_PROC_BIND=close
export OMP_PLACES=cores

if [ ! -x ./src/mhd ]; then
    echo "src/mhd not built -- run 'bash cluster/setup.sh' first"; exit 1
fi

for i in "${IDX[@]}"; do
    MS=${MS_LIST[$i]}; MA=${MA_LIST[$i]}
    TAG=ms${MS}_ma${MA}

    ARGS=(n=$N ms=$MS ma=$MA ntdyn=$NTDYN tsnap0=$TSNAP0 nsnap=$NSNAP \
          edotfac=$EDOTFAC chkint=$CHKINT seed=$((1000+i)) out=$ROOT/$TAG)
    # optional smaller Courant number for a run that went unstable (default 0.4)
    [ -n "${CFL:-}" ] && ARGS+=(cfl=$CFL)

    if [ -n "$UPSAMPLE" ]; then
        SRC=$(ls $UPSAMPLE/$TAG/snap_*.bin 2>/dev/null | tail -1)
        if [ -z "$SRC" ]; then
            echo "no source snapshot for $TAG in $UPSAMPLE"; exit 1
        fi
        ARGS+=(upsample=$SRC)
        echo "starting $TAG from $SRC"
    fi

    echo "--- pair $i: M_S=$MS M_A=$MA  N=$N  threads=$OMP_NUM_THREADS"
    echo "    out=$ROOT/$TAG   ($(date))"

    # Run the binary directly, NOT through srun.  This is a single-task OpenMP
    # job, and recent Slurm no longer propagates --cpus-per-task to srun:
    # `srun ./mhd` would silently bind every thread to one core and run ~16x
    # slower.  A run that is already finished resumes from its checkpoint and
    # exits in seconds, so re-running a chunk after a wall-limit kill is cheap.
    ./src/mhd "${ARGS[@]}"
done
echo "array task $k finished all of: ${IDX[*]}  ($(date))"
