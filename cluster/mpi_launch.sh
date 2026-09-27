# Sourced by the Enzo job scripts.  Launch Enzo on N MPI ranks and verify it.
#
# Why this exists: on FASRC, OpenMPI 5 started with a plain `srun -n 8` does
# NOT form one MPI job -- each of the 8 processes initializes as its own
# 1-process run (Enzo's performance.out said "MPI processes: 1").  One copy
# does all the work at one-core speed and the others die fighting over the
# same output files.  The ranks need Slurm's PMIx plugin (srun --mpi=pmix),
# or OpenMPI's own mpirun inside the allocation.
#
#   MPI_LAUNCH=auto (default) | srun-pmix | mpirun | srun
#   ENZO_SRUN_EXTRA  extra srun flags, e.g. --cpu-bind=verbose,cores

enzo_launcher() {          # prints the launcher that enzo_launch will use
    local mode=${MPI_LAUNCH:-auto}
    if [ "$mode" = auto ]; then
        if command -v srun >/dev/null 2>&1 && srun --mpi=list 2>&1 | grep -qi pmix; then
            mode=srun-pmix
        elif command -v mpirun >/dev/null 2>&1; then
            mode=mpirun
        else
            mode=srun
        fi
    fi
    echo "$mode"
}

enzo_launch() {            # enzo_launch NRANKS enzo-args...
    local nr=$1; shift
    if [ "$nr" -le 1 ] || [ -z "${SLURM_JOB_ID:-}" ]; then
        "$ENZO" "$@"; return
    fi
    # one core per rank; drop CPU settings inherited from an interactive salloc
    unset SLURM_CPUS_PER_TASK SLURM_TRES_PER_TASK
    case "$(enzo_launcher)" in
        srun-pmix) srun --mpi=pmix -n "$nr" --cpus-per-task=1 ${ENZO_SRUN_EXTRA:-} "$ENZO" "$@" ;;
        mpirun)    mpirun -np "$nr" --bind-to core "$ENZO" "$@" ;;
        *)         srun -n "$nr" --cpus-per-task=1 ${ENZO_SRUN_EXTRA:-} "$ENZO" "$@" ;;
    esac
}

enzo_nproc() {             # MPI size Enzo reports in <run>/performance.out (last start); 0 if none
    local n
    n=$(awk -F: '/MPI processes:/ {gsub(/ /, "", $2); n = $2} END {print n + 0}' "$1" 2>/dev/null)
    echo "${n:-0}"
}
