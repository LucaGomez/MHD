#!/bin/bash
# Build and validate.  Run once, on a login node, after rsync'ing the project.
#   bash cluster/setup.sh        (or: CC=icx bash cluster/setup.sh)
set -e
cd "$(dirname "$0")/.."

# slurm_logs must exist BEFORE the first sbatch: slurm opens the -o/-e files
# itself and the job fails outright if the directory is missing.
mkdir -p slurm_logs

CC_BIN=${CC:-gcc}
echo "== building with $CC_BIN (no external libraries: C99 + OpenMP + libm)"
# -fno-finite-math-only is NOT optional: plain -ffast-math lets the compiler
# assume NaN never happens, which deletes the density-floor guard and the
# blow-up check, so one NaN silently contaminates the whole run.
# Build to a temporary name and rename into place.  Overwriting the binary
# directly fails with "Text file busy" while any job is still running it, and
# rename() is atomic: running jobs keep the old inode, new jobs get the new one.
$CC_BIN -O3 -march=native -ffast-math -fno-finite-math-only -fopenmp \
        -DUSE_FLOAT -o src/mhd.new src/mhd.c -lm
mv -f src/mhd.new src/mhd

echo
echo "== circularly polarised Alfven wave, one period (expect order ~1.97)"
OMP_NUM_THREADS=8 ./src/mhd test=cpaw

echo
echo "== divergence cleaning (expect |div B| dx/|B| ~1e-2)"
OMP_NUM_THREADS=8 ./src/mhd test=divb | tail -3

echo
echo "== isothermal shock tube vs the exact Riemann solution"
OMP_NUM_THREADS=8 ./src/mhd test=tube n=64 tend=0.1 mhd=0 out=/tmp/tube64.txt
python3 python/exact_iso_tube.py /tmp/tube64.txt || \
    echo "(needs numpy+scipy; skip if python is not set up on the login node)"

echo
echo "== E/B decomposition unit tests"
python3 python/test_eb.py || echo "(needs numpy)"

echo
echo "== ready.  Submit the first rung with:"
echo "   sbatch --export=ALL,N=128,NTDYN=12,TSNAP0=7.5,NSNAP=19 cluster/run_suite.slurm"
