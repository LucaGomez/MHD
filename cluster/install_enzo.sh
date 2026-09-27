#!/bin/bash
# Build Enzo (MPI) on FASRC for the crosscheck and the production ensemble.
# Run once, on a login node, from ~/MHD:
#
#     bash cluster/install_enzo.sh
#
# Enzo has no OpenMP: at 256^3-512^3 it MUST be an MPI build, which is why
# this is not the serial build used on the laptop.
#
# The compile logic is cluster/install_enzo_core.sh, the script that already
# built and ran Enzo successfully (gfortran >= 10 flags, machine-file quirks,
# opt-high, RAM-aware -j).  This wrapper only loads the FASRC toolchain.
#
# Module names: override with ENZO_MODULES if these don't resolve, e.g.
#     module spider openmpi        # see what exists
#     ENZO_MODULES="gcc/13.2.0-fasrc01 openmpi/5.0.2-fasrc01" bash cluster/install_enzo.sh
#
# HDF5 must be SERIAL.  With no HDF5_ROOT set, the core script finds one
# ($HDF5_HOME, h5cc on PATH, or the system /usr/include/hdf5.h).  The system
# HDF5 on FASRC Rocky nodes is serial and works; no conda env is needed here.
set -euo pipefail
cd "$(dirname "$0")/.."

ENZO_MODULES="${ENZO_MODULES:-gcc openmpi}"
echo "== loading modules: $ENZO_MODULES"
# shellcheck disable=SC1091
command -v module >/dev/null 2>&1 || source /etc/profile.d/modules.sh 2>/dev/null || true
if ! module load $ENZO_MODULES; then
    echo "module load failed.  Find the right names with:  module spider openmpi"
    echo "then:  ENZO_MODULES=\"gcc/<ver> openmpi/<ver>\" bash cluster/install_enzo.sh"
    exit 1
fi
module list 2>&1 | sed 's/^/   /'

if [ -z "${HDF5_ROOT:-}" ] && [ -n "${CONDA_PREFIX:-}" ] && [ -f "$CONDA_PREFIX/include/hdf5.h" ]; then
    export HDF5_ROOT="$CONDA_PREFIX"
fi
echo "== HDF5_ROOT=${HDF5_ROOT:-<auto-detect>}"

PREFIX="${PREFIX:-$HOME/soft}"
ENZO_DIR="$PREFIX/enzo-dev/src/enzo"
ENV_FILE="${ENZO_ENV_FILE:-cluster/enzo_env.sh}"

# --env-only: Enzo is already compiled; just (re)write and check enzo_env.sh
if [ "${1:-}" != "--env-only" ]; then
    PREFIX="$PREFIX" NPROC="${NPROC:-8}" USE_MPI=yes bash cluster/install_enzo_core.sh
fi
[ -x "$ENZO_DIR/enzo.exe" ] || { echo "no $ENZO_DIR/enzo.exe -- the build did not finish"; exit 1; }

# The HDF5 Enzo was actually built against is recorded in its machine file
# (the core script may have auto-detected it, e.g. /usr).
H5=$(awk '/^LOCAL_HDF5_INSTALL/ {print $3; exit}' "$ENZO_DIR/Make.mach.cluster")
H5LIB=""
for d in "$H5/lib64" "$H5/lib"; do
    ls "$d"/libhdf5.so* >/dev/null 2>&1 && { H5LIB="$d"; break; }
done
# System directories are on the default search path already.  Putting one of
# them in FRONT of LD_LIBRARY_PATH makes the loader pick the system's old
# libstdc++/libgfortran over the compiler's -- exactly the failure seen on
# FASRC:  /usr/lib64/libstdc++.so.6: version `GLIBCXX_3.4.32' not found
is_system_dir() { case "$1" in /usr/lib|/usr/lib64|/lib|/lib64|/usr/lib/x86_64-linux-gnu|/lib/x86_64-linux-gnu) return 0 ;; esac; return 1; }

# Runtime libraries of the compiler and MPI that built enzo.exe, recorded
# explicitly so the job works even if module state differs in batch.
libdir_of() { local f; f=$("$1" -print-file-name="$2" 2>/dev/null); [ -f "$f" ] && dirname "$(readlink -f "$f")"; }
GCCLIB=$(libdir_of gcc libstdc++.so)
GFLIB=$(libdir_of gfortran libgfortran.so)
MPILIB="$(dirname "$(dirname "$(command -v mpicc)")")/lib"
[ -d "$MPILIB" ] || MPILIB=""
LIBPATH=""
for d in "$GCCLIB" "$GFLIB" "$MPILIB" "$H5LIB"; do
    [ -n "$d" ] || continue
    is_system_dir "$d" && continue      # default search path already; never in front
    case ":$LIBPATH:" in *":$d:"*) ;; *) LIBPATH="${LIBPATH:+$LIBPATH:}$d" ;; esac
done
echo "== runtime library path for enzo.exe: $LIBPATH"

# Pin the exact versions loaded now, so a later change of cluster defaults
# can't silently load a compiler/MPI that doesn't match the binary.
PINNED=$(module -t list 2>&1 | grep -E '^(gcc|openmpi)/' | tr '\n' ' ')
[ -n "$PINNED" ] || PINNED="$ENZO_MODULES"

cat > "$ENV_FILE" <<ENV
# written by cluster/install_enzo.sh on $(date)
# Safe to source from scripts running with 'set -u' (Lmod is not).
case \$- in *u*) _enzo_env_nounset=1; set +u ;; esac
module load $PINNED
export ENZO=$ENZO_DIR/enzo.exe
# compiler + MPI runtime of the build first, never a system directory
export LD_LIBRARY_PATH=$LIBPATH\${LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH}
if [ -n "\${_enzo_env_nounset:-}" ]; then set -u; fi
unset _enzo_env_nounset
ENV
echo
echo "== wrote $ENV_FILE; job scripts source it."
cat "$ENV_FILE"

# Check the way a batch job sees it: modules purged, then only the env file.
# Any unresolved library or symbol version is a failure.
CHECK=$( ( module purge >/dev/null 2>&1; set +u; source "$ENV_FILE"; "$ENZO_DIR/enzo.exe" -h 2>&1; ldd "$ENZO_DIR/enzo.exe" 2>&1 ) )
if echo "$CHECK" | grep -qiE "not found|error while loading|cannot open shared"; then
    echo "$CHECK" | grep -iE "not found|error while loading|cannot open shared" | sort -u | head -8
    echo "!! enzo.exe does not start with a clean module environment."
    exit 1
fi
set +u
# shellcheck disable=SC1090
source "$ENV_FILE"
set -u
echo "== enzo.exe starts OK"

if python3 -c "import numpy, scipy, h5py, matplotlib" 2>/dev/null; then
    echo "== python: numpy/scipy/h5py/matplotlib present"
else
    echo "!! the post-processing needs numpy, scipy, h5py and matplotlib in venv_gpu:"
    echo "   conda activate venv_gpu && pip install h5py scipy matplotlib"
fi
echo
echo "Next: smoke test on the test partition (see the header of cluster/enzo_box.slurm)"
