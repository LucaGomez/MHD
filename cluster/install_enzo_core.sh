#!/usr/bin/env bash
#
# Compila Enzo para el ensamble de turbulencia MHD.
#
#   bash install/install_enzo.sh                       # con los defaults
#   PREFIX=~/soft NPROC=16 bash install/install_enzo.sh
#   USE_MPI=no bash install/install_enzo.sh            # build serial
#
# Todo lo que hace esta script salio de compilar y correr Enzo de verdad en
# una maquina Ubuntu 24.04 con gcc 13 -- incluidos los tres tropiezos que
# hubo en el camino, que estan resueltos abajo y marcados con [1] [2] [3].
#
# Al terminar imprime las dos lineas de entorno que hay que exportar para
# poder correr enzo.exe.
#
set -euo pipefail

PREFIX="${PREFIX:-$HOME/soft}"
NPROC="${NPROC:-4}"
USE_MPI="${USE_MPI:-yes}"
ENZO_BRANCH="${ENZO_BRANCH:-main}"
SRC="${PREFIX}/enzo-dev"

log() { printf '\n\033[1m==> %s\033[0m\n' "$*"; }
die() { printf '\n\033[31mERROR: %s\033[0m\n' "$*" >&2; exit 1; }

log "Prefijo ${PREFIX}, ${NPROC} procesos, MPI=${USE_MPI}"
mkdir -p "${PREFIX}"

# ---------------------------------------------------------------------------
# 1. Dependencias
# ---------------------------------------------------------------------------
# En un cluster con modules, cargalos ANTES de correr esto, p.ej.:
#     module load gcc openmpi hdf5
# La script respeta lo que ya este cargado y solo completa lo que falte.

log "Buscando compiladores y HDF5"
for c in gcc g++ gfortran make git; do
    command -v "$c" >/dev/null || die "falta '$c' (module load gcc? apt install build-essential gfortran?)"
done
echo "  gcc      $(gcc -dumpversion)"
echo "  gfortran $(gfortran -dumpversion)"

if [ "${USE_MPI}" = "yes" ]; then
    command -v mpicc >/dev/null || die "falta mpicc (module load openmpi / mpich), o corre con USE_MPI=no"
    echo "  mpicc    $(mpicc -show 2>/dev/null | head -c 60)"
fi

# HDF5 **serial** (Enzo no quiere el paralelo). Se busca en este orden:
#   $HDF5_HOME / $HDF5_DIR (lo que suelen exportar los modules)
#   h5cc en el PATH
#   la instalacion de conda
HDF5_ROOT="${HDF5_ROOT:-${HDF5_HOME:-${HDF5_DIR:-}}}"
if [ -z "${HDF5_ROOT}" ] && command -v h5cc >/dev/null; then
    HDF5_ROOT="$(dirname "$(dirname "$(command -v h5cc)")")"
fi
if [ -z "${HDF5_ROOT}" ] && [ -n "${CONDA_PREFIX:-}" ]; then
    HDF5_ROOT="${CONDA_PREFIX}"
fi
[ -f "${HDF5_ROOT:-/nonexistent}/include/hdf5.h" ] \
    || die "no encuentro hdf5.h bajo '${HDF5_ROOT:-<vacio>}'. Exporta HDF5_ROOT=/ruta/a/hdf5-serial"
echo "  hdf5     ${HDF5_ROOT}"

# ---------------------------------------------------------------------------
# 2. Fuente
# ---------------------------------------------------------------------------
if [ -d "${SRC}/.git" ]; then
    log "Enzo ya estaba en ${SRC}, lo reuso"
else
    log "Clonando enzo-dev"
    git clone --depth 1 --branch "${ENZO_BRANCH}" https://github.com/enzo-project/enzo-dev.git "${SRC}"
fi

# ---------------------------------------------------------------------------
# 3. Machine file
# ---------------------------------------------------------------------------
log "Escribiendo Make.mach.cluster"
cd "${SRC}/src/enzo"

# [1] El machine file de ejemplo tiene un MACH_TEXT de varias lineas unidas
#     con '\'. Si lo reemplazas con sed linea por linea, las continuaciones
#     quedan huerfanas y make tira "missing separator" en la linea 29. Por
#     eso aca solo se tocan lineas que son de una sola linea.
#
# [2] gfortran >= 10 rechaza codigo F77 viejo que pasa argumentos con tipos
#     inconsistentes, que Enzo tiene a montones. Hacen falta
#     -fallow-argument-mismatch y -fallow-invalid-boz.
sed -e "s|^LOCAL_HDF5_INSTALL .*|LOCAL_HDF5_INSTALL    = ${HDF5_ROOT}|" \
    -e "s|^MACH_FILE  = Make.mach.linux-gnu|MACH_FILE  = Make.mach.cluster|" \
    -e "s|^MACH_FFLAGS   = .*|MACH_FFLAGS   = -fno-second-underscore -ffixed-line-length-132 -fallow-argument-mismatch -fallow-invalid-boz|" \
    -e "s|^MACH_F90FLAGS = .*|MACH_F90FLAGS = -fno-second-underscore -fallow-argument-mismatch -fallow-invalid-boz|" \
    Make.mach.linux-gnu > Make.mach.cluster
grep -q "^LOCAL_HDF5_INSTALL    = ${HDF5_ROOT}$" Make.mach.cluster || die "el sed del machine file no aplico"

# [4] GCC >= 14 turns implicit function declarations, implicit int, int
#     conversions and incompatible pointer types into hard ERRORS, and GCC 15
#     defaults C to C23.  Enzo's bundled uuid/ code (and some older C files)
#     rely on all of those -- found on FASRC, whose default gcc is 15.2:
#       uuid/gen_uuid.c: error: implicit declaration of function 'getpid'
#     Build C as gnu11 and downgrade those errors back to warnings.
GCC_MAJOR="$(gcc -dumpversion | cut -d. -f1)"
if [ "${GCC_MAJOR}" -ge 14 ]; then
    CFIX="-std=gnu11 -fpermissive -Wno-error=implicit-function-declaration -Wno-error=implicit-int -Wno-error=int-conversion -Wno-error=incompatible-pointer-types"
    # "MACH_CXXFLAGS =" has nothing after the '=' in the stock file, so match
    # any spacing rather than "= .*"
    sed -i -e "s|^MACH_CFLAGS *=.*|MACH_CFLAGS   = ${CFIX}|" \
           -e "s|^MACH_CXXFLAGS *=.*|MACH_CXXFLAGS = -fpermissive|" Make.mach.cluster
    grep -q "^MACH_CFLAGS   = -std=gnu11" Make.mach.cluster || die "no pude poner MACH_CFLAGS para gcc ${GCC_MAJOR}"
    grep -q "^MACH_CXXFLAGS = -fpermissive" Make.mach.cluster || die "no pude poner MACH_CXXFLAGS para gcc ${GCC_MAJOR}"
    echo "  gcc ${GCC_MAJOR}: C como gnu11 con -fpermissive (Enzo no compila si no)"
fi

# ---------------------------------------------------------------------------
# 4. Compilar
# ---------------------------------------------------------------------------
log "Configurando"
( cd "${SRC}" && ./configure >/dev/null )
make machine-cluster >/dev/null
make "use-mpi-${USE_MPI}" >/dev/null

# [3] El default de Enzo es opt-debug (-g), que corre como 3 veces mas
#     lento. Para produccion hay que pedir explicitamente opt-high (-O2).
make opt-high >/dev/null
make show-config | grep -E "CONFIG_USE_MPI|CONFIG_OPT" || true

# Enzo tiene fuentes C++ grandes y g++ come bastante memoria en cada una.
# Con poca RAM libre, un -j alto hace que el OOM killer mate compiladores y
# make aborta con un error poco claro. Regla practica: ~1.5 GB por proceso.
AVAIL_GB="$(awk '/MemAvailable/ {printf "%d", $2/1024/1024}' /proc/meminfo 2>/dev/null || echo 8)"
MAX_BY_MEM=$(( AVAIL_GB * 2 / 3 )); [ "${MAX_BY_MEM}" -lt 1 ] && MAX_BY_MEM=1
if [ "${NPROC}" -gt "${MAX_BY_MEM}" ]; then
    echo "  bajo -j${NPROC} a -j${MAX_BY_MEM}: hay ${AVAIL_GB} GB libres (~1.5 GB por g++)"
    NPROC="${MAX_BY_MEM}"
fi

log "Compilando con -j${NPROC} (tarda; ~10-15 min)"
BUILD_LOG="${SRC}/src/enzo/build.log"
# Sin '| tail', porque si falla queremos ver el error de verdad y no las
# ultimas tres lineas de ruido.
if ! make "-j${NPROC}" > "${BUILD_LOG}" 2>&1; then
    echo
    echo "--- ultimas 30 lineas de ${BUILD_LOG} ---"
    tail -30 "${BUILD_LOG}"
    die "fallo la compilacion (log completo en ${BUILD_LOG})"
fi
tail -2 "${BUILD_LOG}"

[ -x "${SRC}/src/enzo/enzo.exe" ] || die "la compilacion termino sin errores pero no hay enzo.exe (ver ${BUILD_LOG})"

# ---------------------------------------------------------------------------
# 5. Verificacion
# ---------------------------------------------------------------------------
log "Verificando"
export LD_LIBRARY_PATH="${HDF5_ROOT}/lib:${LD_LIBRARY_PATH:-}"
"${SRC}/src/enzo/enzo.exe" -h >/dev/null 2>&1 || true
if ! ldd "${SRC}/src/enzo/enzo.exe" | grep -q "not found"; then
    echo "  todas las bibliotecas resuelven"
else
    ldd "${SRC}/src/enzo/enzo.exe" | grep "not found"
    die "faltan bibliotecas en tiempo de ejecucion"
fi

cat <<EOF

$(printf '\033[1m')Listo.$(printf '\033[0m')  enzo.exe -> ${SRC}/src/enzo/enzo.exe

Agrega esto a tu ~/.bashrc o al sbatch (sin la segunda linea enzo.exe no
arranca: "error while loading shared libraries: libhdf5.so"):

    export ENZO=${SRC}/src/enzo/enzo.exe
    export LD_LIBRARY_PATH=${HDF5_ROOT}/lib:\$LD_LIBRARY_PATH

Siguiente paso:
    bash install/install_python.sh
EOF
