#!/bin/bash
# Make Enzo's Alfven-speed limiter usable at scale.
#
#   bash cluster/patch_enzo_floor.sh        (login node is fine: one file + relink)
#
# UseFloor = 1 + MaximumAlvenSpeed = X (hydro_rk/Grid_SetFloor.C) raises the
# density of cells whose Alfven speed B/sqrt(rho) exceeds X.  Stock Enzo prints
# one line per such cell per step -- at 256^3 that is gigabytes of enzo.log.
# This wraps that printf in `if (debug)` and rebuilds.  Without UseFloor in the
# parameter file the binary behaves exactly as before.
set -eu
cd "$(dirname "$0")/.."
set +u; source cluster/enzo_env.sh; set -u
SRC=$(dirname "$ENZO")
F=$SRC/hydro_rk/Grid_SetFloor.C
[ -f "$F" ] || { echo "no $F"; exit 1; }
if grep -q 'if (debug) printf("floor set based on MaximumAlvenSpeed' "$F"; then
    echo "already patched: $F"
else
    sed -i 's|^\(\s*\)printf("floor set based on MaximumAlvenSpeed|\1if (debug) printf("floor set based on MaximumAlvenSpeed|' "$F"
    grep -q 'if (debug) printf("floor set based on MaximumAlvenSpeed' "$F" || { echo "patch did not apply"; exit 1; }
    echo "patched: $F"
fi
cd "$SRC"
make -j8 > build_floor.log 2>&1 || { tail -20 build_floor.log; exit 1; }
ls -la --time-style=+%F_%T enzo.exe hydro_rk/Grid_SetFloor.o
echo "rebuilt $ENZO"
