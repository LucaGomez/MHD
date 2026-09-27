"""Write Enzo parameter files matched to the paper and to src/mhd.c.

    python python/enzo_setup.py --ms 4.7 --ma 1.5 --n 64 --seed 1 --out runs_enzo/x

Conventions shared with our solver (so both codes run the same physics):
  L = 1, rho0 = 1, c_s ~ 1, B0 along x,
  t_dyn = L0 / M_S with L0 = 1/2, driven NTDYN t_dyn, dumps every DUMP t_dyn.
  v_A = B/sqrt(rho) in code units  ->  B0_code = M_S / M_A,
  DrivenFlowMagField (Gaussian) = sqrt(4 pi) * B0_code.
"""
import argparse, math
from pathlib import Path

TEMPLATE = Path(__file__).resolve().parents[1] / "enzo" / "crosscheck.enzo.in"

# Measured for this template at the Planck point (M_S 4.7, M_A 1.5), 32^3:
#   input DrivenFlowMach 4.7  -> achieved M_S 5.30 (with measured c_s) -> x1.132
#   mean c_s over 5-10 t_dyn  = 1.10 and still rising (~1.08 -> 1.12)
# With (calib 1.132, cs 1.12) the achieved point was M_S 4.67, M_A 1.46.
# Numerical dissipation depends on N, so RE-MEASURE at the production
# resolution: enzo2snap.py prints the achieved M_S, M_A and c_s per dump.
MACH_CALIBRATION = 1.132
CS_EXPECTED = 1.12


def write(out, ms, ma, n, seed, ntdyn=10.0, dump_tdyn=0.5, calib=None, cs=None,
          va_cap=None):
    """calib: achieved/input Mach ratio of Enzo's fixed-amplitude forcing.
    cs: expected mean sound speed over the analysis window.  gamma = 1.001 is
    only nearly isothermal and Enzo heats (c_s ~ 1.1 by 5-10 t_dyn at M_S ~ 5),
    so for M_A = v_rms / v_A the field must be B0 = M_S c_s / M_A."""
    calib = MACH_CALIBRATION if calib is None else calib
    cs = CS_EXPECTED if cs is None else cs
    tdyn = 0.5 / ms
    b0_code = ms * cs / ma
    txt = (TEMPLATE.read_text()
           .replace("__N__", str(n))
           .replace("__STOP__", f"{ntdyn * tdyn:.8g}")
           .replace("__DTDUMP__", f"{dump_tdyn * tdyn:.8g}")
           .replace("__SEED__", str(int(seed)))
           .replace("__MACHIN__", f"{ms / calib:.6g}")
           .replace("__BIN__", f"{math.sqrt(4 * math.pi) * b0_code:.8g}"))
    if va_cap:
        # Alfven-speed limiter (hydro_rk/Grid_SetFloor.C): where B/sqrt(rho)
        # exceeds the cap, rho is raised to B^2/cap^2.  Only near-empty cells
        # are touched (they carry ~1e-6 of the mass, so the dust maps do not
        # change), but without it their Alfven speed sets the time step of the
        # whole box and the run crawls (512^3: 131k cycles per t_dyn at 2.4).
        txt += (f"\n# Alfven-speed limiter, cap = {va_cap:g} code units\n"
                f"UseFloor                   = 1\n"
                f"MaximumAlvenSpeed          = {va_cap:g}\n")
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    (out / "turbulence.enzo").write_text(txt)
    (out / "setup.txt").write_text(
        f"ms_target {ms}\nma_target {ma}\nn {n}\nseed {seed}\ntdyn {tdyn}\n"
        f"ntdyn {ntdyn}\ndump_tdyn {dump_tdyn}\nmach_calibration {calib}\n"
        f"b0_code {b0_code}\ncs_expected {cs}\nva_cap {va_cap or 0}\n")
    return tdyn


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ms", type=float, required=True)
    ap.add_argument("--ma", type=float, required=True)
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--ntdyn", type=float, default=10.0)
    ap.add_argument("--dump-tdyn", type=float, default=0.5)
    ap.add_argument("--calib", type=float, default=None)
    ap.add_argument("--cs", type=float, default=None,
                    help="expected mean sound speed in the analysis window")
    ap.add_argument("--va-cap", type=float, default=None,
                    help="cap on the Alfven speed B/sqrt(rho), code units (UseFloor)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    t = write(a.out, a.ms, a.ma, a.n, a.seed, a.ntdyn, a.dump_tdyn, a.calib, a.cs,
              va_cap=a.va_cap)
    print(f"wrote {a.out}/turbulence.enzo  (t_dyn = {t:.6g})")
