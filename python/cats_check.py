"""Dust statistics of a public CATS MHD snapshot (Cho-ENO, 256^3, FITS files).

    python3 python/cats_check.py DENS.fits.gz MAGX.fits.gz MAGY.fits.gz MAGZ.fits.gz [--label s1]

The CATS "machine_learning" sets store only density and magnetic field (no
velocity), which is all a dust map needs.  Without velocity the Mach numbers
are estimated indirectly:
  B0      = |<B>|: Cho-ENO runs use B0 = 1 (sub-Alfvenic, M_A~0.7) or 0.1 (M_A~2)
  dB/B0   = rms(B - <B>) / B0, which tracks M_A
  sigma   = std(rho)/<rho>, which tracks M_S (solenoidal driving: ~M_S/3)
The line of sight is taken perpendicular to the mean field, as in the paper.
"""
import argparse, os, sys
import numpy as np
from astropy.io import fits
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analysis import project, qu_to_eb, spec2d_cross, fit_powerlaw


def load(p):
    with fits.open(p) as h:
        return np.asarray(h[0].data, dtype=np.float64)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dens"); ap.add_argument("bx"); ap.add_argument("by"); ap.add_argument("bz")
    ap.add_argument("--label", default="")
    ap.add_argument("--kmin", type=float, default=3)
    ap.add_argument("--vel", nargs=3, default=None, help="VELX VELY VELZ (documented runs only)")
    ap.add_argument("--p", type=float, default=None, help="isothermal pressure = c_s^2 (from the run name)")
    a = ap.parse_args()
    rho = load(a.dens)
    B = np.array([load(a.bx), load(a.by), load(a.bz)])
    n = rho.shape[0]
    Bm = B.mean(axis=(1, 2, 3))
    B0 = float(np.linalg.norm(Bm))
    dB = float(np.sqrt(((B - Bm[:, None, None, None])**2).sum(axis=0).mean()))
    sig = float(rho.std() / rho.mean())
    sig_ln = float(np.log(rho).std())
    along = int(np.argmax(np.abs(Bm)))                 # FITS axis order is as stored
    kmax = n // 20
    out = []
    for los in range(3):
        if los == along:
            continue
        T, Q, U = project(rho, B, axis=los)
        E, Bm_ = qu_to_eb(Q, U)
        k, CEE, _ = spec2d_cross(E); _, CBB, _ = spec2d_cross(Bm_)
        _, CTT, _ = spec2d_cross(T); _, CTE, _ = spec2d_cross(T, E)
        m = (k >= a.kmin) & (k <= kmax)
        R = np.abs(np.mean(np.exp(1j * np.arctan2(U, Q))))
        out.append(dict(
            bbee=CBB[m].sum() / CEE[m].sum(),
            rte=CTE[m].sum() / np.sqrt(CTT[m].sum() * CEE[m].sum()),
            aEE=fit_powerlaw(k, CEE, a.kmin, kmax)[1],
            aBB=fit_powerlaw(k, CBB, a.kmin, kmax)[1],
            aTT=fit_powerlaw(k, CTT, a.kmin, kmax)[1],
            spread=np.degrees(0.5 * np.sqrt(-2 * np.log(max(R, 1e-12))))))
    avg = {key: np.mean([o[key] for o in out]) for key in out[0]}
    mach = ""
    if a.vel and a.p:
        v = np.array([load(f) for f in a.vel])
        vrms = float(np.sqrt((v**2).sum(axis=0).mean()))
        vmean = float(np.sqrt((v**2).sum(axis=0)).mean())
        cs = np.sqrt(a.p)
        mach = (f"M_S(rms)={vrms/cs:.2f} M_S(mean|v|)={vmean/cs:.2f} "
                f"M_A(rms)={vrms*np.sqrt(rho.mean())/B0:.2f} | ")
    print(f"{a.label:<6} {mach}n={n} <rho>={rho.mean():.3f} B0={B0:.3f} (along axis {along}) "
          f"dB/B0={dB/B0:.2f} std(rho)/<rho>={sig:.2f} std(ln rho)={sig_ln:.2f} | "
          f"window [{a.kmin:g},{kmax}]  BB/EE={avg['bbee']:.2f} r_TE={avg['rte']:.2f} "
          f"a_TT={avg['aTT']:.2f} a_EE={avg['aEE']:.2f} a_BB={avg['aBB']:.2f} "
          f"angle spread={avg['spread']:.0f} deg")


if __name__ == "__main__":
    main()
