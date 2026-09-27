"""How much the dust statistics depend on the line-of-sight integration depth.

    python3 python/los_depth.py runs/n128/ms5_ma2 runs/n128/ms6_ma1 ... [--depths 128 64 32 16 8]

For every snapshot, the cube is cut along the line of sight (axis 2, perpendicular
to the mean field) into non-overlapping slabs of each depth; each slab is
projected over the FULL transverse sky (so E/B stay clean and periodic), and
the statistics are averaged over slabs and snapshots:

  BB/EE   band power ratio over k in [kmin, kmax]
  r_TE    T-E correlation over the same band
  p_med   median polarization fraction sqrt(Q^2+U^2)/T (no p0 factor)
  a_EE    E-mode spectral slope over the band
"""
import argparse, glob, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analysis import read_snapshot, project, qu_to_eb, spec2d_cross, fit_powerlaw, finite_snapshot


def stats(rho, H, kmin, kmax):
    T, Q, U = project(rho, H, axis=2)
    E, B = qu_to_eb(Q, U)
    k, CEE, _ = spec2d_cross(E); _, CBB, _ = spec2d_cross(B)
    _, CTT, _ = spec2d_cross(T); _, CTE, _ = spec2d_cross(T, E)
    m = (k >= kmin) & (k <= kmax)
    return dict(bbee=CBB[m].sum() / CEE[m].sum(),
                rte=CTE[m].sum() / np.sqrt(CTT[m].sum() * CEE[m].sum()),
                p=float(np.median(np.sqrt(Q**2 + U**2) / T)),
                aee=fit_powerlaw(k, CEE, kmin, kmax)[1])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--depths", type=int, nargs="+", default=None)
    ap.add_argument("--kmin", type=float, default=3)
    ap.add_argument("--kmax", type=float, default=None)
    ap.add_argument("--max-snaps", type=int, default=6)
    a = ap.parse_args()
    acc, meta = {}, []
    for run in a.runs:
        snaps = sorted(glob.glob(os.path.join(run, "snap_*.bin")))[-a.max_snaps:]
        used = 0
        for sf in snaps:
            s = read_snapshot(sf)
            if not finite_snapshot(s):
                continue
            used += 1
            n = s["n"]
            depths = a.depths or [n, n // 2, n // 4, n // 8, n // 16]
            kmax = a.kmax or n // 10
            for d in depths:
                for k0 in range(0, n - d + 1, d):
                    r = stats(s["rho"][:, :, k0:k0 + d], s["H"][:, :, :, k0:k0 + d], a.kmin, kmax)
                    for key, v in r.items():
                        acc.setdefault((d, key), []).append(v)
        if not used:
            meta.append(f"{os.path.basename(run.rstrip('/'))} (SKIPPED: no finite snapshot)")
            continue
        meta.append(f"{os.path.basename(run.rstrip('/'))} (M_S {float(s['meta'].get('ms', 'nan')):.2f}, "
                    f"M_A {float(s['meta'].get('ma', 'nan')):.2f}, {used} snaps)")
    print("runs:", "; ".join(meta))
    if not acc:
        sys.exit("no finite snapshots in any run")
    print(f"box N={n}, band k=[{a.kmin:g},{kmax:g}], line of sight perpendicular to B0\n")
    print(f"  {'depth':>6} {'% of box':>8} {'slabs':>6} {'BB/EE':>13} {'r_TE':>13} {'p median':>13} {'alpha_EE':>13}")
    for d in depths:
        cnt = len(acc[(d, 'bbee')])
        row = f"  {d:6d} {100*d/n:7.1f}% {cnt:6d}"
        for key in ("bbee", "rte", "p", "aee"):
            v = np.array(acc[(d, key)])
            row += f" {np.mean(v):6.3f}±{np.std(v):<5.3f}"
        print(row)
    print("\nPaper/Planck at M_S>4, full-box integration: BB/EE ~0.55 (Planck 0.53), r_TE ~0.3 (0.355), alpha_EE -2.42")


if __name__ == "__main__":
    main()
