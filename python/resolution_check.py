"""Is a 256^3 box enough for 128x128 CNN training maps?

    python3 python/resolution_check.py --well runs_well --out figures/resolution

Three measurements, all from the 256^3 states of The Well (well_stream.py):

1. Where the cascade ends.  The local slope of the 3d spectra of rho, v and B
   is flat over the inertial range and steepens where the scheme's numerical
   dissipation takes over.  k_diss is the wavenumber where the local slope has
   moved more than TOL from its inertial value; the usable range is k = 3 to
   k_diss, and a 128^2 training map has to fit inside it.

2. What a 128^2 map keeps.  The same projected face is turned into a 128^2 map
   two ways -- 2x2-binned whole face (the map spans the box, so patch k = box k)
   and a native-resolution quarter tile (half the box width, so a patch k
   corresponds to box k x 2: the tile's band sits at SMALLER scales, past the
   dissipation knee) -- and the dust
   statistics of both are compared with the native 256^2 face.  Binning halves
   the pixel count; if the statistics survive it, 128^2 pixels are enough to
   carry them.

3. What that means on the sky.  A 128^2 map of angular size THETA has its band
   k = [3, k_diss] at ell = 360 k / THETA; the paper's dust measurements live at
   ell ~ 40-600, so THETA fixes whether the simulated range covers them.
"""
import argparse, glob, json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analysis import qu_to_eb, local_slope
from patch_stats import patch_stats
from well_stream import load_state, state_map

FIELDS = [("Crho", "density"), ("Cv", "velocity"), ("CH", "magnetic field")]


def bin2(a):
    n = a.shape[0] // 2
    return a.reshape(n, 2, n, 2).mean(axis=(1, 3))


def knee(k, C, kfit=(4, 10), tol=0.25):
    """First k where the local slope has steepened by tol from its value over
    kfit.  The search starts just above kfit[0], not above kfit[1], so a knee
    inside the reference window can still be found (Codex review round 1, F09).
    """
    sl = local_slope(k, C)
    m = (k >= kfit[0]) & (k <= kfit[1])
    ref = np.median(sl[m])
    for i in np.where(k > kfit[0])[0]:
        if np.isfinite(sl[i]) and sl[i] < ref - tol:
            return float(k[i]), float(ref)
    return float(k[-1]), float(ref)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--well", default="runs_well")
    ap.add_argument("--pairs", nargs="+", default=["Ma_0.7_Ms_7", "Ma_2_Ms_7",
                                                   "Ma_0.7_Ms_2", "Ma_2_Ms_2"])
    ap.add_argument("--tol", type=float, default=0.25)
    ap.add_argument("--ladder", nargs="*", default=[],
                    help="run suites at other resolutions (e.g. runs/n64 runs/n128) whose "
                         "analysis_ax2.npz files give the knee at that N")
    ap.add_argument("--theta", type=float, nargs="+", default=[20, 10, 5],
                    help="angular sizes (deg) a 128^2 map could represent")
    ap.add_argument("--out", default="figures/resolution")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    L, summary = [], {}

    # ---------------- 1. where the cascade ends
    L.append("1. END OF THE CASCADE IN A 256^3 BOX")
    L.append(f"   local slope measured against its median over k = 4-10; k_diss = first k "
             f"where it has steepened by more than {a.tol:g}\n")
    L.append(f"   {'run':<14}{'M_S':>6}{'M_A':>6}" + "".join(
        f"{n[:9]:>12}{'k_diss':>8}" for _, n in FIELDS))
    kdiss_all = []
    for pair in a.pairs:
        fs = sorted(glob.glob(os.path.join(a.well, pair, "state_*.npz")))
        if not fs:
            continue
        rs = [load_state(f, 2) for f in fs]
        row = f"   {pair:<14}{np.mean([r['Ms'] for r in rs]):6.2f}{np.mean([r['Ma'] for r in rs]):6.2f}"
        for key, _ in FIELDS:
            C = np.mean([r[key] for r in rs], axis=0)
            kd, ref = knee(rs[0]["k3"], C, tol=a.tol)
            row += f"{ref:12.2f}{kd:8.0f}"
            kdiss_all.append(kd)
        L.append(row)
    kd_typ = float(np.median(kdiss_all))
    L.append(f"\n   median k_diss = {kd_typ:.1f}: a 256^3 box is physical over k = 3-{kd_typ:.1f}, "
             f"a factor {kd_typ/3:.1f} in scale")
    L.append("\n   sensitivity of that median to the criterion (reference window, tolerance):")
    for kf in ((3, 8), (4, 10), (5, 12)):
        row = f"      reference k {kf[0]}-{kf[1]}: "
        for tol in (0.15, 0.25, 0.4):
            kk = []
            for pair in a.pairs:
                fs = sorted(glob.glob(os.path.join(a.well, pair, "state_*.npz")))
                if not fs:
                    continue
                rs = [load_state(f, 2) for f in fs]
                for key, _ in FIELDS:
                    kk.append(knee(rs[0]["k3"], np.mean([r[key] for r in rs], axis=0),
                                   kfit=kf, tol=tol)[0])
            row += f"tol {tol:.2f} -> {np.median(kk):5.1f}   "
        L.append(row)
    L.append(f"\n   a 512^3 box of the same code would reach k_diss ~ {2*kd_typ:.0f} if the knee "
             f"scales with resolution.\n   That is an extrapolation, not a measurement: no 512^3 "
             f"spectra were analysed here.  The ladder below tests the scaling where data exist.")
    summary["k_diss_256"] = kd_typ

    if a.ladder:
        L.append("\n   measured knee against resolution (our own solver, analysis_ax2.npz):")
        L.append(f"      {'suite':<22}{'runs':>6}{'N':>6}{'k_diss (median)':>18}{'k_diss/N':>12}")
        for root in a.ladder:
            ks, ns = [], []
            for f in sorted(glob.glob(os.path.join(root, "*", "analysis_ax2.npz"))):
                d = dict(np.load(f, allow_pickle=True))
                n = int(d["n"]); ns.append(n)
                for key in ("Crho", "Cv", "CH"):
                    ks.append(knee(d["k3"], d[key], tol=a.tol)[0])
            if ks:
                L.append(f"      {os.path.basename(root.rstrip('/')):<22}"
                         f"{len(ns)//3 if ns else 0:6d}{ns[0]:6d}{np.median(ks):18.1f}"
                         f"{np.median(ks)/ns[0]:12.3f}")

    # ---------------- 2. what a 128^2 map keeps
    L.append("\n2. DOES A 128^2 MAP KEEP THE STATISTICS OF THE 256^2 FACE?")
    L.append(f"   band k = [3, {kd_typ:.0f}] in units of the map side; 'face 256' is the native "
             f"projected face,\n   'face 128' the same face 2x2-binned, 'tile 128' a native-"
             f"resolution quarter, half the box width, so its band is box k "
             f"{2*3:.0f}-{2*kd_typ:.0f} -- past the knee)\n")
    keys = ["aTT", "aEE", "aBB", "bbee", "rte", "p_med", "S_med"]
    L.append(f"   {'run':<14}{'product':<10}" + "".join(f"{k:>10}" for k in keys))
    prods = {}
    for pair in a.pairs:
        fs = sorted(glob.glob(os.path.join(a.well, pair, "state_*.npz")))
        if not fs:
            continue
        acc = {p: [] for p in ("face 256", "face 128", "tile 128")}
        for f in fs:
            z = np.load(f, allow_pickle=True)
            meta = json.loads(str(z["meta"]))
            for ax in (1, 2):
                T, Q, U = (state_map(z, meta, ax, m) for m in "TQU")
                E, B = qu_to_eb(Q, U)
                acc["face 256"].append(patch_stats(T, Q, U, kmin=3, kmax=kd_typ, E=E, B=B))
                Tb, Qb, Ub = bin2(T), bin2(Q), bin2(U)
                Eb, Bb = qu_to_eb(Qb, Ub)
                acc["face 128"].append(patch_stats(Tb, Qb, Ub, kmin=3, kmax=kd_typ, E=Eb, B=Bb))
                h = T.shape[0] // 2
                for i in (0, h):
                    for j in (0, h):
                        s = (slice(i, i + h), slice(j, j + h))
                        acc["tile 128"].append(patch_stats(T[s], Q[s], U[s], kmin=3, kmax=kd_typ,
                                                           E=E[s], B=B[s]))
        for p, v in acc.items():
            L.append(f"   {pair if p == 'face 256' else '':<14}{p:<10}" + "".join(
                f"{np.median([x[k] for x in v]):10.2f}" for k in keys))
            prods.setdefault(p, []).extend(v)
    L.append("\n   medians over all runs and snapshots:")
    for p, v in prods.items():
        L.append(f"   {'':<14}{p:<10}" + "".join(f"{np.median([x[k] for x in v]):10.2f}" for k in keys))
    summary["products"] = {p: {k: float(np.median([x[k] for x in v])) for k in keys}
                           for p, v in prods.items()}
    d = {k: abs(summary["products"]["face 128"][k] - summary["products"]["face 256"][k])
         for k in keys}
    L.append(f"\n   binning 256 -> 128 pixels changes the statistics by at most "
             f"{max(d.values()):.2f} ({max(d, key=d.get)}); the slopes by "
             f"{max(d[k] for k in ('aTT', 'aEE', 'aBB')):.2f}")

    # ---------------- 3. what that means on the sky
    L.append("\n3. THE SAME BAND ON THE SKY")
    L.append(f"   a 128^2 map covering THETA degrees carries k = 3-{kd_typ:.0f} at these multipoles")
    L.append(f"   {'THETA':>8}{'ell(k=3)':>12}{'ell(k_diss)':>14}{'pixel':>10}")
    for th in a.theta:
        L.append(f"   {th:6.0f} deg{360/th*3:12.0f}{360/th*kd_typ:14.0f}{th*60/128:9.1f}'")
    L.append("\n   Planck measures dust over ell ~ 40-600 (Planck 2018 XI); our PySM3 patch "
             "comparison\n   used 20 deg maps, i.e. ell = 54-234 for this band.")
    summary["theta_bands"] = {th: [360 / th * 3, 360 / th * kd_typ] for th in a.theta}

    txt = "\n".join(L)
    print(txt)
    open(os.path.join(a.out, "table.txt"), "w").write(txt + "\n")
    json.dump(summary, open(os.path.join(a.out, "summary.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
