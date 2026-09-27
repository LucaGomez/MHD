"""Process one simulation directory: analyse every snapshot, average over
snapshots, fit slopes, and store the result as an .npz file."""
import sys, os, glob, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analysis import read_snapshot, analyse, fit_powerlaw, local_slope

# fit window: the paper uses k/kmin in [4,25] at N=512, i.e. from 4x the driving
# scale down to ~20 cells per wavelength.  We scale the dissipation end with
# resolution and keep the driving end fixed.
def fit_window(n, kdrive=3, cells_per_wave=20.0):
    return kdrive, max(kdrive + 2, int(round(n / cells_per_wave)))


def process(rundir, axis=2, kdrive=3, cpw=20.0, out=None):
    snaps = sorted(glob.glob(os.path.join(rundir, "snap_*.bin")))
    if not snaps:
        raise SystemExit(f"no snapshots in {rundir}")
    res, first = [], None
    for s in snaps:
        a = analyse(read_snapshot(s), axis=axis)
        res.append(a)
        if first is None:
            first = a
    # Guard against snapshots polluted by a numerical blow-up.  Must never
    # empty the list: a single non-finite Ms makes the median NaN, which makes
    # every comparison False and would otherwise discard everything.
    msv = np.array([r["Ms"] for r in res], dtype=float)
    # a snapshot is usable only if every derived quantity is finite -- a NaN
    # anywhere in rho, v or H otherwise leaks silently into Ma, r_TE, ...
    finite = np.array([
        np.isfinite(r["Ms"]) and np.isfinite(r["Ma"])
        and all(np.isfinite(r[k]).all()
                for k in ("Crho", "Cv", "CH", "CTT", "CEE", "CBB"))
        for r in res])
    if not finite.any():
        raise SystemExit(f"{rundir}: every snapshot has a non-finite Ms "
                         f"-- the run blew up or the files are corrupt")
    good = finite & (msv < 3 * np.median(msv[finite]))
    if not good.any():
        good = finite          # spread too wide to judge; keep what is finite
    if not good.all():
        bad = [os.path.basename(s) for s, g in zip(snaps, good) if not g]
        print(f"  !! {rundir}: discarding {int((~good).sum())} anomalous "
              f"snapshots ({', '.join(bad)})")
        res = [r for r, g in zip(res, good) if g]
    if out is None:
        out = os.path.join(rundir, f"analysis_ax{axis}.npz")
    return summarize(res, out, label=rundir, kdrive=kdrive, cpw=cpw)


def summarize(res, out, label="", kdrive=3, cpw=20.0):
    """Average per-snapshot results (dicts from analysis.analyse) and fit."""
    first = res[0]
    n = first["n"]
    k1, k2 = fit_window(n, kdrive, cpw)
    keys3 = ["Crho", "Cv", "CH"]
    keys2 = ["CTT", "CEE", "CBB", "CTE", "CTB", "CEB"]
    D = {"k3": first["k3"], "k2": first["k2"], "kfit": (k1, k2), "n": n,
         "nsnap": len(res)}
    for kk in keys3 + keys2:
        arr = np.array([r[kk] for r in res])
        D[kk] = arr.mean(axis=0)
        D[kk + "_std"] = arr.std(axis=0)
    D["Ms"] = np.mean([r["Ms"] for r in res])
    D["Ms_std"] = np.std([r["Ms"] for r in res])
    D["Ma"] = np.mean([r["Ma"] for r in res])
    D["Ma_std"] = np.std([r["Ma"] for r in res])
    D["sigma_lnrho"] = np.mean([r["sigma_lnrho"] for r in res])

    # slopes: fit each snapshot, vary the window, take mean and scatter
    kd_list = [kdrive - 1, kdrive, kdrive + 1]
    kx_list = [k2, k2 + 1, k2 + 2, k2 + 3]
    autos = ["Crho", "Cv", "CH", "CTT", "CEE", "CBB"]
    for kk, kax in [(k, "k3") for k in keys3] + [(k, "k2") for k in keys2]:
        if kk not in autos:
            continue
        A, al = [], []
        for r in res:
            for a1 in kd_list:
                for a2 in kx_list:
                    aa, ss = fit_powerlaw(r[kax], r[kk], a1, a2)
                    if np.isfinite(ss):
                        A.append(aa); al.append(ss)
        D["A_" + kk] = np.mean(A);  D["A_" + kk + "_std"] = np.std(A)
        D["al_" + kk] = np.mean(al); D["al_" + kk + "_std"] = np.std(al)

    # band-power ratios inside the fit window.  The fit amplitudes A are
    # extrapolations to k=1 and therefore sensitive to slope differences;
    # the band powers are the same quantity measured where the data live.
    mw = (first["k2"] >= k1) & (first["k2"] <= k2)
    for nm, (xx, yy) in [("EE_TT", ("CEE", "CTT")), ("BB_TT", ("CBB", "CTT")),
                         ("BB_EE", ("CBB", "CEE"))]:
        v = np.array([r[xx][mw].sum() / r[yy][mw].sum() for r in res])
        D["bp_" + nm] = v.mean(); D["bp_" + nm + "_std"] = v.std()

    # correlation coefficients averaged over the fit window
    for xy, (xx, yy) in [("TE", ("CTT", "CEE")), ("TB", ("CTT", "CBB")),
                         ("EB", ("CEE", "CBB"))]:
        rk = np.array([r["C" + xy] / np.sqrt(r[xx] * r[yy]) for r in res])
        m = (first["k2"] >= k1) & (first["k2"] <= k2)
        D["rk" + xy] = rk.mean(axis=0)
        D["r" + xy] = rk[:, m].mean(axis=1).mean()
        D["r" + xy + "_std"] = rk[:, m].mean(axis=1).std()
    # maps of the last snapshot for imaging
    for m in ["T", "E", "B"]:
        D["map" + m] = res[-1]["maps"][m]
    np.savez_compressed(out, **D)
    print(f"{label}: Ms={D['Ms']:.2f} Ma={D['Ma']:.2f} window=[{k1},{k2}] "
          f"aEE={D['al_CEE']:.2f} aBB={D['al_CBB']:.2f} "
          f"B/E={D['A_CBB']/D['A_CEE']:.2f} rTE={D['rTE']:.2f}")
    return D


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("rundirs", nargs="+")
    p.add_argument("--axis", type=int, default=2)
    p.add_argument("--kdrive", type=int, default=3)
    p.add_argument("--cpw", type=float, default=20.0)
    a = p.parse_args()
    for r in a.rundirs:
        process(r, axis=a.axis, kdrive=a.kdrive, cpw=a.cpw)
