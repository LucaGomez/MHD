"""Does Enzo's Alfven-speed limiter change the turbulence?

    python3 python/floor_test.py runs_enzo/floor_test --tmin 2 --out figures/floor_test

Two 64^3 Enzo runs with identical parameters and forcing seed, one with
UseFloor = 1 / MaximumAlvenSpeed = 50 and one without.  In the rare cells where
B/sqrt(rho) exceeds the cap the limiter raises rho to B^2/cap^2; those cells
hold ~1e-6 of the mass but otherwise set the time step of the whole box.

Turbulence is chaotic, so the two runs decorrelate a few t_dyn after the
limiter first fires: the comparison has to be statistical.  For every dump
after TMIN t_dyn this prints the cost (cycles), the mass added, and the dust
statistics, then the mean and scatter over dumps for each run.
"""
import argparse, glob, os, re, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from enzo2snap import read_dump, dump_params
from analysis import project, qu_to_eb, spec2d_cross, spec3d, fit_powerlaw

KEYS = ["Ms", "Ma", "sigma_lnrho", "aTT", "aEE", "aBB", "bbee", "rte", "arho", "av", "aH"]


def dumps(run):
    for d in sorted(glob.glob(os.path.join(run, "DD[0-9]*"))):
        p = os.path.join(d, "data" + d.rsplit("DD", 1)[1])
        if os.path.exists(p):
            yield p


def cycles(run, tdyn):
    pat = re.compile(r"TopGrid dt = ([0-9.eE+-]+)\s+time = ([0-9.eE+-]+)")
    t = [float(m.group(2)) / tdyn for l in open(os.path.join(run, "enzo.log"), errors="replace")
         for m in [pat.search(l)] if m]
    return np.array(t)


def stats(path, kmin, kmax):
    f, par = read_dump(path)
    # read_dump returns [z, y, x]; transpose to [x, y, z] as enzo2snap does
    rho = np.ascontiguousarray(f["rho"].T)
    v = np.array([f["vx"].T, f["vy"].T, f["vz"].T])
    B = np.array([f["bx"].T, f["by"].T, f["bz"].T])
    cs = float(np.sqrt((rho * 1.001 * 0.001 * f["GasEnergy"].T).sum() / rho.sum())) \
        if "GasEnergy" in f else 1.0
    vrms = float(np.sqrt((v**2).sum(axis=0).mean()))
    b0 = float(np.linalg.norm(B.mean(axis=(1, 2, 3))))
    out = dict(t=float(par["InitialTime"]), Ms=vrms / cs, Ma=vrms / (b0 / np.sqrt(rho.mean())),
               sigma_lnrho=float(np.std(np.log(rho))), mass=float(rho.mean()),
               rho_min=float(rho.min()), va_max=float(np.sqrt((B**2).sum(axis=0) / rho).max()))
    k3, Crho, _ = spec3d([rho]); _, Cv, _ = spec3d(list(v)); _, CH, _ = spec3d(list(B))
    for nm, C in (("arho", Crho), ("av", Cv), ("aH", CH)):
        out[nm] = fit_powerlaw(k3, C, kmin, kmax)[1]
    T, Q, U = project(rho, B, axis=2)
    E, Bm = qu_to_eb(Q, U)
    k, CTT, _ = spec2d_cross(T); _, CEE, _ = spec2d_cross(E); _, CBB, _ = spec2d_cross(Bm)
    _, CTE, _ = spec2d_cross(T, E)
    m = (k >= kmin) & (k <= kmax)
    out.update(aTT=fit_powerlaw(k, CTT, kmin, kmax)[1], aEE=fit_powerlaw(k, CEE, kmin, kmax)[1],
               aBB=fit_powerlaw(k, CBB, kmin, kmax)[1], bbee=CBB[m].sum() / CEE[m].sum(),
               rte=CTE[m].sum() / np.sqrt(CTT[m].sum() * CEE[m].sum()))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", nargs="?", default="runs_enzo/floor_test")
    ap.add_argument("--runs", nargs="+", default=["nofloor", "vA50"])
    ap.add_argument("--ms", type=float, default=4.7)
    ap.add_argument("--tmin", type=float, default=2.0, help="first dump to use, in t_dyn")
    ap.add_argument("--kmin", type=float, default=3)
    ap.add_argument("--kmax", type=float, default=6)
    ap.add_argument("--out", default="figures/floor_test")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    tdyn = 0.5 / a.ms
    L, res = [], {}
    L.append(f"Enzo 64^3, M_S {a.ms}, identical forcing seed, one with the Alfven-speed "
             f"limiter (cap 50) and one without.")
    L.append(f"Dust statistics over k = [{a.kmin:g}, {a.kmax:g}], dumps from {a.tmin:g} t_dyn on.\n")
    for run in a.runs:
        d = os.path.join(a.root, run)
        t = cycles(d, tdyn)
        rows = [stats(p, a.kmin, a.kmax) for p in dumps(d)
                if float(dump_params(p)["InitialTime"]) / tdyn >= a.tmin - 1e-6]
        res[run] = rows
        L.append(f"== {run}: {len(t)} cycles to t = {t.max():.2f} t_dyn, {len(rows)} dumps used")
        L.append("   cycles per t_dyn: " + " ".join(
            f"{lo:.0f}-{lo+1:.0f}: {int(((t >= lo) & (t < lo + 1)).sum())}"
            for lo in range(0, int(t.max()) + 1)))
        L.append(f"   lowest density {min(r['rho_min'] for r in rows):.2e}, "
                 f"highest Alfven speed {max(r['va_max'] for r in rows):.0f}, "
                 f"mean density {np.mean([r['mass'] for r in rows]):.6f}")
    if len(a.runs) == 2 and all(res[r] for r in a.runs):
        m0 = np.mean([r["mass"] for r in res[a.runs[0]]])
        m1 = np.mean([r["mass"] for r in res[a.runs[1]]])
        L.append(f"\nmass added by the limiter: {(m1 / m0 - 1):.2e} of the total")
    L.append("\n" + f"{'quantity':<14}" + "".join(f"{r:>26}" for r in a.runs))
    for k in KEYS:
        row = f"{k:<14}"
        for r in a.runs:
            v = np.array([x[k] for x in res[r]])
            row += f"{v.mean():12.3f} +- {v.std():<11.3f}"
        L.append(row)
    L.append("\nThe two runs are different realizations once the limiter fires (turbulence is "
             "chaotic),\nso agreement within the scatter over dumps is the strongest statement "
             "available.")
    txt = "\n".join(L)
    print(txt)
    open(os.path.join(a.out, "table.txt"), "w").write(txt + "\n")


if __name__ == "__main__":
    main()
