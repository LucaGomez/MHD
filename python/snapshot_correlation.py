"""How long does a driven-turbulence box remember its state?

    python3 python/snapshot_correlation.py runs/n64/ms5_ma2 runs/n64/ms6_ma1 ...

For every pair of snapshots in the same run, and every pair of snapshots taken
from DIFFERENT runs (independent forcing seeds: the "independent" baseline),
computes the correlation of the fields in bands of k/k_min:

    r = sum Re(X1 X2*) / sqrt(sum|X1|^2 sum|X2|^2)       (mean removed)

for the projected maps T, Q, U (line of sight perpendicular to B0) and the 3D
density.  Results are binned by the time separation in units of t_dyn.
Snapshots can be treated as independent realizations once r at the scales you
care about is back at the level of the different-run baseline.
"""
import argparse, glob, os, sys, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analysis import read_snapshot, project, _fft, finite_snapshot

BANDS = [(1, 2.5, "k 1-2 (driving)"), (2.5, 8, "k 3-8"), (8, 16, "k 8-16")]


def band_masks(n, dim):
    k1 = np.fft.fftfreq(n) * n
    g = np.meshgrid(*([k1] * dim), indexing="ij")
    kk = np.sqrt(sum(x**2 for x in g))
    return [(kk >= lo) & (kk < hi) for lo, hi, _ in BANDS]


def fourier(a):
    return _fft.fftn(a - a.mean())       # threaded scipy.fft (MHD_FFT_WORKERS)


def corr(c1, c2):
    num = np.real(c1 * np.conj(c2)).sum()
    den = np.sqrt((np.abs(c1)**2).sum() * (np.abs(c2)**2).sum())
    return num / den if den > 0 else np.nan


def load_run(run, max_snaps, masks2, masks3):
    """Keep only the Fourier coefficients inside each band: a full 512^3 complex
    cube is ~2 GB per snapshot, the band coefficients a few MB."""
    tdyn = float(dict(l.split(None, 1) for l in open(os.path.join(run, "params.txt"))
                      if l.strip())["tdyn"])
    out = []
    for sf in sorted(glob.glob(os.path.join(run, "snap_*.bin")))[-max_snaps:]:
        s = read_snapshot(sf)
        if not finite_snapshot(s):
            continue
        T, Q, U = project(s["rho"], s["H"], axis=2)
        F = {}
        for name, arr, ms in (("T", T, masks2), ("Q", Q, masks2), ("U", U, masks2),
                              ("rho3d", s["rho"].astype(np.float64), masks3)):
            full = fourier(arr)
            F[name] = [full[m] for m in ms]
            del full
        out.append(dict(t=float(s["meta"]["time"]) / tdyn, F=F))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--max-snaps", type=int, default=19)
    ap.add_argument("--lags", type=float, nargs="+", default=[0.25, 0.5, 1, 1.5, 2, 3, 4])
    a = ap.parse_args()
    # all runs share one resolution; read it from the first snapshot's metadata
    first = sorted(glob.glob(os.path.join(a.runs[0], "snap_*.txt")))[0]
    n = int(float(dict(l.split(None, 1) for l in open(first) if l.strip())["n"]))
    m2, m3 = band_masks(n, 2), band_masks(n, 3)
    data = {r: load_run(r, a.max_snaps, m2, m3) for r in a.runs}
    for r in [r for r, v in data.items() if len(v) < 2]:
        print(f"skipping {os.path.basename(r.rstrip('/'))}: fewer than 2 finite snapshots")
        del data[r]
    if not data:
        sys.exit("no usable runs")
    fields = ["T", "Q", "U", "rho3d"]
    same = {}     # (field, band, lagbin) -> list of r
    diff = {}     # (field, band) -> list of r
    edges = a.lags
    def lagbin(dt):
        i = int(np.argmin([abs(dt - L) for L in edges]))
        return edges[i] if abs(dt - edges[i]) <= 0.13 * max(1, edges[i]) else None
    for r, snaps in data.items():
        for s1, s2 in itertools.combinations(snaps, 2):
            L = lagbin(abs(s2["t"] - s1["t"]))
            if L is None:
                continue
            for f in fields:
                for bi in range(len(BANDS)):
                    same.setdefault((f, bi, L), []).append(corr(s1["F"][f][bi], s2["F"][f][bi]))
    for r1, r2 in itertools.combinations(list(data), 2):
        for s1 in data[r1][-4:]:
            for s2 in data[r2][-4:]:
                for f in fields:
                    for bi in range(len(BANDS)):
                        diff.setdefault((f, bi), []).append(corr(s1["F"][f][bi], s2["F"][f][bi]))

    print("runs:", ", ".join(os.path.basename(x.rstrip("/")) for x in data), f"  (N={n})")
    print("mean correlation between snapshots of the SAME run vs separation;")
    print("last column: snapshots of DIFFERENT runs (independent seeds) = rms, the 'independent' level\n")
    for f in fields:
        print(f"== {f}")
        print(f"   {'band':<17}" + "".join(f"{L:>7g}" for L in edges) + "   t_dyn  | indep. rms")
        for bi, (_, _, name) in enumerate(BANDS):
            row = f"   {name:<17}"
            for L in edges:
                v = same.get((f, bi, L))
                row += f"{np.mean(v):7.2f}" if v else "      -"
            d = np.array(diff.get((f, bi), [np.nan]))
            row += f"          |  {np.sqrt(np.mean(d**2)):.2f}" if len(data) > 1 else "          |   -"
            print(row)
        print()


if __name__ == "__main__":
    main()
