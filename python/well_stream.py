"""Analyse The Well's MHD_256 turbulence boxes without downloading the files.

    python3 python/well_stream.py --out runs_well [--pairs Ma_2_Ms_7 ...]
                                  [--splits valid test] [--steps 0 20 40 60 80]
    python3 python/well_stream.py --out runs_well --summarize-only

Data: The Well (Ohana et al. 2024), MHD_256 = Burkhart et al. (2020, CATS)
isothermal MHD boxes, Cho & Lazarian ENO code, 256^3, solenoidal driving at
k~2.5, 10 (M_S, M_A) pairs x 10 trajectories x 100 steps (dt = 0.01), CC BY 4.0.
Each HDF5 file is read over HTTP (fsspec range requests): one state
(rho, v, B) is ~450 MB, ~1.5 min at 5 MB/s.

For every requested state the full analysis of python/analysis.py is run
(3d spectra, T/Q/U/E/B maps and their spectra) for the line of sight along the
mean field (axis 0) and perpendicular to it (axes 1 and 2), and saved to
<out>/<pair>/state_<split><traj>_<step>.npz (spectra + 256^2 maps, ~15 MB).
Existing state files are skipped, so the job can be interrupted and resumed.
--summarize-only averages them into analysis_ax{0,2}.npz, the same files our
own runs produce, so python/figures.py and compare_res.py work unchanged.

Units.  The ENO code uses rho0 = 1 and B in units where v_A = B/sqrt(rho),
like ours, but c_s = sqrt(p) varies between runs and is not stored.  The label
M_S of each file is the *mean* |v|/c_s of the CATS README; c_s is taken as the
CATS pressure (1, .32, .1, .032, .01) that best matches it, and v and B are
divided by c_s so that the saved states are in our c_s = 1 units (dust maps are
unchanged by this; M_S and M_A come out as rms quantities, like the paper's).
The labelled M_A uses the total local field; measured against the mean field,
the "Ma 2" runs have M_A ~ 8.
"""
import argparse, glob, json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analysis import analyse
from process_run import summarize

BASE = "https://sdsc-users.flatironinstitute.org/~polymathic/data/the_well/datasets/MHD_256/data"
PAIRS = [f"Ma_{a}_Ms_{s}" for a in ("0.7", "2") for s in ("0.5", "0.7", "1.5", "2", "7")]
# CATS README lists p = 1, .32, .1, .032, .01; The Well's M_S=0.5 files need a
# colder run still (p = 2), so the list is extended to cover the labels
PRESSURES = [3.2, 2.0, 1.0, 0.32, 0.1, 0.032, 0.01]
KEEP = ["n", "k3", "Crho", "Cv", "CH", "k2", "CTT", "CEE", "CBB", "CTE", "CTB", "CEB",
        "Ms", "Ma", "sigma_lnrho"]


def pick_cs(v, label_ms):
    """c_s from the CATS pressure list whose mean |v|/c_s is closest to the label."""
    vmean = float(np.sqrt((v.astype(np.float64)**2).sum(axis=0)).mean())
    p = min(PRESSURES, key=lambda p: abs(np.log(vmean / np.sqrt(p) / label_ms)))
    return float(np.sqrt(p)), vmean


def divb_check(B):
    """|div B| relative to |grad B| with the arrays taken as (x, y, z)."""
    d = lambda a, ax: np.roll(a, -1, ax) - np.roll(a, 1, ax)
    ok = np.abs(d(B[0], 0) + d(B[1], 1) + d(B[2], 2)).mean()
    swapped = np.abs(d(B[0], 2) + d(B[1], 1) + d(B[2], 0)).mean()
    return ok / swapped


def save_state(path, snap, meta):
    out = {}
    for ax in (0, 1, 2):
        r = analyse(snap, axis=ax)
        for k in KEEP:
            out[f"ax{ax}_{k}"] = r[k]
        for m in ("T", "Q", "U", "E", "B"):
            out[f"ax{ax}_map{m}"] = r["maps"][m].astype(np.float32)
    out["meta"] = json.dumps(meta)
    tmp = path + ".tmp.npz"
    np.savez_compressed(tmp, **out)
    os.replace(tmp, path)


def state_map(z, meta, ax, name):
    """A stored T/Q/U/E/B map, in the frame analysis.project now returns.

    States written before the axis=1 frame fix (Codex review round 1, F01) hold
    the axis=1 maps transposed; they are corrected here so cached states stay
    usable.  New states carry maps_frame_fixed and are returned as written."""
    a = z[f"ax{ax}_map{name}"].astype(np.float64)
    if ax == 1 and not meta.get("maps_frame_fixed"):
        a = a.T
    return a


def load_state(path, ax):
    z = np.load(path)
    meta = json.loads(str(z["meta"]))
    r = {k: z[f"ax{ax}_{k}"] for k in KEEP}
    r["n"] = int(r["n"])
    for k in ("Ms", "Ma", "sigma_lnrho"):
        r[k] = float(r[k])
    r["maps"] = {m: state_map(z, meta, ax, m) for m in ("T", "Q", "U", "E", "B")}
    r["meta"] = meta
    # c_s was chosen when the state was written; if the pressure list has since
    # been extended, rescale M_S (v and B were divided by c_s, so M_A and the
    # dust maps are unaffected)
    cs_new = min(PRESSURES, key=lambda q: abs(np.log(meta["vmean_code"] / np.sqrt(q)
                                                     / meta["label_ms"])))**0.5
    r["Ms"] *= meta["cs_code"] / cs_new
    return r


def stream(pair, splits, steps, trajs, outdir):
    import fsspec, h5py
    label_ms = float(pair.split("_Ms_")[1])
    d = os.path.join(outdir, pair)
    os.makedirs(d, exist_ok=True)
    for split in splits:
        url = f"{BASE}/{split}/MHD_{pair}.hdf5"
        todo = []
        for tr in trajs:
            for st in steps:
                path = os.path.join(d, f"state_{split}{tr}_{st:03d}.npz")
                if not os.path.exists(path):
                    todo.append((tr, st, path))
        if not todo:
            continue
        cs_run = None            # c_s is a property of the run, not of a snapshot
        with fsspec.open(url, "rb", block_size=64 * 2**20) as fh, h5py.File(fh, "r") as f:
            ntraj = f["t0_fields/density"].shape[0]
            for tr, st, path in todo:
                if tr >= ntraj:
                    continue
                t0 = time.time()
                rho = np.ascontiguousarray(f["t0_fields/density"][tr, st])
                B = np.ascontiguousarray(np.moveaxis(f["t1_fields/magnetic_field"][tr, st], -1, 0))
                v = np.ascontiguousarray(np.moveaxis(f["t1_fields/velocity"][tr, st], -1, 0))
                tread = time.time() - t0
                cs, vmean = pick_cs(v, label_ms)
                # fix it from the first snapshot of the run and keep it: an
                # isothermal run has one sound speed (Codex review round 1, F08)
                cs = cs_run = cs if cs_run is None else cs_run
                v /= cs; B /= cs
                meta = dict(pair=pair, split=split, traj=tr, step=st, cs=1.0, cs_code=cs,
                            p_code=cs**2, vmean_code=vmean, label_ms=label_ms,
                            label_ma=float(pair.split("_")[1]),
                            time=float(f["dimensions/time"][st]),
                            B0=[float(x) for x in B.mean(axis=(1, 2, 3))],
                            divb_ratio=float(divb_check(B)) if st == steps[0] else None,
                            maps_frame_fixed=True)
                snap = dict(rho=rho, v=v, H=B, meta=meta, n=rho.shape[0])
                save_state(path, snap, meta)
                r = load_state(path, 2)
                print(f"{pair} {split}{tr} step {st:3d}: c_s={cs:.3f} M_S={r['Ms']:.2f} "
                      f"M_A={r['Ma']:.2f} B0={np.round(meta['B0'], 3).tolist()} "
                      + (f"divB(as read/swapped)={meta['divb_ratio']:.3f} " if meta['divb_ratio'] else "")
                      + f"read {tread:.0f}s total {time.time()-t0:.0f}s", flush=True)
                del rho, B, v, snap


def summarize_pair(pair, outdir):
    d = os.path.join(outdir, pair)
    files = sorted(glob.glob(os.path.join(d, "state_*.npz")))
    if not files:
        return
    for ax in (2, 0):
        res = [load_state(f, ax) for f in files]
        summarize(res, os.path.join(d, f"analysis_ax{ax}.npz"),
                  label=f"{pair} ax{ax} ({len(res)} states)")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="runs_well")
    ap.add_argument("--pairs", nargs="+", default=PAIRS)
    ap.add_argument("--splits", nargs="+", default=["valid", "test"])
    ap.add_argument("--trajs", type=int, nargs="+", default=[0])
    ap.add_argument("--steps", type=int, nargs="+", default=[0, 20, 40, 60, 80])
    ap.add_argument("--summarize-only", action="store_true")
    a = ap.parse_args()
    for pair in a.pairs:
        if not a.summarize_only:
            stream(pair, a.splits, a.steps, a.trajs, a.out)
        summarize_pair(pair, a.out)


if __name__ == "__main__":
    main()
