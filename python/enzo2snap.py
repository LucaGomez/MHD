"""Convert Enzo data dumps into the snapshot format of src/mhd.c, so Enzo runs
go through exactly the same analysis (python/process_run.py) as our solver.

    python python/enzo2snap.py <enzo_run_dir> <out_dir> --tdyn 0.10638 --tmin-tdyn 5

Output, per kept dump: snap_NNNN.bin = 7 float32 cubes (rho, vx, vy, vz, bx,
by, bz), C order with index [x, y, z]; snap_NNNN.txt with the metadata
python/analysis.py reads (n, time, cs, ms, ma, ...); plus params.txt.

Enzo specifics handled here (each one checked against a real dump):
  * HDF5 datasets are [z, y, x]; we transpose to [x, y, z].  The mean field is
    along x, and --check-axes confirms that the density structure is indeed
    elongated along the axis we call x.
  * Only active cells are written (no ghost zones) in current Enzo, but the
    shape is checked rather than assumed.
  * BaryonFileName in the .hierarchy is relative to the run root.
  * With ParallelRootGridIO the root grid is split into several .cpuNNNN
    tiles; they are assembled by GridLeftEdge.
  * Enzo writes a final dump at StopTime that can duplicate the last scheduled
    one; exact duplicates are skipped.
  * Gamma = 1.001 is only nearly isothermal, so c_s is MEASURED from GasEnergy
    (c_s^2 = gamma (gamma-1) e, mass-weighted) instead of assumed to be 1.
  * B is stored in code units with v_A = B / sqrt(rho).  The input parameter
    DrivenFlowMagField is Gaussian (Enzo divides it by sqrt(4 pi)).
"""
import argparse, os, re, sys
from pathlib import Path
import numpy as np

FIELDS = [("rho", ("Density",)),
          ("vx", ("x-velocity",)), ("vy", ("y-velocity",)), ("vz", ("z-velocity",)),
          ("bx", ("Bx", "MagneticField_C_1")), ("by", ("By", "MagneticField_C_2")),
          ("bz", ("Bz", "MagneticField_C_3"))]


def find_dumps(run):
    out = []
    for d in sorted(Path(run).glob("DD[0-9]*")):
        for c in sorted(d.glob("data[0-9]*")):
            if c.is_file() and c.suffix == "":
                out.append(c)
    return out


def dump_params(dump):
    p = {}
    for line in Path(dump).read_text(errors="replace").splitlines():
        m = re.match(r"^\s*(\w+)\s*=\s*(.+?)\s*$", line)
        if m:
            p[m.group(1)] = m.group(2)
    return p


def parse_hierarchy(dump):
    grids, cur = [], {}
    for line in Path(str(dump) + ".hierarchy").read_text().splitlines():
        line = line.strip()
        if line.startswith("Grid = "):
            if cur:
                grids.append(cur)
            cur = {}
        for key in ("GridLeftEdge", "GridRightEdge", "GridStartIndex",
                    "GridEndIndex", "GridDimension"):
            if line.startswith(key):
                val = line.split("=", 1)[1].split()
                cur[key] = [float(v) for v in val] if "Edge" in key else [int(v) for v in val]
        if line.startswith("BaryonFileName"):
            cur["file"] = line.split("=", 1)[1].strip()
    if cur:
        grids.append(cur)
    return [g for g in grids if "file" in g]


def read_dump(dump, extra=("GasEnergy",)):
    """Return (fields dict in [x,y,z] order, params)."""
    import h5py
    dump = Path(dump)
    par = dump_params(dump)
    grids = parse_hierarchy(dump)
    ntop = [int(v) for v in par["TopGridDimensions"].split()]
    if len(set(ntop)) != 1:
        raise SystemExit(f"non-cubic top grid {ntop}")
    n = ntop[0]
    names = [(k, a) for k, a in FIELDS] + [(e, (e,)) for e in extra]
    full = {k: None for k, _ in names}
    run_root = dump.parent.parent
    for g in grids:
        rel = g["file"]
        rel = rel[2:] if rel.startswith("./") else rel
        path = next((c for c in (run_root / rel, dump.parent / Path(rel).name)
                     if c.exists()), None)
        if path is None:
            raise FileNotFoundError(g["file"])
        s, e = g["GridStartIndex"], g["GridEndIndex"]
        act = [e[i] - s[i] + 1 for i in range(3)]          # x, y, z
        off = [int(round(g["GridLeftEdge"][i] * n)) for i in range(3)]
        with h5py.File(path, "r") as f:
            grp = f[next(k for k in f if k.startswith("Grid"))]
            for k, aliases in names:
                key = next((a for a in aliases if a in grp), None)
                if key is None:
                    continue
                arr = np.asarray(grp[key])
                if arr.shape == (act[2], act[1], act[0]):
                    pass
                elif arr.shape == tuple(g["GridDimension"][::-1]):
                    arr = arr[s[2]:e[2]+1, s[1]:e[1]+1, s[0]:e[0]+1]
                else:
                    raise RuntimeError(f"{path}:{key} shape {arr.shape} matches neither "
                                       f"active {act[::-1]} nor full {g['GridDimension'][::-1]}")
                if full[k] is None:
                    full[k] = np.zeros((n, n, n), dtype=np.float64)   # [z, y, x]
                full[k][off[2]:off[2]+arr.shape[0], off[1]:off[1]+arr.shape[1],
                        off[0]:off[0]+arr.shape[2]] = arr
    missing = [k for k, _ in FIELDS if full[k] is None]
    if missing:
        raise SystemExit(f"{dump}: missing fields {missing} (is HydroMethod = 4?)")
    out = {k: np.ascontiguousarray(v.transpose(2, 1, 0)) for k, v in full.items()
           if v is not None}                                          # -> [x, y, z]
    return out, par


def divb_axis_test(bx, by, bz):
    """Definitive check of the component <-> array-axis pairing.

    div B = d(bx)/d(axis0) + d(by)/d(axis1) + d(bz)/d(axis2) is small (Dedner
    cleaning keeps it at truncation level) ONLY if array axis i really is the
    direction of component i.  A wrong transpose pairs components with the
    wrong derivatives and gives an O(|grad B|) residual.  Returns the relative
    divergence for the assumed pairing and for the reversed one.
    """
    n = bx.shape[0]
    k = 2j * np.pi * np.fft.fftfreq(n) * n
    F = [np.fft.fftn(c - c.mean()) for c in (bx, by, bz)]
    def rel(order):
        d = sum(F[c] * k.reshape([n if a == ax else 1 for a in range(3)])
                for c, ax in enumerate(order))
        g = sum(np.abs(F[c] * k.reshape([n if a == ax else 1 for a in range(3)]))**2
                for c in range(3) for ax in range(3))
        return float(np.sqrt((np.abs(d)**2).sum() / (g.sum() / 3 + 1e-300)))
    return rel((0, 1, 2)), rel((2, 1, 0))


def convert(run, out, tdyn, tmin_tdyn=5.0, ms_target=None, ma_target=None,
            check_axes=True, min_spacing_tdyn=0.25):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    dumps = find_dumps(run)
    if not dumps:
        raise SystemExit(f"no Enzo dumps in {run}")
    # Scheduled dump cadence, from the setup.txt python/enzo_setup.py writes
    # next to the run.  Dumps off that cadence are restart points (StopTime,
    # StopCPUTime) and are dropped; without setup.txt fall back to spacing.
    interval = None
    st = Path(run) / "setup.txt"
    if st.exists():
        kv = dict(l.split(None, 1) for l in st.read_text().splitlines() if l.strip())
        if "dump_tdyn" in kv:
            interval = float(kv["dump_tdyn"]) * tdyn
    kept, tprev, gamma = 0, None, None
    for d in dumps:
        par = dump_params(d)
        t = float(par["InitialTime"])
        if t < tmin_tdyn * tdyn - 1e-3 * tdyn:          # dump times carry rounding
            continue
        # Keep one snapshot per dump interval.  Enzo adds extra dumps that are
        # only restart points: one at StopTime (can coincide with the last
        # scheduled dump) and two or three a few cycles apart when it stops on
        # StopCPUTime.  They are nearly the same state, so keeping them would
        # over-weight one moment in every time average.
        if interval is not None:
            k = t / interval
            if abs(k - round(k)) > 0.01:
                print(f"  {d.parent.name}: t={t:.4f} ({t/tdyn:.3f} t_dyn) is off the "
                      f"dump cadence -- restart dump, skipped")
                continue
            if tprev is not None and round(t / interval) == round(tprev / interval):
                print(f"  {d.parent.name}: repeats the dump at t={tprev:.4f}, skipped")
                continue
        elif tprev is not None and (t - tprev) < min_spacing_tdyn * tdyn:
            print(f"  {d.parent.name}: t={t:.4f} only {(t-tprev)/tdyn:.3f} t_dyn after "
                  f"the previous kept snapshot (restart dump), skipped")
            continue
        tprev = t
        f, par = read_dump(d)
        gamma = float(par.get("Gamma", 1.001))
        rho = f["rho"]
        v2 = f["vx"]**2 + f["vy"]**2 + f["vz"]**2
        vrms = float(np.sqrt(v2.mean()))
        if "GasEnergy" in f:
            cs2 = gamma * (gamma - 1.0) * f["GasEnergy"]
            cs = float(np.sqrt((rho * cs2).sum() / rho.sum()))
        else:
            cs = 1.0
        b0 = float(np.sqrt(f["bx"].mean()**2 + f["by"].mean()**2 + f["bz"].mean()**2))
        rhobar = float(rho.mean())
        ms, ma = vrms / cs, vrms / (b0 / np.sqrt(rhobar))
        n = rho.shape[0]
        with open(out / f"snap_{kept:04d}.bin", "wb") as fh:
            for k, _ in FIELDS:
                f[k].astype(np.float32).tofile(fh)
        with open(out / f"snap_{kept:04d}.txt", "w") as fh:
            fh.write(f"n {n}\ntime {t:.8e}\nfields rho vx vy vz bx by bz\n"
                     f"dtype float32\nms {ms:.6f}\nma {ma:.6f}\nvrms {vrms:.6f}\n"
                     f"b0 {b0:.6f}\ncs {cs:.6f}\ncode enzo\nsource {d}\n")
        line = (f"  {d.parent.name} t={t:.4f} ({t/tdyn:5.2f} tdyn) -> snap_{kept:04d}"
                f"  Ms={ms:.3f} Ma={ma:.3f} cs={cs:.4f} <Bx>={f['bx'].mean():.4f}"
                f" <By>={f['by'].mean():+.1e} <Bz>={f['bz'].mean():+.1e}")
        if check_axes and f["by"].std() > 0:
            ok, bad = divb_axis_test(f["bx"], f["by"], f["bz"])
            line += f"\n      axis check: |div B| rel = {ok:.2e} as mapped, {bad:.2e} if reversed"
            if not ok < 0.3 * bad:
                raise SystemExit("axis mapping FAILED the div B test -- transpose is wrong")
            check_axes = False
        print(line)
        kept += 1
    tmax = float(dump_params(dumps[-1])["InitialTime"])
    with open(out / "params.txt", "w") as fh:
        fh.write(f"n {n}\nms_target {ms_target}\nma_target {ma_target}\n"
                 f"tdyn {tdyn}\ntmax {tmax}\ncode enzo\nsource {Path(run).resolve()}\n")
    print(f"{kept} snapshots written to {out}")
    return kept


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run"); ap.add_argument("out")
    ap.add_argument("--tdyn", type=float, required=True, help="L0/M_S = 0.5/M_S")
    ap.add_argument("--tmin-tdyn", type=float, default=5.0)
    ap.add_argument("--ms", type=float); ap.add_argument("--ma", type=float)
    ap.add_argument("--min-spacing-tdyn", type=float, default=0.25,
                    help="drop dumps closer than this to the previous kept one")
    a = ap.parse_args()
    convert(a.run, a.out, a.tdyn, a.tmin_tdyn, a.ms, a.ma,
            min_spacing_tdyn=a.min_spacing_tdyn)
