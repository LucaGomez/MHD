"""Cut 128x128 dust maps out of one simulation snapshot.

    python python/extract_maps.py <snap_NNNN.bin> --out maps.npz [--tiles 4] [--axis 2]

Integrates through the FULL box along the line of sight (the depth for which
the paper's Planck-consistent statistics hold), then tiles the N x N sky into
tiles x tiles maps: a 512^3 box with --tiles 4 gives 16 maps of 128 x 128.

The line of sight defaults to axis 2 (z), perpendicular to the mean field
(along x), as in the paper.  Stokes parameters use the physical dust
convention (polarization perpendicular to H) of python/analysis.py.

E and B are computed on the full periodic N x N map and then tiled: E/B are
non-local, and doing the decomposition on a non-periodic 128^2 tile would
leak E into B at the tile edges.

Output .npz (float32):  T, Q, U, E, B  each (n_maps, n, n);  plus metadata
and a per-map table of mean T, polarization fraction and BB/EE band power.
"""
import argparse, os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analysis import read_snapshot, project, qu_to_eb


def tile(a, t):
    n = a.shape[0]; m = n // t
    return np.stack([a[i*m:(i+1)*m, j*m:(j+1)*m] for i in range(t) for j in range(t)])


def band_bb_ee(Qt, Ut):
    """BB/EE of a single tile, with a cosine taper so edges don't dominate.
    A quick-look number only; use the full-map analysis for real statistics."""
    n = Qt.shape[0]
    w = np.outer(np.hanning(n), np.hanning(n))
    E, B = qu_to_eb(Qt * w, Ut * w)
    return float((B**2).sum() / max((E**2).sum(), 1e-300))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("snap")
    ap.add_argument("--out", required=True)
    ap.add_argument("--tiles", type=int, default=4)
    ap.add_argument("--axis", type=int, default=2)
    ap.add_argument("--png", default=None)
    a = ap.parse_args()

    s = read_snapshot(a.snap)
    n = s["n"]
    if n % a.tiles:
        raise SystemExit(f"N={n} not divisible by tiles={a.tiles}")
    T, Q, U = project(s["rho"], s["H"], axis=a.axis)
    E, B = qu_to_eb(Q, U)                    # on the full periodic map
    maps = {k: tile(v, a.tiles).astype(np.float32)
            for k, v in dict(T=T, Q=Q, U=U, E=E, B=B).items()}
    nm, m = maps["T"].shape[0], maps["T"].shape[1]

    rows = []
    for i in range(nm):
        p = np.sqrt(maps["Q"][i]**2 + maps["U"][i]**2) / np.maximum(maps["T"][i], 1e-30)
        rows.append(dict(map=i, tile_row=i // a.tiles, tile_col=i % a.tiles,
                         T_mean=float(maps["T"][i].mean()),
                         p_median=float(np.median(p)),
                         BB_EE_tile=band_bb_ee(maps["Q"][i], maps["U"][i])))
    meta = dict(source=os.path.abspath(a.snap), n_box=n, tiles=a.tiles, map_size=m,
                axis=a.axis, los_depth_cells=n, time=float(s["meta"].get("time", np.nan)),
                Ms=float(s["meta"].get("ms", np.nan)), Ma=float(s["meta"].get("ma", np.nan)),
                cs=float(s["meta"].get("cs", 1.0)), code=str(s["meta"].get("code", "mhd")))
    np.savez_compressed(a.out, **maps, meta=json.dumps(meta), table=json.dumps(rows))
    print(f"{nm} maps of {m}x{m} from {os.path.basename(a.snap)} "
          f"(N={n}, line of sight axis {a.axis}, full depth {n} cells) -> {a.out}")
    print(f"  M_S={meta['Ms']:.2f}  M_A={meta['Ma']:.2f}  c_s={meta['cs']:.3f}  t={meta['time']:.4g}")
    print(f"  {'map':>3} {'<T>':>10} {'p median':>9} {'BB/EE':>7}")
    for r in rows:
        print(f"  {r['map']:3d} {r['T_mean']:10.4g} {r['p_median']:9.4f} {r['BB_EE_tile']:7.3f}")
    bb = np.array([r["BB_EE_tile"] for r in rows])
    print(f"  BB/EE over tiles: median {np.median(bb):.3f}, 16-84% "
          f"[{np.percentile(bb,16):.3f}, {np.percentile(bb,84):.3f}]   (Planck ~0.53)")

    if a.png:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(3, 1, figsize=(10, 7.5))
        for ax, key, cmap in zip(axs, ["T", "Q", "U"], ["inferno", "RdBu_r", "RdBu_r"]):
            strip = np.concatenate([maps[key][i] for i in range(min(nm, 8))], axis=1)
            if key == "T":
                ax.imshow(np.log(strip).T, origin="lower", cmap=cmap)
            else:
                v = np.percentile(np.abs(strip), 99)
                ax.imshow(strip.T, origin="lower", cmap=cmap, vmin=-v, vmax=v)
            ax.set_yticks([]); ax.set_xticks([m * (i + 0.5) for i in range(min(nm, 8))])
            ax.set_xticklabels(range(min(nm, 8))); ax.set_ylabel(key)
        fig.suptitle(f"first {min(nm,8)} of {nm} maps, {m}x{m}, M_S={meta['Ms']:.2f} M_A={meta['Ma']:.2f}")
        fig.tight_layout(); fig.savefig(a.png, dpi=110); plt.close(fig)
        print(f"  quick look: {a.png}")


if __name__ == "__main__":
    main()
