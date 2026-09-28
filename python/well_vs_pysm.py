"""Are The Well MHD dust maps close to PySM3 dust?  Same estimator on both.

    python3 python/well_vs_pysm.py --well runs_well --pysm pysm_patches \
        --pairs Ma_2_Ms_7 Ma_0.7_Ms_7 --out figures/well_vs_pysm

PySM: 128^2 patches of 20 deg at 353 GHz, 35 <= |b| <= 75, smoothed with a
2-pixel FWHM Gaussian, E/B decomposed on the full sphere (pysm_patches.py).

MHD (The Well state files from well_stream.py), line of sight perpendicular to
the mean field (both such axes), two ways of making a 128^2 map:
  face  the whole 256^2 projected face, 2x2-binned: the patch spans the box,
        so patch k = box k and the fit band k=[3,13] is the box's inertial range
  tile  the 256^2 face cut into 4 native-resolution 128^2 tiles.  A half-width
        tile resolves a given patch k at TWICE the box wavenumber, so patch
        k 3-13 is box k 6-26: the upper half of that lies beyond the box's
        dissipation knee (k ~ 12), which is why tiles measure steeper slopes
In both, the map is smoothed with the same 2-pixel beam, E/B are decomposed on
the full periodic face, and the patch is then tapered exactly like PySM.

All statistics (python/patch_stats.py) over patch k in [KMIN, KMAX], beam
divided out.  For a 20 deg patch, k = 3..13 is ell ~ 54..234.
"""
import argparse, glob, json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analysis import qu_to_eb
from patch_stats import patch_stats, KEYS
from well_stream import state_map

PLANCK = dict(aEE=-2.42, aBB=-2.54, bbee=0.53, rte=0.36)     # Planck 2018 XI, ell 40-600
# colour follows the source, never the position; the native-resolution "tile"
# product is the same source as "face" and is told apart by hatching
SOURCE_COLOR = {"PySM d1": "#2a78d6", "PySM d10": "#eb6834", "PySM d12": "#1baf7a",
                "Ma2": "#eda100", "Ma0.7": "#e87ba4", "enzo": "#008300"}


def color_of(name):
    if name.startswith("PySM"):
        return SOURCE_COLOR.get(name, "#8a8980")
    if name.startswith("enzo") or name.startswith("sim"):
        return SOURCE_COLOR["enzo"]
    return SOURCE_COLOR["Ma2" if " Ma2 " in name else "Ma0.7"]


def smooth_periodic(a, fwhm_px):
    n = a.shape[0]
    k = np.fft.fftfreq(n)
    kk2 = k[:, None]**2 + k[None, :]**2
    sig = fwhm_px / 2.3548
    return np.real(np.fft.ifft2(np.fft.fft2(a) * np.exp(-2 * np.pi**2 * sig**2 * kk2)))


def bin2(a):
    n = a.shape[0] // 2
    return a.reshape(n, 2, n, 2).mean(axis=(1, 3))


def mhd_stats(pair_dir, product, kmin, kmax, fwhm_px=2.0):
    out = []
    for f in sorted(glob.glob(os.path.join(pair_dir, "state_*.npz"))):
        z = np.load(f, allow_pickle=True)
        meta = json.loads(str(z["meta"]))
        for ax in (1, 2):
            T, Q, U = (state_map(z, meta, ax, m) for m in "TQU")
            if product == "face":
                T, Q, U = bin2(T), bin2(Q), bin2(U)
            T, Q, U = (smooth_periodic(x, fwhm_px) for x in (T, Q, U))
            E, B = qu_to_eb(Q, U)                       # periodic face: no leakage
            n = T.shape[0]
            bsig = (fwhm_px / 2.3548) / 128.0            # in units of a 128-px patch
            if product == "face":
                cuts = [(slice(0, n), slice(0, n))]
            else:
                h = n // 2
                cuts = [(slice(i, i + h), slice(j, j + h)) for i in (0, h) for j in (0, h)]
            for sx, sy in cuts:
                out.append(patch_stats(T[sx, sy], Q[sx, sy], U[sx, sy], kmin=kmin, kmax=kmax,
                                       E=E[sx, sy], B=B[sx, sy], beam_sigma=bsig))
    return out


def pysm_stats(path, kmin, kmax):
    z = np.load(path)
    bs = float(z["beam_sigma"])
    return [patch_stats(z["T"][i], z["Q"][i], z["U"][i], kmin=kmin, kmax=kmax,
                        E=z["E"][i], B=z["B"][i], beam_sigma=bs) for i in range(len(z["T"]))], z


def untile(a, t):
    """Reassemble the t x t tiles written by extract_maps.py into the full face
    (they are stored row-major: index i*t + j covers rows i, columns j)."""
    m = a.shape[-1]
    return np.block([[a[i * t + j] for j in range(t)] for i in range(t)])


def enzo_series(path, kmin, kmax, fwhm_px=2.0, product="face"):
    """Maps written by extract_maps.py, with E/B already decomposed on the full
    periodic face before tiling.

    product 'face': the tiles are reassembled into the whole face and binned to
    128^2, so patch k = box k and the fit band sits inside the cascade -- the
    like-for-like counterpart of the Well 'face' product.
    product 'tile': the native-resolution tiles as written; their patch k is
    box k / (map size / box size), i.e. four times higher for 128 of 512, which
    is past the dissipation knee.  Kept only to show that difference.
    """
    z = np.load(path, allow_pickle=True)
    meta = json.loads(str(z["meta"]))
    t = int(round(np.sqrt(len(z["T"]))))
    st, maps = [], {k: z[k] for k in ("T", "Q", "U", "E", "B")}
    if product == "face":
        full = {k: untile(v, t) for k, v in maps.items()}
        while full["T"].shape[0] > 128:
            full = {k: bin2(v.astype(np.float64)) for k, v in full.items()}
        items = [tuple(full[k] for k in ("T", "Q", "U", "E", "B"))]
    else:
        items = [tuple(maps[k][i].astype(np.float64) for k in ("T", "Q", "U", "E", "B"))
                 for i in range(len(maps["T"]))]
    for T, Q, U, E, B in items:
        bsig = (fwhm_px / 2.3548) / T.shape[-1]
        T, Q, U, E, B = (smooth_periodic(x, fwhm_px) for x in (T, Q, U, E, B))
        st.append(patch_stats(T, Q, U, kmin=kmin, kmax=kmax, E=E, B=B, beam_sigma=bsig))
    name = (f"{meta.get('code', 'sim')} Ms{meta['Ms']:.1f} Ma{meta['Ma']:.1f} "
            f"{product}({items[0][0].shape[0]}/{meta['n_box']})")
    return name, st


def summary(stats):
    return {k: np.percentile([s[k] for s in stats], [16, 50, 84]) for k in KEYS}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--well", default="runs_well")
    ap.add_argument("--pysm", default="pysm_patches")
    ap.add_argument("--models", nargs="+", default=["d1", "d10", "d12"])
    ap.add_argument("--pairs", nargs="+", default=["Ma_2_Ms_7", "Ma_0.7_Ms_7"])
    ap.add_argument("--enzo", nargs="*", default=[],
                    help="maps_*.npz from python/extract_maps.py (Enzo or our solver)")
    ap.add_argument("--kmin", type=float, default=3)
    ap.add_argument("--kmax", type=float, default=13)
    ap.add_argument("--out", default="figures/well_vs_pysm")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    series, refs = [], {}
    for m in a.models:
        p = os.path.join(a.pysm, f"{m}.npz")
        if os.path.exists(p):
            st, z = pysm_stats(p, a.kmin, a.kmax)
            series.append((f"PySM {m}", st))
            refs[m] = {k: float(z["ref_" + k]) for k in ("aTT", "aEE", "aBB", "bbee", "rte")}
    for pair in a.pairs:
        d = os.path.join(a.well, pair)
        if not glob.glob(os.path.join(d, "state_*.npz")):
            continue
        ms = pair.split("_Ms_")[1]; ma = pair.split("_")[1]
        for prod in ("face", "tile"):
            series.append((f"Well Ms{ms} Ma{ma} {prod}", mhd_stats(d, prod, a.kmin, a.kmax)))

    for f in a.enzo:
        if os.path.exists(f):
            for prod in ("face", "tile"):
                series.append(enzo_series(f, a.kmin, a.kmax, product=prod))

    # ---- table
    lines = [f"Statistics over patch k = [{a.kmin:g}, {a.kmax:g}] (20 deg patch: ell ~ "
             f"{360/20*a.kmin:.0f}-{360/20*a.kmax:.0f}); median [16%, 84%] over maps",
             f"Planck 2018 XI (ell 40-600): {PLANCK}",
             "p_med: PySM = P/I; MHD = (P/T)/p0, i.e. still to be multiplied by the "
             "intrinsic polarization fraction p0", ""]
    for m, r in refs.items():
        lines.append(f"PySM {m} full sky |b|>35, ell 40-600: " +
                     " ".join(f"{k}={v:.2f}" for k, v in r.items()))
    lines.append("")
    # one block per statistic keeps the table readable in a terminal
    for k in KEYS:
        lines.append(f"{k}:")
        for name, st in series:
            q = summary(st)[k]
            lines.append(f"    {name:<26} N={len(st):4d}   {q[1]:8.3f}   [{q[0]:8.3f}, {q[2]:8.3f}]")
    txt = "\n".join(lines)
    print(txt)
    open(os.path.join(a.out, "table.txt"), "w").write(txt + "\n")
    json.dump({name: {k: [float(x) for x in v] for k, v in summary(st).items()}
               for name, st in series}, open(os.path.join(a.out, "summary.json"), "w"), indent=1)

    # ---- figure: one small multiple per statistic, one box per series
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    show = [("aEE", r"$\alpha_{EE}$"), ("aBB", r"$\alpha_{BB}$"), ("aTT", r"$\alpha_{TT}$"),
            ("bbee", "BB/EE"), ("rte", r"$r_{TE}$"), ("S_med", r"median $S$ at 2 px [deg]"),
            ("S_p_slope", r"d ln$S$ / d ln$p$"), ("skew_lnT", "skewness of ln T"),
            ("kurt_Q", "excess kurtosis of Q")]
    ink, muted, grid = "#1f1f1e", "#6b6a64", "#e4e3dd"
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": muted, "axes.labelcolor": ink,
                         "xtick.color": muted, "ytick.color": muted, "text.color": ink})
    fig, axs = plt.subplots(3, 3, figsize=(12, 10.5))
    names = [n for n, _ in series]
    for ax, (key, lab) in zip(axs.flat, show):
        data = [[s[key] for s in st] for _, st in series]
        bp = ax.boxplot(data, vert=True, patch_artist=True, widths=0.55, showfliers=False,
                        medianprops=dict(color=ink, lw=1.5), whiskerprops=dict(color=muted),
                        capprops=dict(color=muted), boxprops=dict(lw=0))
        for name, b in zip(names, bp["boxes"]):
            b.set_facecolor(color_of(name))
            if name.endswith("tile"):
                b.set_hatch("////"); b.set_edgecolor("#fcfcfb")
        if key in PLANCK:
            ax.axhline(PLANCK[key], color=muted, lw=1, ls="--")
            ax.text(len(series) + 0.45, PLANCK[key], "Planck", color=muted, va="bottom",
                    ha="right", fontsize=8)
        ax.set_title(lab, fontsize=10, color=ink, loc="left")
        ax.set_xticks(range(1, len(series) + 1))
        ax.set_xticklabels(names, rotation=40, ha="right", fontsize=8, color=ink)
        ax.grid(axis="y", color=grid, lw=0.6); ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    fig.suptitle(f"Same estimator on PySM3 dust patches and The Well MHD maps "
                 f"(patch k = {a.kmin:g}-{a.kmax:g}; boxes: 25-75%, whiskers: 1.5 IQR)",
                 fontsize=11, color=ink)
    fig.tight_layout()
    fig.savefig(os.path.join(a.out, "stats.png"), dpi=130)
    plt.close(fig)

    # ---- example maps: one PySM d10 patch and one MHD face, same layout
    ex = []
    zp = os.path.join(a.pysm, "d10.npz")
    if os.path.exists(zp):
        z = np.load(zp)
        i = int(np.argmin(np.abs(z["centers"][:, 1] - 60)))
        ex.append((f"PySM d10 (l={z['centers'][i,0]:.0f}, b={z['centers'][i,1]:.0f})",
                   z["T"][i], z["Q"][i], z["U"][i]))
    for pair in a.pairs:
        fs = sorted(glob.glob(os.path.join(a.well, pair, "state_*.npz")))
        if fs:
            z = np.load(fs[-1])
            T, Q, U = (smooth_periodic(bin2(z[f"ax2_map{m}"].astype(float)), 2) for m in "TQU")
            ex.append((f"Well {pair} (face)", T, Q, U))
    if ex:
        fig, axs = plt.subplots(len(ex), 3, figsize=(9, 3 * len(ex)), squeeze=False)
        for r, (title, T, Q, U) in enumerate(ex):
            for c, (m, cm) in enumerate([(T, "inferno"), (Q, "RdBu_r"), (U, "RdBu_r")]):
                ax = axs[r, c]
                if c == 0:
                    ax.imshow(np.log(np.maximum(m, m.max() * 1e-4)).T, origin="lower", cmap=cm)
                else:
                    s = np.percentile(np.abs(m), 99)
                    ax.imshow(m.T, origin="lower", cmap=cm, vmin=-s, vmax=s)
                ax.set_xticks([]); ax.set_yticks([])
                ax.set_title(f"{title}: {['ln T', 'Q', 'U'][c]}", fontsize=8, color=ink)
        fig.tight_layout()
        fig.savefig(os.path.join(a.out, "maps.png"), dpi=120)
        plt.close(fig)
    print(f"\nwrote {a.out}/table.txt, stats.png, maps.png")


if __name__ == "__main__":
    main()
