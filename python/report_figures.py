"""Figures for the project report: resolution limit and simulation cost.

    python3 python/report_figures.py --out figures/report

fig_resolution.png  local slope of the 3d spectra of a 256^3 box against k:
                    flat over the cascade, steepening where the scheme's
                    dissipation takes over.  The shaded band is what a training
                    map can use.
fig_cost.png        cycles needed per dynamical time, with and without the
                    Alfven-speed limiter (Enzo 64^3, same seed) and for the
                    unlimited 512^3 box, which is why that run stalled.
"""
import argparse, glob, os, re, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analysis import local_slope
from well_stream import load_state

# categorical slots 1-3 (validated palette), text and grid tokens
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
INK, MUTED, GRID = "#1f1f1e", "#6b6a64", "#e4e3dd"


def style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
                         "figure.facecolor": "white", "axes.facecolor": "white"})
    return plt


def fig_resolution(plt, well, pair, out, kdiss=12):
    fs = sorted(glob.glob(os.path.join(well, pair, "state_*.npz")))
    rs = [load_state(f, 2) for f in fs]
    k = rs[0]["k3"]
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.axvspan(3, kdiss, color=GRID, zorder=0)
    for i, (key, lab) in enumerate([("Crho", "density"), ("Cv", "velocity"),
                                    ("CH", "magnetic field")]):
        sl = local_slope(k, np.mean([r[key] for r in rs], axis=0))
        m = (k >= 2) & (k <= 60) & np.isfinite(sl)
        ax.plot(k[m], sl[m], color=C[i], lw=2, label=lab)
        j = np.argmin(np.abs(k[m] - 34))
        ax.annotate(lab, (k[m][j], sl[m][j]), color=C[i], fontsize=9,
                    xytext=(4, -2), textcoords="offset points")
    ax.axvline(kdiss, color=MUTED, lw=1, ls="--")
    ax.text(kdiss * 1.05, ax.get_ylim()[0] + 0.25, f"end of the cascade, k = {kdiss}",
            color=MUTED, fontsize=9)
    ax.text(3.2, ax.get_ylim()[0] + 0.25, "usable by a training map", color=MUTED, fontsize=9)
    ax.set_xscale("log")
    ax.set_ylim(top=0, bottom=-8.2)          # room for the direct labels
    ax.set_xlabel("k  (wavenumber in units of the box)")
    ax.set_ylabel("local slope of the power spectrum")
    ax.set_title(f"A 256$^3$ box is turbulent over a factor {kdiss/3:.0f} in scale\n"
                 f"{pair.replace('_', ' ')}, 10 snapshots", fontsize=11, loc="left")
    ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout(); fig.savefig(out, dpi=140); plt.close(fig)
    return out


def cycles_per_tdyn(log, tdyn, tmax=5):
    pat = re.compile(r"TopGrid dt = ([0-9.eE+-]+)\s+time = ([0-9.eE+-]+)")
    t = np.array([float(m.group(2)) / tdyn for l in open(log, errors="replace")
                  for m in [pat.search(l)] if m])
    return np.array([((t >= i) & (t < i + 1)).sum() for i in range(tmax)])


def fig_cost(plt, root, out, ms=4.7):
    tdyn = 0.5 / ms
    runs = [("64$^3$, no limiter", "nofloor", C[1]), ("64$^3$, Alfven cap 50", "vA50", C[2])]
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    for lab, d, col in runs:
        y = cycles_per_tdyn(os.path.join(root, d, "enzo.log"), tdyn)
        x = np.arange(len(y)) + 0.5
        ax.plot(x, y, color=col, lw=2, marker="o", ms=5, label=lab)
        ax.annotate(lab, (x[-1], y[-1]), color=col, fontsize=9, xytext=(6, -3),
                    textcoords="offset points")
    # the 512^3 production box, measured from its own log (peek report)
    y512 = [3458, 47672, 64051]
    ax.plot(np.arange(len(y512)) + 0.5, y512, color=C[0], lw=2, marker="s", ms=5,
            label="512$^3$, no limiter")
    ax.annotate("512$^3$, no limiter\n(stalled here)", (2.5, y512[-1]), color=C[0], fontsize=9,
                xytext=(6, -14), textcoords="offset points")
    ax.set_yscale("log")
    ax.set_xlabel("dynamical times completed")
    ax.set_ylabel("time steps needed per dynamical time")
    ax.set_title("Without a cap on the Alfven speed, the cost per dynamical time keeps growing",
                 fontsize=11, loc="left")
    ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True)
    ax.set_xlim(0, 6.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout(); fig.savefig(out, dpi=140); plt.close(fig)
    return out


def fig_dust_summary(plt, summary_json, out):
    """Compact print version of the comparison: median and 16-84% range for the
    four statistics that carry the conclusion, dust models against simulations."""
    import json
    d = json.load(open(summary_json))
    names = [n for n in d if n.startswith("PySM")] + \
            [n for n in d if not n.startswith("PySM") and "face" in n]
    short = {n: n.replace("PySM ", "PySM3 ").replace(" face", "")
              .replace("enzo ", "Enzo ").replace("face(128/512)", "").strip() for n in names}
    panels = [("aEE", r"$\alpha_{EE}$", -2.42), ("bbee", "BB/EE", 0.53),
              ("rte", r"$r_{TE}$", 0.36), ("S_med", r"$S$ at 2 px [deg]", None)]
    fig, axs = plt.subplots(1, 4, figsize=(11, 3.7), sharey=True)
    y = np.arange(len(names))[::-1]
    for ax, (key, lab, planck) in zip(axs, panels):
        for yy, n in zip(y, names):
            lo, md, hi = d[n][key]
            col = C[0] if n.startswith("PySM") else (C[2] if n.startswith("enzo") else C[3])
            ax.plot([lo, hi], [yy, yy], color=col, lw=2.4, solid_capstyle="round", alpha=.55)
            ax.plot([md], [yy], "o", color=col, ms=6.5, zorder=3)
        if planck is not None:
            ax.axvline(planck, color=MUTED, lw=1, ls="--")
            ax.annotate("Planck", xy=(planck, 1), xycoords=("data", "axes fraction"),
                        xytext=(3, -11), textcoords="offset points",
                        color=MUTED, fontsize=8)
        ax.set_title(lab, fontsize=11, color=INK, loc="left")
        ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True)
        for sp in ("top", "right", "left"):
            ax.spines[sp].set_visible(False)
    axs[0].set_yticks(y); axs[0].set_yticklabels([short[n] for n in names], fontsize=9, color=INK)
    fig.suptitle("Dust models (blue) against simulated maps: median and 16-84% range "
                 "over maps", fontsize=10, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out, dpi=150); plt.close(fig)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--well", default="runs_well")
    ap.add_argument("--pair", default="Ma_0.7_Ms_7")
    ap.add_argument("--floor", default="runs_enzo/floor_test")
    ap.add_argument("--kdiss", type=float, default=12)
    ap.add_argument("--out", default="figures/report")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    plt = style()
    print(fig_resolution(plt, a.well, a.pair, os.path.join(a.out, "fig_resolution.png"), a.kdiss))
    print(fig_cost(plt, a.floor, os.path.join(a.out, "fig_cost.png")))
    sj = "figures/well_vs_pysm/summary.json"
    if os.path.exists(sj):
        print(fig_dust_summary(plt, sj, os.path.join(a.out, "fig_dust_summary.png")))


if __name__ == "__main__":
    main()
