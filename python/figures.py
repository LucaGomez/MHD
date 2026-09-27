"""Reproduce the figures and Table 1 of Stalpes, Collins & Huffenberger (2024)."""
import os, sys, glob, json, argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, SymLogNorm
from matplotlib.cm import ScalarMappable
import paper_values as PV

# Planck 2018 XI values quoted in the paper
PL = dict(aEE=-2.42, aBB=-2.54, BB_EE=0.53, rTE=0.355, rTB=0.055)

CMAP = plt.get_cmap("viridis")
NORM = Normalize(vmin=0, vmax=7)

def color(ms):  return CMAP(NORM(ms))
def size(ma):   return 25 + 55 * (ma / 2.0)
def lstyle(ma): return ":" if ma < 0.7 else ("--" if ma < 1.2 else "-")

def load(root, axis=2):
    """Load one projection axis only.  Globbing analysis_ax*.npz would mix the
    perpendicular (axis 2) and parallel (axis 0) projections into one sample."""
    out = []
    pat = os.path.join(root, "*", f"analysis_ax{axis}.npz")
    for f in sorted(glob.glob(pat)):
        d = dict(np.load(f, allow_pickle=True))
        d["tag"] = os.path.basename(os.path.dirname(f))
        out.append(d)
    out.sort(key=lambda d: (float(d["Ms"]), float(d["Ma"])))
    if not out:
        raise SystemExit(f"no files matching {pat}")
    return out


def sc(ax, x, y, runs, yerr=None):
    for i, r in enumerate(runs):
        ax.scatter(x[i], y[i], color=color(float(r["Ms"])), s=size(float(r["Ma"])),
                   edgecolor="k", linewidth=0.3, zorder=3)
    if yerr is not None:
        ax.errorbar(x, y, yerr=yerr, fmt="none", ecolor="0.6", lw=0.8, zorder=2)


def fig_spectra(runs, keys, kk, comp, labels, fn, ylab):
    fig, axes = plt.subplots(1, len(keys), figsize=(4.1*len(keys), 3.6))
    for ax, key, lab in zip(np.atleast_1d(axes), keys, labels):
        for r in runs:
            k, C = r[kk], r[key]
            ax.loglog(k, C * k**comp, color=color(float(r["Ms"])),
                      ls=lstyle(float(r["Ma"])), lw=1.1)
        k1, k2 = r["kfit"]
        ax.axvline(k1, color="0.8", lw=0.8); ax.axvline(k2, color="0.8", lw=0.8)
        ax.set_xlabel("$k/k_{min}$"); ax.set_title(lab)
        ax.set_ylabel(ylab % comp)
    fig.colorbar(ScalarMappable(norm=NORM, cmap=CMAP), ax=np.atleast_1d(axes).tolist(),
                 label="$M_S$")
    fig.savefig(fn, dpi=140, bbox_inches="tight"); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="runs/n64")
    ap.add_argument("--outdir", default="figures/n64")
    ap.add_argument("--label", default="")
    ap.add_argument("--axis", type=int, default=2,
                    help="2 = perpendicular to the mean field (the paper's "
                         "default), 0 = along it (their section 3.6)")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    runs = load(a.root, axis=a.axis)
    if not runs:
        raise SystemExit("no analysis files in " + a.root)
    MS = np.array([float(r["Ms"]) for r in runs])
    MA = np.array([float(r["Ma"]) for r in runs])
    print(f"{len(runs)} runs (axis={a.axis}, "
          f"{'perpendicular to' if a.axis == 2 else 'along'} the mean field), "
          f"Ms in [{MS.min():.2f},{MS.max():.2f}], "
          f"Ma in [{MA.min():.2f},{MA.max():.2f}]")

    # ---- Fig 2: legend of achieved Mach numbers
    fig, ax = plt.subplots(figsize=(4.2, 3.4))
    sc(ax, MA, MS, runs)
    ax.set_xlabel("$M_A$"); ax.set_ylabel("$M_S$"); ax.set_title("achieved Mach numbers")
    fig.savefig(f"{a.outdir}/fig02_mach_legend.png", dpi=140, bbox_inches="tight"); plt.close(fig)

    # ---- Fig 4 / 6: spectra
    fig_spectra(runs, ["Crho", "Cv", "CH"], "k3", 11/3.,
                [r"$\rho$", r"$v$", r"$H$"], f"{a.outdir}/fig04_fluid_spectra.png",
                r"$k^{%.2f}C_k$")
    fig_spectra(runs, ["CTT", "CEE", "CBB"], "k2", 2.5,
                [r"$T$", r"$E$", r"$B$"], f"{a.outdir}/fig06_projected_spectra.png",
                r"$k^{%.1f}C_k$")

    # ---- Fig 5: fluid slopes
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.4))
    for ax, key, lab, xv, xlab in [
            (axs[0], "al_Crho", r"$\alpha_{\rho\rho}$", MS, "$M_S$"),
            (axs[1], "al_Cv",  r"$\alpha_{vv}$",       MS, "$M_S$"),
            (axs[2], "al_CH",  r"$\alpha_{HH}$",       MA, "$M_A$")]:
        y = np.array([float(r[key]) for r in runs])
        e = np.array([float(r[key+"_std"]) for r in runs])
        sc(ax, xv, y, runs, e); ax.set_xlabel(xlab); ax.set_ylabel(lab)
    axs[1].axhline(-11/3, color="0.7", ls="--", lw=0.8)
    axs[1].axhline(-2.5, color="0.7", ls=":", lw=0.8)
    fig.tight_layout(); fig.savefig(f"{a.outdir}/fig05_fluid_slopes.png", dpi=140); plt.close(fig)

    # ---- Fig 7: projected slopes
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.4))
    for ax, key, lab, ref in [(axs[0], "al_CTT", r"$\alpha_{TT}$", None),
                              (axs[1], "al_CEE", r"$\alpha_{EE}$", PL["aEE"]),
                              (axs[2], "al_CBB", r"$\alpha_{BB}$", PL["aBB"])]:
        y = np.array([float(r[key]) for r in runs])
        e = np.array([float(r[key+"_std"]) for r in runs])
        sc(ax, MS, y, runs, e); ax.set_xlabel("$M_S$"); ax.set_ylabel(lab)
        if ref: ax.axhline(ref, color="0.6", lw=1)
    fig.tight_layout(); fig.savefig(f"{a.outdir}/fig07_projected_slopes.png", dpi=140); plt.close(fig)

    # ---- Fig 3: alpha_EE vs A_BB/A_EE
    rat = np.array([float(r["A_CBB"])/float(r["A_CEE"]) for r in runs])
    aEE = np.array([float(r["al_CEE"]) for r in runs])
    fig, ax = plt.subplots(figsize=(5, 4))
    sc(ax, rat, aEE, runs)
    ax.axhline(PL["aEE"], color="0.6"); ax.axvline(PL["BB_EE"], color="0.6")
    ax.set_xlabel(r"$A_{BB}/A_{EE}$"); ax.set_ylabel(r"$\alpha_{EE}$")
    fig.colorbar(ScalarMappable(norm=NORM, cmap=CMAP), ax=ax, label="$M_S$")
    fig.savefig(f"{a.outdir}/fig03_summary.png", dpi=140, bbox_inches="tight"); plt.close(fig)

    # ---- Fig 8: amplitude ratios
    fig, axs = plt.subplots(3, 1, figsize=(5, 7), sharex=True)
    for ax, num, den, lab in [(axs[0], "A_CEE", "A_CTT", r"$A_{EE}/A_{TT}$"),
                              (axs[1], "A_CBB", "A_CTT", r"$A_{BB}/A_{TT}$"),
                              (axs[2], "A_CBB", "A_CEE", r"$A_{BB}/A_{EE}$")]:
        y = np.array([float(r[num])/float(r[den]) for r in runs])
        sc(ax, MS, y, runs); ax.set_yscale("log"); ax.set_ylabel(lab)
        ax.axhline(0.5, color="0.7")
    axs[2].set_xlabel("$M_S$")
    fig.tight_layout(); fig.savefig(f"{a.outdir}/fig08_amplitude_ratios.png", dpi=140); plt.close(fig)

    # ---- Fig 9: correlations
    fig, axs = plt.subplots(2, 3, figsize=(11, 6))
    for j, (xy, ref, yl) in enumerate([("TE", PL["rTE"], 1.0), ("TB", PL["rTB"], 0.25),
                                       ("EB", None, 0.25)]):
        for r in runs:
            axs[0, j].semilogx(r["k2"], r["rk"+xy], color=color(float(r["Ms"])),
                               ls=lstyle(float(r["Ma"])), lw=1.0)
        k1, k2 = runs[0]["kfit"]
        axs[0, j].axvline(k1, color="0.85"); axs[0, j].axvline(k2, color="0.85")
        axs[0, j].set_ylim(-yl, yl); axs[0, j].set_xlabel("$k/k_{min}$")
        axs[0, j].set_ylabel(f"$r_k^{{{xy}}}$")
        y = np.array([float(r["r"+xy]) for r in runs])
        e = np.array([float(r["r"+xy+"_std"]) for r in runs])
        sc(axs[1, j], MS, y, runs, e)
        axs[1, j].set_xlabel("$M_S$"); axs[1, j].set_ylabel(f"$\\langle r^{{{xy}}}\\rangle$")
        if ref is not None: axs[1, j].axhline(ref, color="0.6")
        axs[1, j].axhline(0, color="0.85", lw=0.6)
    fig.tight_layout(); fig.savefig(f"{a.outdir}/fig09_correlations.png", dpi=140); plt.close(fig)

    # ---- Fig 1: maps
    want = [(0.5, 0.5), (0.5, 1), (0.5, 2), (3, 0.5), (3, 1), (3, 2), (6, 0.5), (6, 1), (6, 2)]
    def nearest(ms, ma):
        d = (MS/ms - 1)**2 + (MA/ma - 1)**2
        return runs[int(np.argmin(d))]
    for mp, cmapn in [("T", "inferno"), ("E", "RdBu_r"), ("B", "RdBu_r")]:
        fig, axs = plt.subplots(3, 3, figsize=(8.4, 8.4))
        for i, (ms, ma) in enumerate(want):
            r = nearest(ms, ma); ax = axs[i//3, i%3]
            m = r["map"+mp]
            if mp == "T":
                ax.imshow(np.log(m).T, origin="lower", cmap=cmapn)
            else:
                s = 3*np.std(m)
                ax.imshow(m.T, origin="lower", cmap=cmapn,
                          norm=SymLogNorm(linthresh=0.1*s, vmin=-s, vmax=s))
            ax.set_xticks([]); ax.set_yticks([])
            ax.set_title(f"$M_S$={float(r['Ms']):.1f}, $M_A$={float(r['Ma']):.1f}", fontsize=9)
        fig.suptitle({"T": r"$\ln T$", "E": "E", "B": "B"}[mp])
        fig.tight_layout(); fig.savefig(f"{a.outdir}/fig01_maps_{mp}.png", dpi=140); plt.close(fig)

    # ---- Table 1: linear fits  q = a + b MS + c MA
    A = np.column_stack([np.ones_like(MS), MS, MA])
    rows, tab = [], {}
    for key, nm in [("al_Crho", "alpha_rho"), ("al_Cv", "alpha_v"), ("al_CH", "alpha_H"),
                    ("al_CTT", "alpha_TT"), ("al_CEE", "alpha_EE"), ("al_CBB", "alpha_BB")]:
        y = np.array([float(r[key]) for r in runs])
        c, *_ = np.linalg.lstsq(A, y, rcond=None)
        res = y - A @ c
        rows.append((nm, c[0], c[1], c[2], c[2]/c[1] if c[1] else np.nan, res.std()))
        tab[nm] = dict(a=c[0], b=c[1], c=c[2], scatter=res.std())
    txt = [f"Linear fits  alpha_q = a + b M_S + c M_A     ({a.label or a.root}, "
           f"{len(runs)} runs, N={int(runs[0]['n'])}, axis={a.axis}, "
           f"fit window k=[{runs[0]['kfit'][0]},"
           f"{runs[0]['kfit'][1]}])",
           "",
           "                    this work                 |        paper (Table 1)",
           "  quantity        a       b       c    scatter |      a       b       c",
           "  " + "-"*72]
    for nm, a0, b0, c0, cb, rs in rows:
        pv = PV.TABLE1.get(nm)
        pstr = (f"  {pv['a']:7.2f} {pv['b']:7.3f} {pv['c']:7.3f}" if pv else "")
        txt.append(f"  {nm:<12s} {a0:7.2f} {b0:7.3f} {c0:7.3f} {rs:8.3f}  |{pstr}")
    # inferred ISM Mach numbers from alpha_EE = -2.4, alpha_BB = -2.5
    M = np.array([[tab["alpha_EE"]["b"], tab["alpha_EE"]["c"]],
                  [tab["alpha_BB"]["b"], tab["alpha_BB"]["c"]]])
    rhs = np.array([PL["aEE"] - tab["alpha_EE"]["a"], PL["aBB"] - tab["alpha_BB"]["a"]])
    try:
        sol = np.linalg.solve(M, rhs)
        txt.append(f"\n  inferred from alpha_EE={PL['aEE']}, alpha_BB={PL['aBB']}: "
                   f"M_S = {sol[0]:.2f}, M_A = {sol[1]:.2f}   "
                   f"(paper: {PV.INFERRED['Ms']}, {PV.INFERRED['Ma']})")
        tab["inferred"] = dict(Ms=float(sol[0]), Ma=float(sol[1]))
    except np.linalg.LinAlgError:
        pass
    hi = MS > 4
    if hi.sum():
        txt.append(f"  <A_BB/A_EE> for M_S>4 : {rat[hi].mean():.3f} +- {rat[hi].std():.3f} "
                   f"(paper {PV.HIGH_MS['A_BB/A_EE'][0]} +- {PV.HIGH_MS['A_BB/A_EE'][1]}, "
                   f"Planck {PL['BB_EE']})")
        ret = np.array([float(r["A_CEE"])/float(r["A_CTT"]) for r in runs])
        rbt = np.array([float(r["A_CBB"])/float(r["A_CTT"]) for r in runs])
        txt.append(f"  <A_EE/A_TT> for M_S>4 : {ret[hi].mean():.3f} +- {ret[hi].std():.3f} "
                   f"(paper {PV.HIGH_MS['A_EE/A_TT']})")
        txt.append(f"  <A_BB/A_TT> for M_S>4 : {rbt[hi].mean():.3f} +- {rbt[hi].std():.3f} "
                   f"(paper {PV.HIGH_MS['A_BB/A_TT']})")
        tab["BB_EE_highMs"] = [float(rat[hi].mean()), float(rat[hi].std())]
        rte = np.array([float(r["rTE"]) for r in runs])
        txt.append(f"  <r_TE>     for M_S>4 : {rte[hi].mean():.3f} +- {rte[hi].std():.3f} "
                   f"(Planck 0.355)")
        tab["rTE_highMs"] = [float(rte[hi].mean()), float(rte[hi].std())]
        for xy in ["TB", "EB"]:
            v = np.array([float(r["r"+xy]) for r in runs])
            txt.append(f"  r_{xy} > 0 : {int((v>0).sum())}/{len(v)} runs "
                       f"({int((v[hi]>0).sum())}/{int(hi.sum())} for M_S>4);  "
                       f"paper {PV.PARITY['r'+xy+'_positive'][0]}/21 "
                       f"({PV.PARITY['r'+xy+'_positive_highMs'][0]}/9);  "
                       f"mean |r_{xy}| = {np.abs(v).mean():.3f}")
            tab["r"+xy+"_pos"] = [int((v>0).sum()), len(v)]
    out = "\n".join(txt)
    print(out)
    open(f"{a.outdir}/table1.txt", "w").write(out + "\n")
    json.dump(tab, open(f"{a.outdir}/summary.json", "w"), indent=1)

    # ---- per-run table
    with open(f"{a.outdir}/runs_table.txt", "w") as f:
        f.write(f"{'tag':<12}{'Ms':>6}{'Ma':>6}{'a_rho':>8}{'a_v':>8}{'a_H':>8}"
                f"{'a_TT':>8}{'a_EE':>8}{'a_BB':>8}{'BB/EE':>8}{'rTE':>7}{'rTB':>7}{'rEB':>7}\n")
        for r in runs:
            f.write(f"{str(r['tag']):<12}{float(r['Ms']):6.2f}{float(r['Ma']):6.2f}"
                    f"{float(r['al_Crho']):8.2f}{float(r['al_Cv']):8.2f}{float(r['al_CH']):8.2f}"
                    f"{float(r['al_CTT']):8.2f}{float(r['al_CEE']):8.2f}{float(r['al_CBB']):8.2f}"
                    f"{float(r['A_CBB'])/float(r['A_CEE']):8.2f}"
                    f"{float(r['rTE']):7.2f}{float(r['rTB']):7.3f}{float(r['rEB']):7.3f}\n")
    print(open(f"{a.outdir}/runs_table.txt").read())

if __name__ == "__main__":
    main()
