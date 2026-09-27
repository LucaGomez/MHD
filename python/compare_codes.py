"""Side-by-side comparison of Enzo and our solver (src/mhd.c) on the same
physical point, analysed by the identical pipeline (python/process_run.py).

    python python/compare_codes.py runs_xcheck/mhd/ms4.7_ma1.5_n32 \
                                   runs_xcheck/enzo/ms4.7_ma1.5_n32 --out figures/xcheck_n32

Prints a table of every quantity the paper reports -- Mach numbers, spectral
slopes of rho, v, H, T, E, B, band-power ratios and T/E/B correlations -- for
both codes and both projections, with the difference in units of the combined
snapshot-to-snapshot scatter, next to the paper and Planck values.  Also
overlays the spectra.
"""
import argparse, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paper_values as PV

ROWS = [("Ms", "M_S", None), ("Ma", "M_A", None),
        ("al_Crho", "alpha_rho", "alpha_rho"), ("al_Cv", "alpha_v", "alpha_v"),
        ("al_CH", "alpha_H", "alpha_H"), ("al_CTT", "alpha_TT", "alpha_TT"),
        ("al_CEE", "alpha_EE", "alpha_EE"), ("al_CBB", "alpha_BB", "alpha_BB"),
        ("bp_BB_EE", "BB/EE band", None), ("bp_EE_TT", "EE/TT band", None),
        ("rTE", "r_TE", None), ("rTB", "r_TB", None), ("rEB", "r_EB", None)]


def paper_expect(key, ms, ma):
    """Value the paper's Table-1 linear fit predicts at the achieved (M_S, M_A)."""
    if key in PV.TABLE1:
        c = PV.TABLE1[key]
        return c["a"] + c["b"] * ms + c["c"] * ma
    return None


def load(d, axis):
    f = os.path.join(d, f"analysis_ax{axis}.npz")
    if not os.path.exists(f):
        raise SystemExit(f"missing {f}; run python/process_run.py {d} --axis {axis}")
    return dict(np.load(f, allow_pickle=True))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mhd"); ap.add_argument("enzo")
    ap.add_argument("--out", default="figures/xcheck")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    lines = []
    for axis, label in [(2, "perpendicular to B0 (the paper's case)"),
                        (0, "along B0")]:
        m, e = load(a.mhd, axis), load(a.enzo, axis)
        n = int(m["n"])
        lines += ["", f"=== projection {label},  N={n},  fit window "
                      f"k=[{m['kfit'][0]},{m['kfit'][1]}],  snapshots: mhd {int(m['nsnap'])}, "
                      f"enzo {int(e['nsnap'])}",
                  f"  {'quantity':<12}{'our solver':>14}{'Enzo':>14}{'diff/sigma':>11}"
                  f"{'paper fit':>11}   (paper Table-1 fit evaluated at Enzo's M_S, M_A)"]
        for key, name, pkey in ROWS:
            if key not in m or key not in e:
                continue
            vm, ve = float(m[key]), float(e[key])
            sm = float(m.get(key + "_std", np.nan)); se = float(e.get(key + "_std", np.nan))
            sig = np.sqrt(np.nansum([sm**2, se**2]))
            z = (vm - ve) / sig if sig > 0 else np.nan
            pe = paper_expect(pkey, float(e["Ms"]), float(e["Ma"])) if pkey else None
            pes = f"{pe:11.2f}" if pe is not None else " " * 11
            fs = lambda v, sd: (f"{v:8.3f}±{sd:<5.2f}" if np.isfinite(sd) else f"{v:8.3f}      ")
            zs = f"{z:11.1f}" if np.isfinite(z) else " " * 11
            lines.append(f"  {name:<12}{fs(vm, sm)}{fs(ve, se)}{zs}{pes}")
        if axis == 2:
            lines.append(f"  Planck: alpha_EE {PV.PLANCK['aEE']}, alpha_BB {PV.PLANCK['aBB']}, "
                         f"BB/EE {PV.PLANCK['BB_EE']}, r_TE {PV.PLANCK['rTE']}")
        # spectra overlay
        fig, axs = plt.subplots(2, 3, figsize=(12, 6.5))
        for ax, key, comp, kk in zip(axs.flat,
                                     ["Crho", "Cv", "CH", "CTT", "CEE", "CBB"],
                                     [11/3]*3 + [2.5]*3, ["k3"]*3 + ["k2"]*3):
            for d, lab, st in [(m, "our solver", "-"), (e, "Enzo", "--")]:
                k = d[kk]; C = d[key]; s = d[key + "_std"]
                ax.fill_between(k, (C - s) * k**comp, (C + s) * k**comp, alpha=0.2)
                ax.loglog(k, C * k**comp, st, label=lab)
            k1, k2 = m["kfit"]
            ax.axvspan(k1, k2, color="0.9", zorder=0)
            ax.set_title(key[1:]); ax.set_xlabel("k / k_min")
        axs[0, 0].legend()
        fig.suptitle(f"N={n}, {label}")
        fig.tight_layout()
        fig.savefig(os.path.join(a.out, f"spectra_ax{axis}.png"), dpi=130)
        plt.close(fig)
    txt = "\n".join(lines)
    print(txt)
    open(os.path.join(a.out, "comparison.txt"), "w").write(txt + "\n")


if __name__ == "__main__":
    main()
