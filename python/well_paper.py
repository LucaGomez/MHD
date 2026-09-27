"""The Well MHD_256 boxes measured the way Stalpes, Collins & Huffenberger (2024)
measure theirs, next to the paper's Table 1 prediction at the measured Mach
numbers.

    python3 python/well_paper.py --well runs_well [--out figures/well_paper]

Uses <well>/<pair>/analysis_ax2.npz (line of sight perpendicular to the mean
field; periodic 256^2 faces, fit window k = [3, 13]) written by
well_stream.py --summarize-only.  The paper's runs span M_S 0.5-7 and M_A
0.5-2; the weak-field Well runs sit at M_A ~ 8, so their "prediction" is an
extrapolation of the linear fit and is flagged.
"""
import argparse, glob, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paper_values as PV

QUANT = [("alpha_rho", "al_Crho"), ("alpha_v", "al_Cv"), ("alpha_H", "al_CH"),
         ("alpha_TT", "al_CTT"), ("alpha_EE", "al_CEE"), ("alpha_BB", "al_CBB")]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--well", default="runs_well")
    ap.add_argument("--out", default="figures/well_paper")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    runs = []
    for f in sorted(glob.glob(os.path.join(a.well, "*", "analysis_ax2.npz"))):
        d = dict(np.load(f, allow_pickle=True))
        d["tag"] = os.path.basename(os.path.dirname(f))
        runs.append(d)
    runs.sort(key=lambda d: (float(d["Ma"]), float(d["Ms"])))
    L = ["The Well MHD_256 (Cho-ENO, 256^3, solenoidal driving) analysed like Stalpes+2024",
         "line of sight perpendicular to B0; slopes over k=[3,13]; 'pred' = paper Table 1 at the",
         "measured (M_S, M_A); '*' = outside the paper's range (M_S<=7, M_A<=2): extrapolation", ""]
    hdr = f"{'pair':<13}{'N':>3}{'M_S':>6}{'M_A':>6} " + "".join(
        f"{q[6:]:>13}" for q, _ in QUANT) + f"{'BB/EE':>8}{'r_TE':>7}{'r_TB':>7}"
    L.append(hdr + "       (each slope: measured / pred)")
    for d in runs:
        ms, ma = float(d["Ms"]), float(d["Ma"])
        flag = "*" if (ms > 7.5 or ma > 2.5) else " "
        row = f"{d['tag']:<13}{int(d['nsnap']):3d}{ms:6.2f}{ma:6.2f}{flag}"
        for q, key in QUANT:
            p = PV.TABLE1[q]
            pred = p["a"] + p["b"] * ms + p["c"] * ma
            row += f"  {float(d[key]):5.2f}/{pred:5.2f}"
        row += f"{float(d['bp_BB_EE']):8.2f}{float(d['rTE']):7.2f}{float(d['rTB']):7.3f}"
        L.append(row)
    L += ["", "paper at M_S > 4: BB/EE (amplitude ratio) 0.55 +- 0.07, r_TE ~ 0.3;  "
          "Planck: alpha_EE -2.42, alpha_BB -2.54, BB/EE 0.53, r_TE 0.355"]

    # M_S dependence at fixed (strong) field, the only family inside the paper's M_A range
    for fam, sel in [("M_A ~ 0.8 family", lambda d: float(d["Ma"]) < 2),
                     ("M_A ~ 8 family", lambda d: float(d["Ma"]) >= 2)]:
        rs = [d for d in runs if sel(d)]
        if len(rs) < 3:
            continue
        ms = np.array([float(d["Ms"]) for d in rs])
        L.append(f"\n{fam}: linear fit alpha = a + b M_S over M_S {ms.min():.1f}-{ms.max():.1f}"
                 f" (paper b in brackets; paper a + c*M_A at the family's mean M_A)")
        mam = np.mean([float(d["Ma"]) for d in rs])
        for q, key in QUANT:
            y = np.array([float(d[key]) for d in rs])
            b, a0 = np.polyfit(ms, y, 1)
            p = PV.TABLE1[q]
            L.append(f"   {q:<10} a={a0:6.2f} [{p['a'] + p['c'] * mam:6.2f}]   "
                     f"b={b:6.3f} [{p['b']:6.3f}]   scatter={np.std(y - (a0 + b * ms)):.3f}")
    txt = "\n".join(L)
    print(txt)
    open(os.path.join(a.out, "table.txt"), "w").write(txt + "\n")


if __name__ == "__main__":
    main()
