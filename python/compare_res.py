"""Compare the same measurements across resolutions (and against the paper)."""
import sys, os, glob, json, argparse
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paper_values as PV
from figures import load

KEYS = [("al_Crho", "alpha_rho"), ("al_Cv", "alpha_v"), ("al_CH", "alpha_H"),
        ("al_CTT", "alpha_TT"), ("al_CEE", "alpha_EE"), ("al_CBB", "alpha_BB")]

def fits(runs):
    MS = np.array([float(r["Ms"]) for r in runs])
    MA = np.array([float(r["Ma"]) for r in runs])
    A = np.column_stack([np.ones_like(MS), MS, MA])
    out = {}
    for key, nm in KEYS:
        y = np.array([float(r[key]) for r in runs])
        c, *_ = np.linalg.lstsq(A, y, rcond=None)
        out[nm] = c
    rat = np.array([float(r["bp_BB_EE"]) for r in runs])
    rte = np.array([float(r["rTE"]) for r in runs])
    hi = MS > 4
    out["_hi"] = (rat[hi].mean() if hi.any() else np.nan,
                  rte[hi].mean() if hi.any() else np.nan, int(hi.sum()))
    out["_n"] = int(runs[0]["n"]); out["_nruns"] = len(runs)
    out["_win"] = tuple(runs[0]["kfit"])
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("roots", nargs="+")
    ap.add_argument("--out", default="figures/convergence.txt")
    ap.add_argument("--axis", type=int, default=2)
    a = ap.parse_args()
    res = []
    for r in a.roots:
        runs = load(r, axis=a.axis)
        if runs: res.append((r, fits(runs)))
    lines = ["Resolution dependence of the Table-1 coefficients "
             "(alpha_q = a + b M_S + c M_A)", ""]
    hdr = "  quantity  coef |" + "".join(f"  N={f['_n']:<4d}" for _, f in res) + "   paper"
    lines += [hdr, "  " + "-" * (len(hdr) - 2)]
    for _, nm in KEYS:
        for i, cf in enumerate("abc"):
            row = f"  {nm if i==0 else '':<10s}{cf:>4s} |"
            for _, f in res:
                row += f" {f[nm][i]:8.3f}"
            row += f"   {PV.TABLE1[nm][cf]:7.3f}"
            lines.append(row)
        lines.append("")
    row = "  <BB/EE> M_S>4  |"
    for _, f in res: row += f" {f['_hi'][0]:8.3f}"
    lines.append(row + f"   {PV.HIGH_MS['A_BB/A_EE'][0]:7.3f}")
    row = "  <r_TE>  M_S>4  |"
    for _, f in res: row += f" {f['_hi'][1]:8.3f}"
    lines.append(row + f"   {PV.HIGH_MS['r_TE']:7.3f}")
    row = "  fit window     |"
    for _, f in res: row += f"  [{f['_win'][0]},{f['_win'][1]}]".ljust(9)
    lines.append(row + "   [4,25]")
    txt = "\n".join(lines)
    print(txt)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    open(a.out, "w").write(txt + "\n")

if __name__ == "__main__":
    main()
