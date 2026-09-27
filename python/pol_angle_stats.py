"""How ordered is the projected field?  Polarization-angle statistics per run.

    python3 python/pol_angle_stats.py runs/n128/ms5_ma1 runs/n128/ms5_ma2 ... [--tiles 4]

For the last snapshots of each run (line of sight perpendicular to B0):
  Q<0, U<0     fraction of pixels with negative Stokes Q / U
  spread       circular standard deviation of the polarization angle over the map
  tile spread  median of the same over each of tiles x tiles sub-maps
A strong mean field (low M_A) gives a small spread and Q of one sign; the
overall sign is set by the orientation of B0 relative to the map axes.
"""
import argparse, glob, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analysis import read_snapshot, project, finite_snapshot


def spread_deg(Q, U):
    R = np.abs(np.mean(np.exp(1j * np.arctan2(U, Q))))      # resultant of 2*psi
    return float(np.degrees(0.5 * np.sqrt(-2 * np.log(max(R, 1e-12)))))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--tiles", type=int, default=4)
    ap.add_argument("--snaps", type=int, default=4)
    a = ap.parse_args()
    print(f"{'run':<14}{'N':>5}{'M_S':>6}{'M_A':>6}{'Q<0':>7}{'U<0':>7}{'spread':>9}{'tile spread':>13}")
    for run in a.runs:
        snaps = sorted(glob.glob(os.path.join(run, "snap_*.bin")))[-a.snaps:]
        if not snaps:
            print(f"{os.path.basename(run.rstrip('/')):<14}  (no snapshots)")
            continue
        qn, un, sp, spt = [], [], [], []
        for sf in snaps:
            s = read_snapshot(sf)
            if not finite_snapshot(s):
                continue
            T, Q, U = project(s["rho"], s["H"], axis=2)
            qn.append((Q < 0).mean()); un.append((U < 0).mean()); sp.append(spread_deg(Q, U))
            m = s["n"] // a.tiles
            spt += [spread_deg(Q[i*m:(i+1)*m, j*m:(j+1)*m], U[i*m:(i+1)*m, j*m:(j+1)*m])
                    for i in range(a.tiles) for j in range(a.tiles)]
        if not qn:
            print(f"{os.path.basename(run.rstrip('/')):<14}  (SKIPPED: no finite snapshot)")
            continue
        print(f"{os.path.basename(run.rstrip('/')):<14}{s['n']:5d}{float(s['meta']['ms']):6.2f}"
              f"{float(s['meta']['ma']):6.2f}{100*np.mean(qn):6.0f}%{100*np.mean(un):6.0f}%"
              f"{np.mean(sp):8.1f}°{np.median(spt):11.1f}°")


if __name__ == "__main__":
    main()
