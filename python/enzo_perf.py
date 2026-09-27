"""Summarize Enzo's performance.out: where the time goes, and how fast per core.

    python3 python/enzo_perf.py <run>/performance.out [--skip 5] [--oneline]

For each cycle Enzo records every routine's time as mean/std/min/max across
MPI ranks, and on the Total line also the cells updated and Enzo's own
"cell updates / processor / sec".  Reported here, averaged over cycles after
--skip (start-up):

  MHDRK2                 the MHD solve itself (compute)
  SetBoundaryConditions  ghost-zone exchange (communication, on MPI runs)
  imbalance              slowest / fastest rank per cycle

Reading it: a rank that spends as long in MHDRK2 as one core on the whole box
means ranks share cores (bad binding).  A large SetBoundaryConditions share
means the pieces are too small for this many ranks.

Python 3.6 compatible (system python3 on the cluster).
"""
import argparse
import statistics


def parse(path):
    nproc, cycles, cur = None, [], None
    for line in open(path, errors="replace"):
        if line.startswith("# Starting performance log. MPI processes:"):
            nproc = int(line.rsplit(":", 1)[1])
            continue
        if line.startswith("Cycle_Number"):
            cur = {"cycle": int(line.split()[1])}
            cycles.append(cur)
            continue
        if cur is None or not line.strip() or line.startswith("#"):
            continue
        p = line.split()
        try:
            cur[p[0]] = [float(x) for x in p[1:]]
        except ValueError:
            pass
    return nproc or 1, cycles


def summarize(path, skip=5):
    nproc, cycles = parse(path)
    use = [c for c in cycles if "Total" in c][skip:]
    if not use:
        raise SystemExit("no timed cycles after --skip in %s" % path)
    mean = lambda key, i: statistics.mean(c[key][i] for c in use if key in c)
    tot_mean, tot_max, tot_min = mean("Total", 0), mean("Total", 3), mean("Total", 2)
    cells = use[-1]["Total"][4]
    out = dict(
        nproc=nproc, cycles=len(use), cells=cells, n=round(cells ** (1 / 3)),
        total_mean=tot_mean, total_max=tot_max,
        imbalance=statistics.mean(c["Total"][3] / max(c["Total"][2], 1e-30) for c in use),
        rk2=mean("MHDRK2", 0) if "MHDRK2" in use[0] else float("nan"),
        bc=mean("SetBoundaryConditions", 0) if "SetBoundaryConditions" in use[0] else float("nan"),
        io=mean("Group_WriteAllData", 0) if "Group_WriteAllData" in use[0] else 0.0,
        enzo_cu_per_proc=statistics.median(c["Total"][6] for c in use),
    )
    # wall time of a cycle is set by the slowest rank
    out["cu_per_core"] = cells / (tot_max * nproc)
    out["cu_total"] = cells / tot_max
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("perf")
    ap.add_argument("--skip", type=int, default=5)
    ap.add_argument("--oneline", action="store_true")
    a = ap.parse_args()
    s = summarize(a.perf, a.skip)
    if a.oneline:
        print("%5d %6d %6d %10.3e %10.3e %7.1f%% %7.1f%% %6.2f" % (
            s["n"], s["nproc"], s["cycles"], s["cu_per_core"], s["cu_total"],
            100 * s["rk2"] / s["total_mean"], 100 * s["bc"] / s["total_mean"], s["imbalance"]))
        return
    print("file                 : %s" % a.perf)
    print("box / ranks          : %d^3 on %d MPI ranks, %d^3 cells per rank" % (
        s["n"], s["nproc"], round((s["cells"] / s["nproc"]) ** (1 / 3))))
    print("cycles averaged      : %d (after skipping %d)" % (s["cycles"], a.skip))
    print("time per cycle       : %.4f s mean rank, %.4f s slowest rank (imbalance %.2f)" % (
        s["total_mean"], s["total_max"], s["imbalance"]))
    print("  MHD solve (MHDRK2) : %.4f s per rank  (%.0f%%)" % (s["rk2"], 100 * s["rk2"] / s["total_mean"]))
    print("  ghost exchange     : %.4f s per rank  (%.0f%%)" % (s["bc"], 100 * s["bc"] / s["total_mean"]))
    print("  output             : %.4f s per rank" % s["io"])
    print("throughput           : %.3e cell-updates/s per core   (%.3e total)" % (
        s["cu_per_core"], s["cu_total"]))
    print("  MHD solve per core : %.3e cell-updates/s (compute only)" % (
        s["cells"] / (s["rk2"] * s["nproc"])))


if __name__ == "__main__":
    main()
