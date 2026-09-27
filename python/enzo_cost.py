"""What an Enzo run actually cost, from the timing log the job script keeps.

    python python/enzo_cost.py <enzo_run_dir> [--boxes 10]

cost.log has one line per job segment (a run that crosses the wall limit has
several):  start_epoch end_epoch ntasks cycle_start cycle_end sim_time_end
"""
import argparse, os, re
from pathlib import Path


def last_cycle(log):
    c, t = None, None
    if not os.path.exists(log):
        return c, t
    for line in open(log, errors="replace"):
        m = re.search(r"CycleNumber\s*=\s*(\d+)", line)
        if m:
            c = int(m.group(1))
        m2 = re.search(r"^Time\s*=\s*([0-9.eE+-]+)", line.strip())
        if m2:
            t = float(m2.group(1))
    return c, t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--boxes", type=int, default=10,
                    help="also project the cost of this many boxes")
    a = ap.parse_args()
    run = Path(a.run)
    rows = [l.split() for l in open(run / "cost.log") if l.strip() and not l.startswith("#")]
    if not rows:
        raise SystemExit("cost.log is empty")
    setup = dict(l.split(None, 1) for l in open(run / "setup.txt") if l.strip())
    n = int(setup["n"]); tdyn = float(setup["tdyn"]); ntdyn = float(setup["ntdyn"])
    wall = sum(float(r[1]) - float(r[0]) for r in rows)
    core_s = sum((float(r[1]) - float(r[0])) * int(r[2]) for r in rows)
    cyc = [(int(r[3]), int(r[4])) for r in rows if r[3] != "-" and r[4] != "-"]
    cycles = sum(max(0, e - b) for b, e in cyc)
    t_end = float(rows[-1][5]) if rows[-1][5] != "-" else float("nan")
    # NOT the RunFinished file: Enzo writes it on every clean exit, including
    # a stop on StopCPUTime.  Completion is the simulation time reached.
    frac = min(1.0, t_end / (ntdyn * tdyn)) if t_end == t_end else float("nan")
    done = frac == frac and frac >= 1 - 1e-4

    print(f"run            : {run}")
    print(f"status         : {'FINISHED' if done else f'in progress, {100*frac:.1f}% of {ntdyn:g} t_dyn'}")
    print(f"segments       : {len(rows)} (ranks per segment: {', '.join(r[2] for r in rows)})")
    print(f"wall time      : {wall/3600:.2f} h  ({wall/86400:.2f} days)")
    ch = core_s / 3600
    print(f"core-hours     : {ch:,.0f}" if ch >= 10 else f"core-hours     : {ch:.3f}")
    if cycles and wall > 0 and core_s > 0:
        print(f"cycles         : {cycles:,}   ({cycles/wall*3600:,.0f} per hour)")
        zu = n**3 * cycles / core_s
        print(f"throughput     : {zu:.3e} zone-updates / s / core")
    elif cycles:
        # a segment shorter than the 1 s clock resolution (tiny test boxes)
        print(f"cycles         : {cycles:,}   (too short to time)")
    if not done and frac == frac and frac > 0 and core_s > 0:
        pt = core_s / 3600 / frac
        print(f"projected total: {(f'{pt:,.0f}' if pt >= 10 else f'{pt:.3f}')} core-hours, "
              f"{wall/3600/frac:.2f} h wall at these ranks")
    tot = core_s / 3600 / (frac if (not done and frac == frac and frac > 0) else 1.0)
    tt = a.boxes * tot
    print(f"\n{a.boxes} independent boxes like this: "
          f"{(f'{tt:,.0f}' if tt >= 10 else f'{tt:.2f}')} core-hours "
          f"-> 16 maps each = {16*a.boxes} maps from full-depth 4x4 tiles, one snapshot per box")


if __name__ == "__main__":
    main()
