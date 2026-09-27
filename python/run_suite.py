"""Launch the (M_S, M_A) suite of driven-turbulence runs with a simple pool."""
import os, sys, subprocess, itertools, argparse, time

MS_LIST = [0.5, 1, 2, 3, 4, 5, 6]
MA_LIST = [0.5, 1, 2]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=64)
    ap.add_argument("--root", default=None)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--nsnap", type=int, default=21)
    ap.add_argument("--ntdyn", type=float, default=10.0)
    ap.add_argument("--edotfac", type=float, default=0.42)
    ap.add_argument("--ms", type=float, nargs="*", default=None)
    ap.add_argument("--ma", type=float, nargs="*", default=None)
    ap.add_argument("--dry", action="store_true")
    ap.add_argument("--binary", default="./src/mhd")
    ap.add_argument("--tsnap0", type=float, default=5.0)
    ap.add_argument("--upsample-from", default=None,
                    help="root of a lower-resolution suite; each run starts from "
                         "the last snapshot of the matching run there")
    a = ap.parse_args()
    root = a.root or f"runs/n{a.n}"
    os.makedirs(root, exist_ok=True)
    ms_list = a.ms if a.ms else MS_LIST
    ma_list = a.ma if a.ma else MA_LIST
    todo = []
    for i, (ms, ma) in enumerate(itertools.product(MS_LIST, MA_LIST)):
        if ms not in ms_list or ma not in ma_list:
            continue
        tag = f"ms{ms:g}_ma{ma:g}"
        out = os.path.join(root, tag)
        if os.path.exists(os.path.join(out, f"snap_{a.nsnap-1:04d}.bin")):
            print("skip (done):", out); continue
        cmd = [a.binary, f"n={a.n}", f"ms={ms}", f"ma={ma}",
               f"ntdyn={a.ntdyn}", f"nsnap={a.nsnap}", f"edotfac={a.edotfac}",
               f"tsnap0={a.tsnap0}", f"seed={1000+i}", f"out={out}"]
        if a.upsample_from:
            import glob as _g
            src = sorted(_g.glob(os.path.join(a.upsample_from, tag, "snap_*.bin")))
            if not src:
                print(f"!! no source snapshot for {tag} in {a.upsample_from}"); continue
            cmd.append(f"upsample={src[-1]}")
        todo.append((tag, cmd, out))
    print(f"{len(todo)} runs, {a.jobs} at a time, {a.threads} threads each", flush=True)
    if a.dry:
        for t, c, _ in todo: print(" ".join(c))
        return
    env = dict(os.environ); env["OMP_NUM_THREADS"] = str(a.threads)
    running, t0 = [], time.time()
    while todo or running:
        while todo and len(running) < a.jobs:
            tag, cmd, out = todo.pop(0)
            os.makedirs(out, exist_ok=True)
            lf = open(os.path.join(out, "stdout.log"), "w")
            p = subprocess.Popen(cmd, stdout=lf, stderr=subprocess.STDOUT, env=env)
            running.append((tag, p, lf, time.time()))
            print(f"[{time.time()-t0:7.0f}s] start {tag} (pid {p.pid})", flush=True)
        time.sleep(5)
        for r in running[:]:
            tag, p, lf, ts = r
            if p.poll() is not None:
                lf.close(); running.remove(r)
                print(f"[{time.time()-t0:7.0f}s] done  {tag} rc={p.returncode} "
                      f"({time.time()-ts:.0f}s)", flush=True)
    print(f"suite complete in {time.time()-t0:.0f}s", flush=True)

if __name__ == "__main__":
    main()
