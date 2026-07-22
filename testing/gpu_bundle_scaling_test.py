"""
For a representative full-CCD case per band, sweep --workers-per-gpu on a
SINGLE GPU (--gpu 1) with --first-bundle 0 --last-bundle N-1 (so N bundles
launch as one wave, all packed onto GPU 0) to find how many concurrent
bundle-fit workers a single A100 can hold under default mixed precision
before OOM/failure or a real compute-contention slowdown.

Polls nvidia-smi for GPU 0 in a background thread during each run to record
peak memory used, independent of whatever the process itself reports.

Usage:
  python testing/gpu_bundle_scaling_test.py --night 20260401 --expid 00344649 \
      --cameras b5,r5,z6 --n-list 5,8,10,13,16,20 \
      --outdir /pscratch/.../gpu_scaling
"""
import os
import sys
import time
import threading
import argparse
import subprocess

CURRENT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(CURRENT_DIR, 'testing'))
from select_test_case import parse_log_line
import glob


def find_case(night, expid, camera):
    log_dir = f"/global/cfs/cdirs/desi/spectro/redux/matterhorn/run/scripts/night/{night}"
    for log in glob.glob(os.path.join(log_dir, "arc*.log")):
        with open(log) as f:
            for line in f:
                if "desi_compute_psf" in line and expid in line and f"-{camera}-" in line:
                    case = parse_log_line(line)
                    if case:
                        return case
    return None


def poll_gpu_mem(gpu_index, stop_event, peak_holder):
    peak = 0
    while not stop_event.is_set():
        try:
            out = subprocess.check_output(
                ["nvidia-smi", f"--query-gpu=memory.used", "--format=csv,noheader,nounits", "-i", str(gpu_index)],
                text=True, timeout=5)
            used = int(out.strip().split("\n")[0])
            peak = max(peak, used)
        except Exception:
            pass
        time.sleep(0.4)
    peak_holder[0] = peak


def run_one(case, n, out_fits, log_path, gpu_index=0):
    if n < 1:
        return None
    last_bundle = n - 1
    cmd = [sys.executable, "-m", "specex.specex",
           "-a", case['image'], "--in-psf", case['input_psf'],
           "--out-psf", out_fits,
           "--first-bundle", "0", "--last-bundle", str(last_bundle),
           "--gpu", "1", "--workers-per-gpu", str(n)]
    if case['broken_fibers']:
        cmd += ["--broken-fibers", case['broken_fibers']]
    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = str(gpu_index)

    stop_event = threading.Event()
    peak_holder = [0]
    t = threading.Thread(target=poll_gpu_mem, args=(gpu_index, stop_event, peak_holder))
    t.start()
    t0 = time.time()
    with open(log_path, 'w') as lf:
        r = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, env=env)
    dt = time.time() - t0
    stop_event.set()
    t.join()
    return dt, r.returncode, peak_holder[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--night", default="20260401")
    ap.add_argument("--expid", default="00344649")
    ap.add_argument("--cameras", required=True, help="one representative camera per band, e.g. b5,r5,z6")
    ap.add_argument("--n-list", default="5,8,10,13,16,20")
    ap.add_argument("--outdir", default="/pscratch/sd/c/cdwarner/specex/testing/gpu_scaling")
    ap.add_argument("--results", default=None)
    ap.add_argument("--gpu-index", type=int, default=0)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    results_path = args.results or os.path.join(args.outdir, "gpu_scaling_results.txt")
    n_list = [int(x) for x in args.n_list.split(",")]

    hdr = f"{'camera':<7} {'n_workers':>9} {'result':<10} {'peak_mem_MiB':>13} {'wall_s':>8}"
    print(hdr, flush=True)
    new_file = not os.path.exists(results_path)
    with open(results_path, 'a') as rf:
        if new_file:
            rf.write(hdr + "\n")

    for cam in args.cameras.split(","):
        case = find_case(args.night, args.expid, cam)
        if case is None:
            print(f"{cam:<7} CASE NOT FOUND", flush=True)
            continue
        for n in n_list:
            out_fits = os.path.join(args.outdir, f"py-{cam}-n{n}.fits")
            log_path = os.path.join(args.outdir, f"py-{cam}-n{n}.log")
            res = run_one(case, n, out_fits, log_path, gpu_index=args.gpu_index)
            if res is None:
                continue
            dt, rc, peak = res
            status = "OK" if rc == 0 else f"FAIL(rc={rc})"
            row = f"{cam:<7} {n:>9} {status:<10} {peak:>13} {dt:>8.1f}"
            print(row, flush=True)
            with open(results_path, 'a') as rf:
                rf.write(row + "\n")
            if rc != 0:
                # once it fails, larger N will too -- stop sweeping this camera
                break


if __name__ == "__main__":
    main()
