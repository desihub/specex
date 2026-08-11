"""
Full-CCD C++ vs Python parity campaign: for a list of cameras, run the real
production C++ wrapper (`srun -n 20 desi_compute_psf --mpi`, true per-band
defaults baked in) and the Python port (`specex.specex`, auto-detected
per-band settings, default mixed precision, --workers-per-gpu 5) on the
*whole* CCD (20 bundles), then compare:
  - wall time for both
  - X/Y trace RMS between the two merged output PSFs, across all 500 fibers
    (broken fibers excluded)
  - wavelength-residual-vs-line-list RMS/std for both pipelines, using C++'s
    own final per-bundle spot selection (cppspots_pass4.txt, written as a
    side effect of desi_compute_psf's per-bundle desi_psf_fit calls) as the
    common measurement set -- same methodology as bundle_parity_suite.py,
    extended to all 20 bundles.

Usage:
  python testing/full_ccd_campaign.py --night 20260401 --expid 00344649 \
      --cameras b5,b4,b2,r3,r5,r1,z1,z6,z9 --outdir /pscratch/.../full_ccd
"""
import os
import sys
import time
import glob
import argparse
import subprocess

import numpy as np
import fitsio

CURRENT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(CURRENT_DIR, 'py'))
sys.path.insert(0, os.path.join(CURRENT_DIR, 'testing'))

from specex.math import Legendre1DPol
from select_test_case import parse_log_line


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


def load_traces(path):
    f = fitsio.FITS(path)
    hdr = f['PSF'].read_header()
    wmin, wmax = hdr['WAVEMIN'], hdr['WAVEMAX']
    xtrace = f['XTRACE'].read().astype(np.float64)
    ytrace = f['YTRACE'].read().astype(np.float64)
    traces = {}
    for fib in range(xtrace.shape[0]):
        traces[fib] = {
            'X_vs_W': Legendre1DPol(deg=xtrace.shape[1] - 1, xmin=wmin, xmax=wmax, coeff=xtrace[fib]),
            'Y_vs_W': Legendre1DPol(deg=ytrace.shape[1] - 1, xmin=wmin, xmax=wmax, coeff=ytrace[fib]),
        }
    return traces, float(wmin), float(wmax)


def load_spots(path):
    fibers, waves, ycs = [], [], []
    with open(path) as f:
        for line in f:
            p = line.strip().split(',')
            if len(p) < 4:
                continue
            fibers.append(int(float(p[0]))); waves.append(float(p[1])); ycs.append(float(p[3]))
    return np.array(fibers), np.array(waves), np.array(ycs)


def start_cpp_full(case, out_fits, log_path):
    cmd = ["srun", "-n", "20", "desi_compute_psf", "--mpi",
           "--input-image", case['image'],
           "--input-psf", case['input_psf'],
           "-o", out_fits, "--extra=--debug-spots"]
    if case['broken_fibers']:
        cmd += ["--broken-fibers", case['broken_fibers']]
    lf = open(log_path, 'w')
    t0 = time.time()
    proc = subprocess.Popen(cmd, stdout=lf, stderr=subprocess.STDOUT)
    return proc, lf, t0


def start_py_full(case, out_fits, log_path):
    cmd = [sys.executable, "-m", "specex.specex",
           "-a", case['image'],
           "--in-psf", case['input_psf'],
           "--out-psf", out_fits,
           "--gpu", "4", "--workers-per-gpu", "5"]
    if case['broken_fibers']:
        cmd += ["--broken-fibers", case['broken_fibers']]
    lf = open(log_path, 'w')
    t0 = time.time()
    proc = subprocess.Popen(cmd, stdout=lf, stderr=subprocess.STDOUT)
    return proc, lf, t0


def run_cpp_and_py_concurrent(case, cpp_fits, py_fits, cpp_log, py_log):
    """C++ (CPU-only, MPI) and Python (GPU-only, mixed precision) don't
    contend for the same resources, so launch both at once instead of
    sequentially -- roughly halves per-camera wall time vs run-then-run.

    Poll both independently rather than cpp_proc.wait() then py_proc.wait()
    in sequence -- now that the persistent JAX compilation cache + power-
    of-2 shape padding (see porting-notes.md) often make Python finish
    *before* C++, a strict wait()-then-wait() ordering silently caps t_py
    at whatever t_cpp was: py_proc.wait() on an already-finished process
    returns instantly, but time.time() - t0_py at that point measures "how
    long since Python started until *C++* finished," not Python's actual
    finish time. Confirmed this was happening: t_py in a recent run matched
    t_cpp to the decimal in every camera, while each Python process's own
    internally-printed "Total CCD Fit Time" was 3-4x shorter.
    """
    cpp_proc, cpp_lf, t0_cpp = start_cpp_full(case, cpp_fits, cpp_log)
    py_proc, py_lf, t0_py = start_py_full(case, py_fits, py_log)
    t_cpp = t_py = rc_cpp = rc_py = None
    while t_cpp is None or t_py is None:
        if t_cpp is None:
            rc = cpp_proc.poll()
            if rc is not None:
                rc_cpp = rc; t_cpp = time.time() - t0_cpp; cpp_lf.close()
        if t_py is None:
            rc = py_proc.poll()
            if rc is not None:
                rc_py = rc; t_py = time.time() - t0_py; py_lf.close()
        if t_cpp is None or t_py is None:
            time.sleep(0.2)
    return t_cpp, rc_cpp, t_py, rc_py


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--night", default="20260401")
    ap.add_argument("--expid", default="00344649")
    ap.add_argument("--cameras", help="comma list, e.g. b5,r3,z1 (fixed night/expid mode)")
    ap.add_argument("--cases-file", help="JSON-lines file of resolved cases "
                     "(camera,night,expid,image,input_psf,broken_fibers) -- "
                     "for multi-night random campaigns, produced by random_case_picker.py")
    ap.add_argument("--outdir", default="/pscratch/sd/c/cdwarner/specex/testing/full_ccd_campaign")
    ap.add_argument("--results", default=None)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    results_path = args.results or os.path.join(args.outdir, "full_ccd_results.txt")
    manifest_path = os.path.join(args.outdir, "input_files_manifest.txt")

    hdr = (f"{'case':<14} {'night':<9} {'expid':<9} {'nspots_cpp':>10} {'nspots_py':>9} {'xrms_px':>8} {'yrms_px':>8} "
           f"{'wrms_cpp_A':>10} {'wrms_py_A':>9} {'wstd_cpp_A':>10} {'wstd_py_A':>9} "
           f"{'t_cpp_s':>8} {'t_py_s':>7}")
    print(hdr, flush=True)
    new_file = not os.path.exists(results_path)
    with open(results_path, 'a') as rf:
        if new_file:
            rf.write(hdr + "\n")

    if args.cases_file:
        import json
        with open(args.cases_file) as f:
            queue = [json.loads(line) for line in f if line.strip()]
    else:
        queue = []
        for cam in args.cameras.split(","):
            case = find_case(args.night, args.expid, cam)
            if case is None:
                print(f"{cam:<7} CASE NOT FOUND", flush=True)
                continue
            case['camera'] = cam; case['night'] = args.night; case['expid'] = args.expid
            queue.append(case)

    for case in queue:
        cam, night, expid = case['camera'], case['night'], case['expid']
        tag = f"{cam}@{night}"
        broken = set(int(x) for x in case['broken_fibers'].split(',') if x.strip()) if case['broken_fibers'] else set()

        with open(manifest_path, 'a') as mf:
            mf.write(f"{case['image']}\n{case['input_psf']}\n")

        cpp_fits = os.path.join(args.outdir, f"cpp-{tag}-{expid}.fits")
        py_fits = os.path.join(args.outdir, f"py-{tag}-{expid}.fits")

        t_cpp, rc_cpp, t_py, rc_py = run_cpp_and_py_concurrent(
            case, cpp_fits, py_fits,
            os.path.join(args.outdir, f"cpp-{tag}.log"),
            os.path.join(args.outdir, f"py-{tag}.log"))
        if rc_cpp != 0:
            print(f"{tag:<14} CPP FAILED rc={rc_cpp}", flush=True)
            continue
        if rc_py != 0:
            print(f"{tag:<14} PY FAILED rc={rc_py}", flush=True)
            continue

        cpp_tr, wmin, wmax = load_traces(cpp_fits)
        py_tr, _, _ = load_traces(py_fits)
        grid = np.linspace(wmin, wmax, 100)
        dx, dy = [], []
        for fib in range(500):
            if fib in broken:
                continue
            dx.append(np.array(cpp_tr[fib]['X_vs_W'].value(grid)) - np.array(py_tr[fib]['X_vs_W'].value(grid)))
            dy.append(np.array(cpp_tr[fib]['Y_vs_W'].value(grid)) - np.array(py_tr[fib]['Y_vs_W'].value(grid)))
        dx, dy = np.array(dx), np.array(dy)
        xr, yr = float(np.sqrt(np.mean(dx**2))), float(np.sqrt(np.mean(dy**2)))

        spot_files = sorted(glob.glob(os.path.join(os.path.dirname(cpp_fits), f"cpp-{tag}-{expid}_*.cppspots_pass4.txt")))
        if not spot_files:
            # desi_compute_psf writes per-bundle debug files next to its own
            # intermediate output naming, not next to our renamed merged
            # cpp_fits -- search the desispec exposures dir as a fallback.
            spot_files = sorted(glob.glob(os.path.join(os.path.dirname(case['input_psf']), f"fit-psf-{cam}-{expid}_*.cppspots_pass4.txt")))

        n_cpp = sum(1 for fp in spot_files for _ in open(fp)) if spot_files else -1
        if spot_files:
            all_resid_cpp, all_resid_py = [], []
            for fp in spot_files:
                fibers, waves_true, ycs = load_spots(fp)
                for fib in np.unique(fibers):
                    mask = fibers == fib
                    all_resid_cpp.append(np.array(cpp_tr[fib]['Y_vs_W'].invert(ycs[mask])) - waves_true[mask])
                    all_resid_py.append(np.array(py_tr[fib]['Y_vs_W'].invert(ycs[mask])) - waves_true[mask])
            resid_cpp = np.concatenate(all_resid_cpp); resid_py = np.concatenate(all_resid_py)
            wrms_cpp, wstd_cpp = float(np.sqrt(np.mean(resid_cpp**2))), float(resid_cpp.std())
            wrms_py, wstd_py = float(np.sqrt(np.mean(resid_py**2))), float(resid_py.std())
        else:
            wrms_cpp = wrms_py = wstd_cpp = wstd_py = float('nan')

        # fit_ccd_native's 20 workers all share the same psf.output_psf_path
        # (the final merged out_fits), so per-bundle .pyspots.txt debug files
        # race/overwrite each other -- not reliable for a total spot count.
        # Sum the per-bundle "Iterative spot selection took ... (N spots)"
        # lines from the run's own stdout log instead.
        import re
        n_py = 0
        with open(os.path.join(args.outdir, f"py-{tag}.log")) as lf:
            for line in lf:
                m = re.search(r"Iterative spot selection took [\d.]+s \((\d+) spots\)", line)
                if m:
                    n_py += int(m.group(1))

        row = (f"{tag:<14} {night:<9} {expid:<9} {n_cpp:>10} {n_py:>9} {xr:>8.4f} {yr:>8.4f} "
               f"{wrms_cpp:>10.4f} {wrms_py:>9.4f} {wstd_cpp:>10.4f} {wstd_py:>9.4f} "
               f"{t_cpp:>8.1f} {t_py:>7.1f}")
        print(row, flush=True)
        with open(results_path, 'a') as rf:
            rf.write(row + "\n")


if __name__ == "__main__":
    main()
