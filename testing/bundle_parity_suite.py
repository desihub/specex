"""
Run C++ and Python PSF fits for a list of (camera, bundle) cases and compare:
  - X/Y trace RMS between the two output PSFs (100-pt wave grid per fiber)
  - wavelength-residual RMS vs the lamp line list for BOTH pipelines
    (invert each pipeline's fitted Y_vs_W trace at the C++ final spots' yc,
     compare to the true line wavelength; report raw RMS and mean-subtracted
     scatter so the shared ~0.5A calibration offset does not dominate)
  - final spot counts and wall times

Usage:
  python testing/bundle_parity_suite.py --night 20260401 --expid 00344649 \
      --cases z8:5,z8:0,z8:18,z5:5 [--outdir /pscratch/.../multi] [--skip-cpp]

Assumes the instrumented C++ build (writes .cppspots_pass4.txt) and the
Python driver are both importable per env_setup.sh.
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
            'X_vs_W': Legendre1DPol(deg=xtrace.shape[1]-1, xmin=wmin, xmax=wmax, coeff=xtrace[fib]),
            'Y_vs_W': Legendre1DPol(deg=ytrace.shape[1]-1, xmin=wmin, xmax=wmax, coeff=ytrace[fib]),
        }
    return traces, float(wmin), float(wmax)


def trace_rms(cpp_fits, py_fits, bundle_id, broken):
    cpp_tr, wmin, wmax = load_traces(cpp_fits)
    py_tr, _, _ = load_traces(py_fits)
    fmin, fmax = bundle_id * 25, (bundle_id + 1) * 25 - 1
    grid = np.linspace(wmin, wmax, 100)
    dx, dy = [], []
    for fib in range(fmin, fmax + 1):
        if fib in broken:
            continue
        dx.append(np.array(cpp_tr[fib]['X_vs_W'].value(grid)) - np.array(py_tr[fib]['X_vs_W'].value(grid)))
        dy.append(np.array(cpp_tr[fib]['Y_vs_W'].value(grid)) - np.array(py_tr[fib]['Y_vs_W'].value(grid)))
    dx, dy = np.array(dx), np.array(dy)
    return float(np.sqrt(np.mean(dx**2))), float(np.sqrt(np.mean(dy**2)))


def load_spots(path):
    fibers, waves, ycs = [], [], []
    with open(path) as f:
        for line in f:
            p = line.strip().split(',')
            if len(p) < 4:
                continue
            fibers.append(int(float(p[0])))
            waves.append(float(p[1]))
            ycs.append(float(p[3]))
    return np.array(fibers), np.array(waves), np.array(ycs)


def wave_residual_stats(fits_path, fibers, waves_true, ycs):
    traces, _, _ = load_traces(fits_path)
    resid = np.zeros(len(fibers))
    for fib in np.unique(fibers):
        mask = fibers == fib
        resid[mask] = np.array(traces[fib]['Y_vs_W'].invert(ycs[mask])) - waves_true[mask]
    rms = float(np.sqrt(np.mean(resid**2)))
    return rms, float(resid.mean()), float(resid.std())


def run_cpp(case, bundle_id, out_fits, log_path):
    fmin, fmax = bundle_id * 25, (bundle_id + 1) * 25 - 1
    cmd = ["desi_psf_fit",
           "-a", case['image'],
           "--in-psf", case['input_psf'],
           "--lamp-lines", os.path.join(CURRENT_DIR, "py/specex/data/specex_linelist_desi.txt"),
           "--out-psf", out_fits,
           "--first-bundle", str(bundle_id), "--last-bundle", str(bundle_id),
           "--first-fiber", str(fmin), "--last-fiber", str(fmax),
           "--legendre-deg-wave", "3", "--fit-continuum"]
    if case['broken_fibers']:
        cmd += ["--broken-fibers", case['broken_fibers']]
    t0 = time.time()
    with open(log_path, 'w') as lf:
        r = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
    return time.time() - t0, r.returncode


def run_py(case, bundle_id, out_fits, log_path):
    fmin, fmax = bundle_id * 25, (bundle_id + 1) * 25 - 1
    cmd = [sys.executable, "-m", "specex.specex",
           "-a", case['image'],
           "--in-psf", case['input_psf'],
           "--lamp-lines", os.path.join(CURRENT_DIR, "py/specex/data/specex_linelist_desi.txt"),
           "--out-psf", out_fits,
           "--first-bundle", str(bundle_id), "--last-bundle", str(bundle_id),
           "--first-fiber", str(fmin), "--last-fiber", str(fmax),
           "--legendre-deg-wave", "3", "--fit-continuum", "--gpu", "1"]
    if case['broken_fibers']:
        cmd += ["--broken-fibers", case['broken_fibers']]
    t0 = time.time()
    with open(log_path, 'w') as lf:
        r = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
    return time.time() - t0, r.returncode


def count_lines(path):
    if not os.path.exists(path):
        return -1
    with open(path) as f:
        return sum(1 for line in f if line.strip())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--night", default="20260401")
    ap.add_argument("--expid", default="00344649")
    ap.add_argument("--cases", required=True, help="comma list of cam:bundle, e.g. z8:5,z5:5")
    ap.add_argument("--outdir", default="/pscratch/sd/c/cdwarner/specex/testing/multi")
    ap.add_argument("--skip-cpp", action="store_true", help="reuse existing C++ outputs")
    ap.add_argument("--skip-py", action="store_true", help="reuse existing Python outputs")
    ap.add_argument("--results", default=None, help="results table path (default <outdir>/parity_results.txt)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    results_path = args.results or os.path.join(args.outdir, "parity_results.txt")
    cases = [c.split(":") for c in args.cases.split(",")]

    hdr = (f"{'case':<8} {'nspots_cpp':>10} {'nspots_py':>9} {'xrms_px':>8} {'yrms_px':>8} "
           f"{'wrms_cpp_A':>10} {'wrms_py_A':>9} {'wstd_cpp_A':>10} {'wstd_py_A':>9} "
           f"{'t_cpp_s':>8} {'t_py_s':>7}")
    print(hdr, flush=True)
    new_file = not os.path.exists(results_path)
    with open(results_path, 'a') as rf:
        if new_file:
            rf.write(hdr + "\n")

    for cam, bstr in cases:
        bid = int(bstr)
        tag = f"{cam}:{bid}"
        case = find_case(args.night, args.expid, cam)
        if case is None:
            print(f"{tag:<8} CASE NOT FOUND", flush=True)
            continue
        broken = set(int(x) for x in case['broken_fibers'].split(',') if x.strip()) if case['broken_fibers'] else set()

        base = f"{cam}-{args.expid}_{bid:02d}"
        cpp_fits = os.path.join(args.outdir, f"cpp-{base}.fits")
        py_fits = os.path.join(args.outdir, f"py-{base}.fits")

        t_cpp = rc = 0
        if not (args.skip_cpp and os.path.exists(cpp_fits)):
            t_cpp, rc = run_cpp(case, bid, cpp_fits, os.path.join(args.outdir, f"cpp-{base}.log"))
            if rc != 0:
                print(f"{tag:<8} CPP FAILED rc={rc} (see cpp-{base}.log)", flush=True)
                continue
        t_py = 0
        if not (args.skip_py and os.path.exists(py_fits)):
            t_py, rc = run_py(case, bid, py_fits, os.path.join(args.outdir, f"py-{base}.log"))
            if rc != 0:
                print(f"{tag:<8} PY FAILED rc={rc} (see py-{base}.log)", flush=True)
                continue

        xr, yr = trace_rms(cpp_fits, py_fits, bid, broken)
        cpp_spots_file = cpp_fits.replace('.fits', '.cppspots_pass4.txt')
        py_spots_file = py_fits.replace('.fits', '.pyspots.txt')
        n_cpp, n_py = count_lines(cpp_spots_file), count_lines(py_spots_file)

        if os.path.exists(cpp_spots_file):
            fibers, waves_true, ycs = load_spots(cpp_spots_file)
            wrms_cpp, _, wstd_cpp = wave_residual_stats(cpp_fits, fibers, waves_true, ycs)
            wrms_py, _, wstd_py = wave_residual_stats(py_fits, fibers, waves_true, ycs)
        else:
            wrms_cpp = wrms_py = wstd_cpp = wstd_py = float('nan')

        row = (f"{tag:<8} {n_cpp:>10} {n_py:>9} {xr:>8.4f} {yr:>8.4f} "
               f"{wrms_cpp:>10.4f} {wrms_py:>9.4f} {wstd_cpp:>10.4f} {wstd_py:>9.4f} "
               f"{t_cpp:>8.1f} {t_py:>7.1f}")
        print(row, flush=True)
        with open(results_path, 'a') as rf:
            rf.write(row + "\n")


if __name__ == "__main__":
    main()
