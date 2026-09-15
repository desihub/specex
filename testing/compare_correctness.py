#!/usr/bin/env python3
"""Compare C++ (desi_proc -n101, staged real preproc via stage_preproc.py)
vs Python (run_night.py --backend python) trace output across a list of
nights: per-camera and per-night xrms/yrms (px), same methodology as
full_ccd_campaign.py (Legendre trace polys evaluated on a 100-point
wavelength grid, broken fibers excluded).

Usage:
  python testing/compare_correctness.py \
      --py-dir /pscratch/.../cpptest/redux/cdwarner \
      --cpp-base /pscratch/.../cpptest/desiproc2 \
      --nights 20260316:00342128,20260401:00344649,20220120:00119496
"""
import os
import sys
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from full_ccd_campaign import load_traces
from run_night import find_cases, ALL_CAMERAS

ap = argparse.ArgumentParser()
ap.add_argument("--py-dir", required=True, help="dir containing fit-psf-<cam>-<expid>.fits from --backend python")
ap.add_argument("--cpp-base", required=True, help="dir containing <night>_<expid>/redux/<specprod>/exposures/<night>/<expid>/fit-psf-<cam>-<expid>.fits")
ap.add_argument("--specprod", default=os.environ.get("USER", "specex"))
ap.add_argument("--nights", required=True, help="comma-separated night:expid pairs, e.g. 20260316:00342128,20260401:00344649")
args = ap.parse_args()

NIGHTS = [tuple(pair.split(":")) for pair in args.nights.split(",")]

all_night_xr = []
all_night_yr = []

for night, expid in NIGHTS:
    cases = find_cases(night, expid, ALL_CAMERAS)
    cpp_dir = os.path.join(args.cpp_base, f"{night}_{expid}", "redux", args.specprod,
                            "exposures", night, expid)
    cam_xr, cam_yr = [], []
    print(f"\n=== {night}/{expid} ===")
    for cam in ALL_CAMERAS:
        py_fits = os.path.join(args.py_dir, f"fit-psf-{cam}-{expid}.fits")
        cpp_fits = os.path.join(cpp_dir, f"fit-psf-{cam}-{expid}.fits")
        if not (os.path.exists(py_fits) and os.path.exists(cpp_fits)):
            print(f"  {cam}: MISSING (py={os.path.exists(py_fits)} cpp={os.path.exists(cpp_fits)})")
            continue
        broken = set()
        bf = cases.get(cam, {}).get("broken_fibers", "")
        if bf:
            broken = set(int(x) for x in bf.split(",") if x.strip())

        py_tr, wmin, wmax = load_traces(py_fits)
        cpp_tr, _, _ = load_traces(cpp_fits)
        grid = np.linspace(wmin, wmax, 100)
        dx, dy = [], []
        for fib in range(500):
            if fib in broken:
                continue
            dx.append(np.array(cpp_tr[fib]['X_vs_W'].value(grid)) - np.array(py_tr[fib]['X_vs_W'].value(grid)))
            dy.append(np.array(cpp_tr[fib]['Y_vs_W'].value(grid)) - np.array(py_tr[fib]['Y_vs_W'].value(grid)))
        dx, dy = np.array(dx), np.array(dy)
        xr, yr = float(np.sqrt(np.mean(dx**2))), float(np.sqrt(np.mean(dy**2)))
        cam_xr.append(xr); cam_yr.append(yr)
        print(f"  {cam}: xrms={xr:.4f}px yrms={yr:.4f}px")
    if cam_xr:
        mxr, myr = np.mean(cam_xr), np.mean(cam_yr)
        maxxr, maxyr = np.max(cam_xr), np.max(cam_yr)
        print(f"  -- night mean: xrms={mxr:.4f}px yrms={myr:.4f}px  (max: xrms={maxxr:.4f} yrms={maxyr:.4f})")
        all_night_xr.append(mxr); all_night_yr.append(myr)

print(f"\n=== OVERALL ({len(all_night_xr)} nights) ===")
print(f"mean xrms={np.mean(all_night_xr):.4f}px  mean yrms={np.mean(all_night_yr):.4f}px")
print(f"night means -- xrms range [{min(all_night_xr):.4f}, {max(all_night_xr):.4f}], "
      f"yrms range [{min(all_night_yr):.4f}, {max(all_night_yr):.4f}]")
