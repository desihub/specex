#!/usr/bin/env python3
"""Per-fiber xrms/yrms breakdown for one camera, C++ vs Python -- to find
whether a camera-level RMS elevation is one/few bad fibers or a systematic
offset across the whole camera.

Usage:
  python testing/per_fiber_breakdown.py --night 20241021 --expid 00259030 \
      --camera z6 --py-dir /pscratch/.../cpptest/redux/cdwarner \
      --cpp-base /pscratch/.../cpptest/desiproc2 --top 20
"""
import os
import sys
import argparse
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from full_ccd_campaign import load_traces
from run_night import find_cases, ALL_CAMERAS

ap = argparse.ArgumentParser()
ap.add_argument("--night", required=True)
ap.add_argument("--expid", required=True)
ap.add_argument("--camera", required=True)
ap.add_argument("--py-dir", required=True, help="dir containing fit-psf-<cam>-<expid>.fits from --backend python")
ap.add_argument("--cpp-base", required=True, help="dir containing <night>_<expid>/redux/<specprod>/exposures/<night>/<expid>/fit-psf-<cam>-<expid>.fits")
ap.add_argument("--specprod", default=os.environ.get("USER", "specex"))
ap.add_argument("--top", type=int, default=15)
args = ap.parse_args()

cases = find_cases(args.night, args.expid, ALL_CAMERAS)
broken = set()
bf = cases.get(args.camera, {}).get("broken_fibers", "")
if bf:
    broken = set(int(x) for x in bf.split(",") if x.strip())
print(f"broken fibers for {args.camera}: {sorted(broken) if broken else 'none'}")

py_fits = os.path.join(args.py_dir, f"fit-psf-{args.camera}-{args.expid}.fits")
cpp_fits = os.path.join(args.cpp_base, f"{args.night}_{args.expid}", "redux", args.specprod,
                         "exposures", args.night, args.expid, f"fit-psf-{args.camera}-{args.expid}.fits")

py_tr, wmin, wmax = load_traces(py_fits)
cpp_tr, _, _ = load_traces(cpp_fits)
grid = np.linspace(wmin, wmax, 100)

rows = []
for fib in range(500):
    if fib in broken:
        continue
    dx = np.array(cpp_tr[fib]['X_vs_W'].value(grid)) - np.array(py_tr[fib]['X_vs_W'].value(grid))
    dy = np.array(cpp_tr[fib]['Y_vs_W'].value(grid)) - np.array(py_tr[fib]['Y_vs_W'].value(grid))
    xr = float(np.sqrt(np.mean(dx**2)))
    yr = float(np.sqrt(np.mean(dy**2)))
    rows.append((fib, xr, yr))

xr_all = np.array([r[1] for r in rows])
yr_all = np.array([r[2] for r in rows])
print(f"\n{args.camera}@{args.night}/{args.expid}: {len(rows)} fibers")
print(f"  overall: xrms={np.sqrt(np.mean(xr_all**2)):.4f}px yrms={np.sqrt(np.mean(yr_all**2)):.4f}px")
print(f"  median per-fiber: xrms={np.median(xr_all):.4f}px yrms={np.median(yr_all):.4f}px")

rows_sorted = sorted(rows, key=lambda r: -r[2])
print(f"\nTop {args.top} fibers by yrms:")
print(f"{'fiber':>6} {'bundle':>7} {'xrms_px':>9} {'yrms_px':>9}")
for fib, xr, yr in rows_sorted[:args.top]:
    print(f"{fib:>6} {fib//25:>7} {xr:>9.4f} {yr:>9.4f}")

rows_sorted_x = sorted(rows, key=lambda r: -r[1])
print(f"\nTop {args.top} fibers by xrms:")
print(f"{'fiber':>6} {'bundle':>7} {'xrms_px':>9} {'yrms_px':>9}")
for fib, xr, yr in rows_sorted_x[:args.top]:
    print(f"{fib:>6} {fib//25:>7} {xr:>9.4f} {yr:>9.4f}")
