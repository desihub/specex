#!/usr/bin/env python3
"""Copy real production preproc + shifted-input-psf files into a private
redux tree at the exact paths desi_proc's own idempotency checks look for
(findfile('preproc', ...) / .../exposures/.../shifted-input-psf-*), so a
`run_night.py --backend cpp` run skips regenerating them and goes straight to
the real -n101 specex fit stage against genuine, complete production input --
avoiding any private-SPECPROD preprocessing/calib-state gaps (like a missing
preproc-<cam> file seen on a from-scratch desi_proc rerun that never
generates its own preproc for every camera).

Usage:
  python testing/stage_preproc.py --night 20260316 --expid 00342128 \
      --redux-dir /pscratch/.../desiproc2/20260316_00342128/redux \
      --specprod cdwarner
  python testing/run_night.py --night 20260316 --expid 00342128 --backend cpp \
      --outdir /pscratch/.../desiproc2/20260316_00342128 \
      --redux-dir /pscratch/.../desiproc2/20260316_00342128/redux --specprod cdwarner
"""
import os
import sys
import shutil
import argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_night import find_cases, ALL_CAMERAS

ap = argparse.ArgumentParser()
ap.add_argument("--night", required=True)
ap.add_argument("--expid", required=True)
ap.add_argument("--redux-dir", required=True, help="private redux root, e.g. <outdir>/redux")
ap.add_argument("--specprod", default=os.environ.get("USER", "specex"))
args = ap.parse_args()

cases = find_cases(args.night, args.expid, ALL_CAMERAS)
missing = set(ALL_CAMERAS) - set(cases)
if missing:
    print(f"ERROR: no arc log entry for {sorted(missing)} -- aborting", file=sys.stderr)
    sys.exit(1)

preproc_dir = os.path.join(args.redux_dir, args.specprod, "preproc", args.night, args.expid)
exposures_dir = os.path.join(args.redux_dir, args.specprod, "exposures", args.night, args.expid)
os.makedirs(preproc_dir, exist_ok=True)
os.makedirs(exposures_dir, exist_ok=True)

n_copied = 0
for cam in ALL_CAMERAS:
    case = cases[cam]
    src_preproc = case["image"]
    src_psf = case["input_psf"]
    dst_preproc = os.path.join(preproc_dir, os.path.basename(src_preproc))
    dst_psf = os.path.join(exposures_dir, os.path.basename(src_psf))
    for src, dst in ((src_preproc, dst_preproc), (src_psf, dst_psf)):
        if not os.path.exists(src):
            print(f"ERROR: {cam}: source file missing: {src}", file=sys.stderr)
            sys.exit(1)
        if not os.path.exists(dst):
            shutil.copy(src, dst)
            n_copied += 1
print(f"staged {len(ALL_CAMERAS)} cameras ({n_copied} files copied, rest already present) "
      f"into {preproc_dir} and {exposures_dir}")
