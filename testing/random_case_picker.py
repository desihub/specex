"""
Pick N random (camera, night) test cases per band by scanning random nights'
arc*.log files for desi_compute_psf invocations, instead of restricting to a
single fixed night/expid like select_test_case.py/bundle_parity_suite.py do.

Writes JSON-lines output consumable by full_ccd_campaign.py's --cases-file.

Usage:
  python testing/random_case_picker.py --band b --n 5 --seed 42 \
      --exclude-file used_cases.jsonl --out b_cases.jsonl
"""
import os
import sys
import glob
import json
import random
import argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from select_test_case import parse_log_line

SCRIPTS_DIR = "/global/cfs/cdirs/desi/spectro/redux/matterhorn/run/scripts/night"


def list_nights():
    return sorted(d for d in os.listdir(SCRIPTS_DIR) if d.isdigit() and len(d) == 8)


def cases_in_night(night, band):
    log_dir = os.path.join(SCRIPTS_DIR, night)
    out = []
    for log in glob.glob(os.path.join(log_dir, "arc*.log")):
        try:
            with open(log) as f:
                for line in f:
                    if "desi_compute_psf" not in line:
                        continue
                    case = parse_log_line(line)
                    if case and case.get('camera', '')[:1] == band:
                        case['night'] = night
                        out.append(case)
        except OSError:
            continue
    return out


def pick(band, n, seed, exclude, tried_nights_limit=400):
    rng = random.Random(seed)
    nights = list_nights()
    rng.shuffle(nights)
    picked = []
    seen = set(exclude)
    tried = 0
    for night in nights:
        if len(picked) >= n or tried >= tried_nights_limit:
            break
        tried += 1
        cands = cases_in_night(night, band)
        rng.shuffle(cands)
        for c in cands:
            key = (c['camera'], c['night'], c.get('expid'))
            if key in seen:
                continue
            seen.add(key)
            picked.append(c)
            break  # at most one case per night per band for diversity
    return picked


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--band", required=True, choices=["b", "r", "z"])
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--exclude-file", help="JSON-lines file(s) of previously-used cases to avoid duplicates (comma-sep)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    exclude = set()
    if args.exclude_file:
        for fp in args.exclude_file.split(","):
            if os.path.exists(fp):
                with open(fp) as f:
                    for line in f:
                        if line.strip():
                            c = json.loads(line)
                            exclude.add((c['camera'], c['night'], c.get('expid')))

    picked = pick(args.band, args.n, args.seed, exclude)
    with open(args.out, 'w') as f:
        for c in picked:
            f.write(json.dumps(c) + "\n")
    print(f"Picked {len(picked)}/{args.n} for band {args.band}: "
          f"{[(c['camera'], c['night']) for c in picked]}")


if __name__ == "__main__":
    main()
