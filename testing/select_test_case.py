import os
import sys
import glob
import re
import argparse
import random

def parse_log_line(line):
    """
    Extracts parameters from a 'desi_compute_psf' log line.
    """
    if "desi_compute_psf" not in line:
        return None
    
    res = {}
    # Extract image path to get camera and expid
    img_match = re.search(r'--input-image\s+(\S+)', line)
    if img_match:
        res['image'] = img_match.group(1)
        # Extract camera (e.g. z8, r3)
        cam_match = re.search(r'preproc-([a-z]\d)-', res['image'])
        if cam_match:
            res['camera'] = cam_match.group(1)
        # Extract expid from the filename itself (preproc-{cam}-{expid}.fits.gz)
        # rather than the path -- the path's FIRST 8-digit segment is the
        # night, which also happens to be 8 digits and would be matched
        # incorrectly by a naive '/(\d{8})/' path search.
        exp_match = re.search(r'preproc-[a-z]\d+-(\d+)\.fits', res['image'])
        if exp_match:
            res['expid'] = exp_match.group(1)

    psf_match = re.search(r'--input-psf\s+(\S+)', line)
    if psf_match:
        res['input_psf'] = psf_match.group(1)
        
    broken_match = re.search(r'--broken-fibers\s+(\S+)', line)
    if broken_match:
        res['broken_fibers'] = broken_match.group(1)
    else:
        res['broken_fibers'] = ""
        
    return res

def main():
    parser = argparse.ArgumentParser(description="Scrape DESI logs for PSF fit test cases.")
    parser.add_argument("--night", type=str, default="20260401", help="YYYYMMDD")
    parser.add_argument("--expid", type=str, help="8-digit exposure ID")
    parser.add_argument("--camera", type=str, help="e.g. z8, r3")
    parser.add_argument("--random", action="store_true", help="Select one case at random")
    parser.add_argument("--list", action="store_true", help="List all found cases")
    
    args = parser.parse_args()
    
    log_dir = f"/global/cfs/cdirs/desi/spectro/redux/matterhorn/run/scripts/night/{args.night}"
    log_pattern = os.path.join(log_dir, "arc*.log")
    
    logs = glob.glob(log_pattern)
    if not logs:
        print(f"No logs found in {log_dir}")
        return

    cases = []
    for log in logs:
        with open(log, 'r') as f:
            for line in f:
                if "desi_compute_psf" in line:
                    case = parse_log_line(line)
                    if case:
                        # Filters
                        if args.expid and args.expid not in line: continue
                        if args.camera and f"-{args.camera}-" not in line: continue
                        cases.append(case)
    
    # De-duplicate
    unique_cases = { (c['camera'], c['image']): c for c in cases }
    cases = list(unique_cases.values())

    if args.list:
        print(f"{'Camera':<8} | {'Expid':<10} | {'Broken Fibers':<20}")
        print("-" * 45)
        for c in sorted(cases, key=lambda x: (x.get('camera',''), x.get('expid',''))):
            print(f"{c.get('camera','?'):<8} | {c.get('expid','?'):<10} | {c.get('broken_fibers','none'):<20}")

    if args.random and cases:
        c = random.choice(cases)
        print("\n[Selected Random Case]")
        for k, v in c.items():
            print(f"  {k}: {v}")
    elif not args.list:
        print(f"Found {len(cases)} cases for night {args.night}.")

if __name__ == "__main__":
    main()
