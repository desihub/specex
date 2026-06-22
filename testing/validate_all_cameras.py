import os
import sys
import time
import argparse
import numpy as np
import jax
import jax.numpy as jnp
import fitsio

# Ensure we use the current workspace code
current_dir = os.getcwd()
sys.path.insert(0, os.path.join(current_dir, 'py'))

from specex.io import read_preproc, load_python_psf, read_lamp_lines, write_python_psf
from specex.fitter import PSF_Fitter, get_bundle_spots

def run_test(arc_file, psf_file, broken_fibers, camera, bundle_id=5):
    """
    Runs a single bundle fit and returns (chi2, time).
    """
    print(f"\n--- Testing Camera {camera} Bundle {bundle_id} ---")
    
    # 1. Load Data
    class Opts:
        def __init__(self, arc, in_psf):
            self.arc_image_filename = arc
            self.input_psf_filename = in_psf
    opts = Opts(arc_file, psf_file)
    
    ddata = read_preproc(opts)
    image = ddata['image'].T
    weight = ddata['ivar'].T
    
    psf_py = load_python_psf(psf_file, opts)
    psf_py.h_size_y = 5 # Standard production size
    
    lamp_lines_file = os.path.join(current_dir, 'py/specex/data/specex_linelist_desi.txt')
    lamp_lines = read_lamp_lines(lamp_lines_file)
    
    # 2. Spot selection
    fmin, fmax = bundle_id * 25, (bundle_id + 1) * 25 - 1
    spots = get_bundle_spots(psf_py, fmin, fmax, lamp_lines, 
                             image=image, weight=weight,
                             min_dist_angstrom=0.0, sn_threshold=3.0)
    
    print(f"Reconstructed {len(spots)} spots.")
    
    fitter = PSF_Fitter(psf_py)
    
    # 3. Fit
    t0 = time.time()
    chi2, pc, tc, cc = fitter.fit(image, weight, spots, bundle_id, max_iter=50)
    t1 = time.time()
    
    elapsed = t1 - t0
    print(f"Done in {elapsed:.2f}s. Chi2: {chi2:.1f}")
    
    return chi2, elapsed

def main():
    parser = argparse.ArgumentParser(description="Batch validator for multiple cameras.")
    parser.add_argument("--night", type=str, default="20260401")
    parser.add_argument("--expid", type=str, default="00344649")
    parser.add_argument("--cameras", type=str, help="Comma-separated list of cameras (e.g. b0,r3,z8)")
    parser.add_argument("--bundle", type=int, default=5)
    parser.add_argument("--output", type=str, default="validation_results.txt")
    
    args = parser.parse_args()
    
    # Use our scraper to get the real paths and broken fibers
    from select_test_case import parse_log_line
    log_dir = f"/global/cfs/cdirs/desi/spectro/redux/matterhorn/run/scripts/night/{args.night}"
    log_pattern = os.path.join(log_dir, "arc*.log")
    import glob
    logs = glob.glob(log_pattern)
    
    case_map = {}
    for log in logs:
        with open(log, 'r') as f:
            for line in f:
                if "desi_compute_psf" in line and args.expid in line:
                    case = parse_log_line(line)
                    if case: case_map[case['camera']] = case
                    
    selected_cams = args.cameras.split(',') if args.cameras else sorted(case_map.keys())
    
    with open(args.output, "w") as f:
        f.write(f"{'Night':<10} | {'Expid':<10} | {'Cam':<5} | {'Spots':<5} | {'Time(s)':<8} | {'Chi2':<12}\n")
        f.write("-" * 65 + "\n")
        
        for cam in selected_cams:
            if cam not in case_map:
                print(f"Skipping {cam}, no log info found.")
                continue
            
            case = case_map[cam]
            try:
                # Re-run spot count to be precise
                # (Actual run inside run_test)
                chi2, dt = run_test(case['image'], case['input_psf'], case['broken_fibers'], cam, args.bundle)
                
                f.write(f"{args.night:<10} | {args.expid:<10} | {cam:<5} | {'?':<5} | {dt:>8.2f} | {chi2:>12.1f}\n")
                f.flush()
            except Exception as e:
                print(f"Error fitting camera {cam}: {e}")
                f.write(f"{args.night:<10} | {args.expid:<10} | {cam:<5} | ERROR | {'-':>8} | {'-':>12}\n")

if __name__ == "__main__":
    main()
