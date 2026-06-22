import os
import sys
import time
import numpy as np
import fitsio

# Setup environment
current_dir = os.getcwd()
sys.path.insert(0, os.path.join(current_dir, 'py'))
os.environ["PYTHONPATH"] = os.path.join(current_dir, 'py') + ":" + os.path.join(current_dir, 'build') + ":" + os.environ.get("PYTHONPATH", "")

from specex.io import read_lamp_lines, load_python_psf, read_image
from specex.fitter import get_bundle_spots, PSF_Fitter
from specex.specex import run_specex

def run_edge_case(name, arc, psf, bundles, broken=None):
    print(f"\n>>> Running Edge Case: {name}")
    print(f"Arc: {arc}")
    print(f"PSF: {psf}")
    
    lamp_lines = "py/specex/data/specex_linelist_desi.txt"
    lines = read_lamp_lines(lamp_lines)
    
    ddata = read_image(arc)
    image = ddata['image'].T; weight = ddata['ivar'].T
    
    class Dummy: pass
    opts = Dummy(); opts.arc_image_filename = arc; opts.input_psf_filename = psf
    psf_obj = load_python_psf(psf, opts)
    
    for b_id in bundles:
        print(f"--- Bundle {b_id} ---")
        fmin, fmax = b_id * 25, (b_id + 1) * 25 - 1
        spots = get_bundle_spots(psf_obj, fmin, fmax, lines, image=image, weight=weight, sn_threshold=3.0, broken_fibers=broken)
        print(f"Spots found: {len(spots)}")
        
        if not spots:
            print("No spots found, skipping fit.")
            continue
            
        fitter = PSF_Fitter(psf_obj)
        t0 = time.time()
        try:
            chi2, pc, tc, cc = fitter.fit(image, weight, spots, b_id, max_iter=15)
            dt = time.time() - t0
            print(f"Fit Successful: Time={dt:.1f}s, Chi2={chi2:.1f}")
        except Exception as e:
            print(f"Fit Failed: {e}")

if __name__ == "__main__":
    # Case 1: Bad Amp (Missing A)
    run_edge_case("Bad Amp (20211028/106399/r8)", 
                  "/global/cfs/cdirs/desi/spectro/redux/daily/preproc/20211028/00106399/preproc-r8-00106399.fits",
                  "/global/cfs/cdirs/desi/spectro/redux/daily/exposures/20211028/00106399/shifted-input-psf-r8-00106399.fits",
                  [0, 10], broken="473,474")

    # Case 2: Overlap
    run_edge_case("Overlap (20250822/307722/z7)",
                  "/global/cfs/cdirs/desi/spectro/redux/daily/preproc/20250822/00307722/preproc-z7-00307722.fits.gz",
                  "/global/cfs/cdirs/desi/spectro/redux/daily/exposures/20250822/00307722/shifted-input-psf-z7-00307722.fits",
                  [10], broken="473,474")

    # Case 3: Known Failure
    run_edge_case("Known Failure (20211028/106396/r8)",
                  "/global/cfs/cdirs/desi/spectro/redux/daily/preproc/20211028/00106396/preproc-r8-00106396.fits",
                  "/global/cfs/cdirs/desi/spectro/redux/daily/exposures/20211028/00106396/stash/shifted-input-psf-r8-00106396.fits",
                  [5, 10], broken="473,474")
