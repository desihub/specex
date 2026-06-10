import os
import sys
import time
import argparse
import numpy as np
import jax
import jax.numpy as jnp

# Ensure we use the current workspace code
current_dir = os.getcwd()
sys.path.insert(0, os.path.join(current_dir, 'py'))

from specex.io import read_preproc, load_python_psf, read_image, read_lamp_lines, write_python_psf
from specex.fitter import PSF_Fitter, get_bundle_spots

class DummyOptions:
    def __init__(self, arc, in_psf):
        self.arc_image_filename = arc
        self.input_psf_filename = in_psf

def run_comparison():
    parser = argparse.ArgumentParser()
    parser.add_argument('--camera', type=str, default='z8')
    parser.add_argument('--bundle', type=int, default=5)
    parser.add_argument('--out-psf', type=str, default='test_py.fits')
    args = parser.parse_args()

    # Construct paths
    night = '20260401'
    expid = '00344649'
    arc_base = '/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/preproc'
    psf_base = '/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/exposures'
    
    arc_file = os.path.join(arc_base, night, expid, f'preproc-{args.camera}-{expid}.fits.gz')
    psf_file = os.path.join(psf_base, night, expid, f'shifted-input-psf-{args.camera}-{expid}.fits')
    lamp_lines_file = os.path.join(current_dir, 'py/specex/data/specex_linelist_desi.txt')
    
    opts = DummyOptions(arc_file, psf_file)

    # Load data
    print(f"Loading data for {args.camera} bundle {args.bundle}...")
    ddata = read_preproc(opts)
    image = ddata['image'].T
    weight = ddata['ivar'].T
    
    # Setup PSF
    psf_py = load_python_psf(opts.input_psf_filename, opts)
    psf_py.h_size_x = 8
    psf_py.h_size_y = 5
    lamp_lines = read_lamp_lines(lamp_lines_file)
    
    # Fiber Range
    fmin, fmax = args.bundle * 25, (args.bundle + 1) * 25 - 1
    # Infer broken fibers (simple heuristic: z-band has broken 367 or 473,474)
    broken = "473,474" if 'z' in args.camera else None
    
    # Run Spot Selection
    spots = get_bundle_spots(psf_py, fmin, fmax, lamp_lines, 
                             image=image, weight=weight,
                             min_dist_angstrom=0.0, sn_threshold=3.15,
                             broken_fibers=broken)
    print(f"Reconstructed {len(spots)} spots (after filtering).")
    
    # Run Fit
    fitter = PSF_Fitter(psf_py)
    fitter.fit(image, weight, spots, args.bundle, max_iter=2) # Warm-up

    t0 = time.time()
    final_chi2, pc, tc, cont = fitter.fit(image, weight, spots, args.bundle, fit_type='full', max_iter=50)
    t1 = time.time()
    
    print(f"\nFinal Results for {args.camera} bundle {args.bundle}:")
    print(f"  Time: {t1 - t0:.2f}s")
    print(f"  Chi2: {final_chi2:.4f}")
    print(f"  Spots: {len(spots)}")

    results = {args.bundle: {'chi2': final_chi2, 'psf_coeffs': pc, 'trace_coeffs': tc, 'continuum': cont}}
    out_file = f'test_py_{args.camera}.fits' if args.out_psf == 'test_py.fits' else args.out_psf
    write_python_psf(out_file, results, psf_file)
    print(f"Saved Python result to {out_file}")

if __name__ == "__main__":
    run_comparison()
