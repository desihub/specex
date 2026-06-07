import os
import sys
import time
import numpy as np
import jax
import jax.numpy as jnp

# Ensure we use the current workspace code
current_dir = os.getcwd()
sys.path.insert(0, os.path.join(current_dir, 'py'))

# from specex._libspecex import (PyOptions, PyIO, PyPrior, PyPSF, PyFitting, VectorString)
from specex.io import read_preproc, load_python_psf, read_image, read_lamp_lines
from specex.fitter import PSF_Fitter, get_bundle_spots

class DummyOptions:
    def __init__(self, arc, in_psf):
        self.arc_image_filename = arc
        self.input_psf_filename = in_psf

def run_comparison():
    # 1. Setup options
    lamp_lines_file = os.path.join(current_dir, 'py/specex/data/specex_linelist_desi.txt')
    arc_file = '/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/preproc/20260401/00344649/preproc-z8-00344649.fits.gz'
    psf_file = '/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/exposures/20260401/00344649/shifted-input-psf-z8-00344649.fits'
    
    opts = DummyOptions(arc_file, psf_file)

    # 2. Load data
    print("Loading data...")
    ddata = read_preproc(opts)
    image = ddata['image'].T
    weight = ddata['ivar'].T
    
    # 3. Running C++ Fit (Baseline) - SKIPPED on compute node to avoid libfabric issues
    print("\n--- SKIPPING C++ Fit (Dependency-free mode) ---")

    # 4. Running Python/JAX Fit
    print("\n--- Running Python/JAX Fit ---")
    print(f"JAX devices: {jax.devices()}")
    
    # psf_py is our pure-python PSF object
    psf_py = load_python_psf(opts.input_psf_filename, opts)
    psf_py.h_size_x = 7 # Trial: match C++ 15-pix footprint
    psf_py.h_size_y = 5 # Production baseline uses HSIZEY=5
    lamp_lines = read_lamp_lines(lamp_lines_file)
    print(f"Loaded {len(lamp_lines)} lamp lines.")
    
    bundle_id = 5
    # C++ Final Stage: Include blended lines (dist=0) and noisier spots (S/N > 3)
    spots = get_bundle_spots(psf_py, 125, 149, lamp_lines, 
                             image=image, weight=weight,
                             min_dist_angstrom=0.0, sn_threshold=3.0)
    print(f"Reconstructed {len(spots)} spots (after filtering).")
    
    fitter = PSF_Fitter(psf_py)
    
    t0 = time.time()
    # Run full non-linear fit for final convergence
    final_chi2, pc, tc, cont = fitter.fit(image, weight, spots, bundle_id, fit_type='full', max_iter=50)
    t1 = time.time()
    print(f"Python/JAX Time: {t1 - t0:.2f}s")
    print(f"Final Python Chi2: {final_chi2:.4f}")

if __name__ == "__main__":
    run_comparison()
