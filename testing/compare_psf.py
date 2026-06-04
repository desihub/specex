import os
import sys
import time
import numpy as np
import jax
import jax.numpy as jnp

# Ensure we use the current workspace code
current_dir = os.getcwd()
sys.path.insert(0, os.path.join(current_dir, 'py'))

from specex._libspecex import (PyOptions, PyIO, PyPrior, PyPSF, PyFitting, VectorString)
from specex.io import read_preproc, read_psf, load_python_psf, read_image
from specex.fitter import PSF_Fitter, read_lamp_lines, get_bundle_spots

def run_comparison():
    # 1. Setup options
    lamp_lines_file = os.path.join(current_dir, 'py/specex/data/specex_linelist_desi.txt')
    
    com = [
        'desi_psf_fit',
        '-a', '/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/preproc/20260401/00344649/preproc-z8-00344649.fits.gz',
        '--in-psf', '/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/exposures/20260401/00344649/shifted-input-psf-z8-00344649.fits',
        '--lamp-lines', lamp_lines_file,
        '--first-bundle', '5',
        '--last-bundle', '5',
        '--first-fiber', '125',
        '--last-fiber', '149',
        '--legendre-deg-wave', '3',
        '--fit-continuum'
    ]

    opts = PyOptions()
    spxargs = VectorString()
    for s in com: spxargs.append(s)
    opts.parse(spxargs)

    # 2. Load data
    print("Loading data...")
    pyps_baseline = PyPSF()
    read_psf(opts, pyps_baseline)
    pyio = PyIO()
    pyio.load_psf(opts, pyps_baseline)
    # read_preproc returns spx.PyImage for C++ compatibility
    pymg_cpp = read_preproc(opts)
    
    # Also load for Python fit
    ddata = read_image(opts.arc_image_filename)
    image = ddata['image'].T
    weight = ddata['ivar'].T
    mask = ddata['mask'].T
    weight[mask != 0] = 0.0

    # 3. Running C++ Fit (Baseline)
    print("\n--- SKIPPING C++ Fit ---")
    pyio = PyIO()
    pypr = PyPrior()
    pyft = PyFitting()
    
    t0 = time.time()
    retval_cpp = pyft.fit_psf(opts, pyio, pypr, pymg_cpp, pyps_baseline)
    t1 = time.time()
    print(f"C++ Time: {t1 - t0:.2f}s, Return: {retval_cpp}")

    # 4. Running Python/JAX Fit
    print("\n--- Running Python/JAX Fit ---")
    print(f"JAX devices: {jax.devices()}")
    
    # psf_py is our pure-python PSF object
    psf_py = load_python_psf(opts.input_psf_filename, opts)
    lamp_lines = read_lamp_lines(lamp_lines_file)
    print(f"Loaded {len(lamp_lines)} lamp lines.")
    
    bundle_id = 5
    # C++ Final Stage: Include blended lines (dist=0) and noisier spots (S/N > 3)
    spots = get_bundle_spots(psf_py, 125, 149, lamp_lines, 
                             image=image, weight=weight,
                             min_dist_angstrom=0.0, sn_threshold=3.0,
                             wave_min=psf_py.fiber_traces[125]['X_vs_W'].xmin, 
                             wave_max=psf_py.fiber_traces[125]['X_vs_W'].xmax)
    print(f"Reconstructed {len(spots)} spots (after filtering).")
    
    fitter = PSF_Fitter(psf_py)
    
    t0 = time.time()
    # Run full non-linear fit for final convergence
    final_chi2 = fitter.fit(image, weight, spots, bundle_id, fit_type='full', max_iter=50)
    t1 = time.time()
    print(f"Python/JAX Time: {t1 - t0:.2f}s")
    print(f"Final Python Chi2: {final_chi2:.4f}")

if __name__ == "__main__":
    run_comparison()
