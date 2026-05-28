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
from specex.io import read_preproc, read_psf, load_python_psf
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
    pyps = PyPSF()
    read_psf(opts, pyps)
    pymg = read_preproc(opts)
    
    from specex.io import read_image
    dsmg = read_image(opts.arc_image_filename)
    image = dsmg.pix.T
    weight = dsmg.ivar.T

    # 3. Running C++ Fit (SKIP for now)
    print("\n--- Skipping C++ Fit ---")
    # pyio = PyIO()
    # pypr = PyPrior()
    # pyft = PyFitting()
    # t0 = time.time()
    # retval_cpp = pyft.fit_psf(opts, pyio, pypr, dsmg_cpp, pyps) # wait, dsmg_cpp is not defined
    # ...
    
    # 4. Running Python/JAX Fit
    print("\n--- Running Python/JAX Fit ---")
    # Fresh PSF setup
    pyps_fresh = PyPSF()
    read_psf(opts, pyps_fresh) # just headers
    
    psf_py = load_python_psf(pyps_fresh, opts) # full load
    print(f"PSF properties: hSizeX={psf_py.h_size_x} hSizeY={psf_py.h_size_y} fiber_min={psf_py.fiber_min}")
    print(f"PSF trace range: wave={pyps_fresh.trace_WAVEMIN:.1f}-{pyps_fresh.trace_WAVEMAX:.1f}")
    
    lamp_lines = read_lamp_lines(lamp_lines_file)
    print(f"Loaded {len(lamp_lines)} lamp lines.")
    
    bundle_id = 5
    spots = get_bundle_spots(psf_py, 125, 149, lamp_lines, 
                             wave_min=pyps_fresh.trace_WAVEMIN, 
                             wave_max=pyps_fresh.trace_WAVEMAX)
    print(f"Reconstructed {len(spots)} spots (after filtering).")
    if len(spots) > 0:
        print(f"First spot: fiber={spots[0]['fiber']} wave={spots[0]['wave']:.2f} xc={spots[0]['xc']:.2f} yc={spots[0]['yc']:.2f}")
        print(f"First spot stamp: i={spots[0]['stamp_imin']}-{spots[0]['stamp_imax']} j={spots[0]['stamp_jmin']}-{spots[0]['stamp_jmax']}")
    
    fitter = PSF_Fitter(psf_py)
    
    t0 = time.time()
    # Initial flux fit
    fitter.fit(image, weight, spots, bundle_id, fit_type='flux')
    # Full fit
    final_chi2 = fitter.fit(image, weight, spots, bundle_id, fit_type='full')
    t1 = time.time()
    print(f"Python/JAX Time: {t1 - t0:.2f}s")
    print(f"Final Python Chi2: {final_chi2:.4f}")

if __name__ == "__main__":
    run_comparison()
