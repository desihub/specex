import os
import sys
import time
import numpy as np

# Force JAX to CPU
os.environ["JAX_PLATFORM_NAME"] = "cpu"

from specex.io import read_preproc, load_python_psf, read_lamp_lines, write_python_psf
from specex.fitter import PSF_Fitter, get_bundle_spots

def run_cpu_one_bundle(bid=5):
    arc_file = '/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/preproc/20260401/00344649/preproc-z8-00344649.fits.gz'
    psf_file = '/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/exposures/20260401/00344649/shifted-input-psf-z8-00344649.fits'
    lamp_file = 'py/specex/data/specex_linelist_desi.txt'
    
    class Opts:
        def __init__(self):
            self.arc_image_filename = arc_file
            self.input_psf_filename = psf_file
    opts = Opts()
    
    ddata = read_preproc(opts)
    image = ddata['image'].T; weight = ddata['ivar'].T
    psf = load_python_psf(psf_file, opts)
    lamp_lines = read_lamp_lines(lamp_file)
    
    fmin, fmax = bid*25, (bid+1)*25 - 1
    spots = get_bundle_spots(psf, fmin, fmax, lamp_lines, image=image, weight=weight, sn_threshold=3.0)
    
    fitter = PSF_Fitter(psf)
    print(f"Running CPU Fit for bundle {bid}...")
    t0 = time.time()
    chi2, pc, tc, cont = fitter.fit(image, weight, spots, bid, max_iter=50)
    t1 = time.time()
    print(f"CPU Fit done in {t1-t0:.2f}s. Chi2: {chi2:.1f}")
    
    results = {bid: {'chi2': chi2, 'psf_coeffs': pc, 'trace_coeffs': tc, 'continuum': cont, 'spots': spots}}
    write_python_psf(f'python-cpu-fit-bundle-{bid}.fits', results, psf_file)

if __name__ == "__main__":
    sys.path.insert(0, os.path.join(os.getcwd(), 'py'))
    run_cpu_one_bundle(5)
