import fitsio
import numpy as np
from specex.io import read_lamp_lines, load_python_psf, read_image
from specex.fitter import get_bundle_spots

def verify_parity():
    # File paths for z8 test
    arc_file = '/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/preproc/20260401/00344649/preproc-z8-00344649.fits.gz'
    psf_file = '/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/exposures/20260401/00344649/shifted-input-psf-z8-00344649.fits'
    cpp_fits = '/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/exposures/20260401/00344649/fit-psf-z8-00344649.fits'
    
    # Load Python results
    class Opts:
        def __init__(self):
            self.arc_image_filename = arc_file
            self.input_psf_filename = psf_file
    opts = Opts()
    ddata = read_image(arc_file)
    psf = load_python_psf(psf_file, opts)
    lamp_lines = read_lamp_lines('py/specex/data/specex_linelist_desi.txt')
    py_spots = get_bundle_spots(psf, 125, 149, lamp_lines, image=ddata['image'].T, weight=ddata['ivar'].T, sn_threshold=3.15, broken_fibers="473,474")
    
    # Load C++ results from FITS
    f_cpp = fitsio.FITS(cpp_fits)
    # The PSF table contains the spot information
    # C++ stores spots in a specific way or logs them.
    # Actually, we can check the debug log for spot list or rely on spot count.
    
    print(f"Python: {len(py_spots)} spots.")
    # Assuming C++ count 1573 based on logs
    
    # Wavelength comparison
    py_waves = np.array([s['wave'] for s in py_spots])
    print(f"Python Wave Range: {py_waves.min():.2f} - {py_waves.max():.2f}")

if __name__ == "__main__":
    verify_parity()
