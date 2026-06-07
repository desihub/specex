import os
import sys
import numpy as np

current_dir = os.getcwd()
sys.path.insert(0, os.path.join(current_dir, 'py'))

from specex.io import load_python_psf, read_lamp_lines, read_preproc
from specex.fitter import get_bundle_spots

def find_params():
    arc_file = '/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/preproc/20260401/00344649/preproc-z8-00344649.fits.gz'
    psf_file = '/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/exposures/20260401/00344649/shifted-input-psf-z8-00344649.fits'
    lamp_file = 'py/specex/data/specex_linelist_desi.txt'
    
    class Opts:
        def __init__(self):
            self.arc_image_filename = arc_file
            self.input_psf_filename = psf_file
    opts = Opts()
    ddata = read_preproc(opts)
    weight = ddata['ivar'].T
    psf = load_python_psf(psf_file, opts)
    lamp_lines = read_lamp_lines(lamp_file)
    
    # 1. Spot Selection
    spots = get_bundle_spots(psf, 125, 149, lamp_lines)
    print(f"Number of spots: {len(spots)}")
    
    # 2. Footprint Grid Search
    target_npix = 119097
    
    for margin in [3, 4, 5, 6, 7]:
        for hsx in [5, 6, 7, 8]:
            for hsy in [3, 4, 5, 6]:
                # Temporary override
                psf.h_size_x = hsx
                psf.h_size_y = hsy
                # Re-select spots with new stamp size
                spots_test = get_bundle_spots(psf, 125, 149, lamp_lines)
                
                # Compute NPix
                xmin_env = np.zeros(4128); xmax_env = np.zeros(4128)
                rows_j = np.arange(4128).astype(float)
                w1 = psf.fiber_traces[125]['Y_vs_W'].invert(rows_j)
                w2 = psf.fiber_traces[149]['Y_vs_W'].invert(rows_j)
                x1 = psf.x_ccd(125, w1); x2 = psf.x_ccd(149, w2)
                xmin_env = np.floor(np.minimum(x1, x2) + 0.5).astype(int) - margin
                xmax_env = np.floor(np.maximum(x1, x2) + 0.5).astype(int) + margin + 1
                
                pixels = set()
                for s in spots_test:
                    j_min, j_max = max(0, s['stamp_jmin']), min(4128, s['stamp_jmax'])
                    for j in range(j_min, j_max):
                        i_s = max(s['stamp_imin'], xmin_env[j]); i_e = min(s['stamp_imax'], xmax_env[j])
                        for i in range(i_s, i_e):
                            if 0 <= i < 4114 and weight[i, j] > 0: pixels.add((i, j))
                
                npix = len(pixels)
                if abs(npix - target_npix) < 1000:
                    print(f"MATCH? margin={margin}, hsx={hsx}, hsy={hsy} -> npix={npix}")

if __name__ == "__main__":
    find_params()
