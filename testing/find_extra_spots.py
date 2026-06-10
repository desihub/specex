import numpy as np
from specex.io import read_lamp_lines, load_python_psf, read_image
from specex.fitter import get_bundle_spots

def find_extra_spots():
    arc_file = '/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/preproc/20260401/00344649/preproc-z8-00344649.fits.gz'
    psf_file = '/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/exposures/20260401/00344649/shifted-input-psf-z8-00344649.fits'
    
    class Opts:
        def __init__(self):
            self.arc_image_filename = arc_file
            self.input_psf_filename = psf_file
    opts = Opts()
    
    ddata = read_image(arc_file)
    psf = load_python_psf(psf_file, opts)
    lamp_lines = read_lamp_lines('py/specex/data/specex_linelist_desi.txt')
    
    # Get spots at SN > 0.0 to see all candidates
    all_spots = get_bundle_spots(psf, 125, 149, lamp_lines, image=ddata['image'].T, weight=ddata['ivar'].T, sn_threshold=0.0, broken_fibers="473,474")
    
    spots_30 = [s for s in all_spots if s['snr'] > 3.0]
    spots_315 = [s for s in all_spots if s['snr'] > 3.15]
    
    keys_315 = set([(s['fiber'], s['wave']) for s in spots_315])
    extra = [s for s in spots_30 if (s['fiber'], s['wave']) not in keys_315]
    
    print(f"Total candidates: {len(all_spots)}")
    print(f"Spots at SN>3.0: {len(spots_30)}")
    print(f"Spots at SN>3.15: {len(spots_315)}")
    print(f"\nDifference (the 10 extra spots):")
    for s in extra:
        i0, i1 = max(0, s['stamp_imin']), min(4114, s['stamp_imax'])
        j0, j1 = max(0, s['stamp_jmin']), min(4128, s['stamp_jmax'])
        w_stamp = ddata['ivar'].T[i0:i1, j0:j1]
        ndead = np.sum(w_stamp == 0)
        npix = w_stamp.size
        print(f"  Fiber {s['fiber']}, Wave {s['wave']:.3f}, S/N {s['snr']:.3f}, ndead {ndead}/{npix}, Pos ({s['xc_init']:.1f}, {s['yc_init']:.1f})")

if __name__ == "__main__":
    find_extra_spots()
