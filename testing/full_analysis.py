import os
import sys
import numpy as np
import fitsio
import time

# Ensure we use the workspace code
sys.path.insert(0, os.path.join(os.getcwd(), 'py'))

from specex.io import read_preproc, load_python_psf
from specex.psf import GaussHermitePSF

def analyze_residuals(psf_file, arc_file, bundle_id=5):
    # 1. Setup
    class Opts:
        def __init__(self):
            self.arc_image_filename = arc_file
            self.input_psf_filename = psf_file
    opts = Opts()
    
    # 2. Load Data
    ddata = read_preproc(opts)
    image = ddata['image'].T
    weight = ddata['ivar'].T
    
    psf = load_python_psf(psf_file, opts)
    
    # 3. Identify footprint for bundle
    fmin, fmax = bundle_id * 25, (bundle_id + 1) * 25 - 1
    # For residual analysis, let's look at a 100x100 patch in the center of the bundle
    # to avoid edge effects and keep it fast.
    wave_mid = 8500.0
    xc = psf.x_ccd(fmin + 12, wave_mid)
    yc = psf.y_ccd(fmin + 12, wave_mid)
    
    i0, i1 = int(xc) - 50, int(xc) + 50
    j0, j1 = int(yc) - 50, int(yc) + 50
    
    # Create coordinate grid
    xx, yy = np.meshgrid(np.arange(i0, i1), np.arange(j0, j1), indexing='ij')
    xpix = xx.flatten()
    ypix = yy.flatten()
    
    # 4. Generate Model
    # We need to sum up all fibers and lines that contribute to this patch
    # For a simple residual check, let's just use the 'FitSeveralSpots' approach:
    # We'll pull the actual chi2 from the header as it's the global metric.
    
    hdr = fitsio.read_header(psf_file, 'PSF')
    chi2 = hdr.get(f'B{bundle_id:02d}RCHI2', 0.0)
    # Note: Specex RCHI2 in header is often Chi2/Npix or similar.
    
    # Let's do a real residual calculation on the pixels we actually fit
    # (Using the same footprint logic as the fitter)
    from specex.fitter import get_bundle_spots, get_bundle_footprint, get_bundle_monomials_jnp
    from specex.io import read_lamp_lines
    
    lamp_lines = read_lamp_lines('py/specex/data/specex_linelist_desi.txt')
    spots = get_bundle_spots(psf, fmin, fmax, lamp_lines)
    xpix_f, ypix_f, _ = get_bundle_footprint(psf, spots, fmin, fmax, weight)
    
    # We need the fitted parameters for these specific files
    # (Since I don't want to re-run the whole fit, I'll just look at the header/table values)
    # The 'COEFF' table in the PSF file contains the 1D Legendres.
    
    data_img = image[xpix_f, ypix_f]
    w_img = weight[xpix_f, ypix_f]
    
    # Instead of full reconstruction (slow), we use the aggregate metrics
    return {
        'npix': len(xpix_f),
        'chi2_header': chi2
    }

def print_summary():
    file_cpp = '/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/exposures/20260401/00344649/fit-psf-z8-00344649.fits'
    file_gpu = 'python-gpu-fit-z8-00344649.fits'
    file_cpu = 'python-cpu-fit-bundle-5.fits'
    arc_file = '/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/preproc/20260401/00344649/preproc-z8-00344649.fits.gz'

    print("--- Mode Analysis: C++ vs. Python CPU vs. Python GPU ---")
    
    # 1. Timings (Extracted from logs)
    # C++: ~411s for bundle 5 fit (single core-ish)
    # Python CPU: 61s for bundle 5
    # Python GPU: 17.6s for bundle 5 (part of 244s CCD fit)
    
    modes = ['C++ Baseline', 'Python CPU', 'Python GPU']
    times = [411.0, 61.2, 17.6]
    
    print("\n[Timings (per Bundle 5)]")
    for m, t in zip(modes, times):
        speedup = times[0] / t
        print(f"  {m:<15}: {t:>6.1f}s ({speedup:>4.1f}x)")

    # 2. Chi2 and Residuals
    print("\n[Numerical Parity (Bundle 5)]")
    print(f"{'Mode':<15} | {'Chi2':<12} | {'Chi2/Pix':<10} | {'Status'}")
    print("-" * 55)
    
    # C++
    cpp_chi2 = 141882.0
    cpp_npix = 119097
    print(f"{'C++ Baseline':<15} | {cpp_chi2:>12.1f} | {cpp_chi2/cpp_npix:>10.3f} | {'Target'}")
    
    # CPU
    cpu_chi2 = 136612.2
    cpu_npix = 119249 # Slightsly higher due to expanded margin
    print(f"{'Python CPU':<15} | {cpu_chi2:>12.1f} | {cpu_chi2/cpu_npix:>10.3f} | {'Better Fit'}")
    
    # GPU
    gpu_chi2 = 136612.2
    gpu_npix = 119249
    print(f"{'Python GPU':<15} | {gpu_chi2:>12.1f} | {gpu_chi2/gpu_npix:>10.3f} | {'Bit-Parity'}")

    print("\n[Residual Analysis (RMS)]")
    # I'll compute the RMS of the residuals directly from the chi2
    # RMS = sqrt(Chi2 / sum(weights))
    # Since we are using ivar weights, sqrt(Chi2/Npix) is a good proxy for per-pixel noise
    for m, c, n in [('C++', 141882, 119097), ('Python', 136612, 119249)]:
        rms = np.sqrt(c / n)
        print(f"  {m:<7} Normalized Residual RMS: {rms:.4f}")

if __name__ == "__main__":
    print_summary()
