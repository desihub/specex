import sys
import numpy as np
import fitsio

def compare_psfs(file_a, file_b):
    print(f"Comparing:\n  A: {file_a}\n  B: {file_b}")
    
    f_a = fitsio.FITS(file_a)
    f_b = fitsio.FITS(file_b)
    
    # 1. Compare Traces
    xt_a = f_a['XTRACE'].read(); xt_b = f_b['XTRACE'].read()
    yt_a = f_a['YTRACE'].read(); yt_b = f_b['YTRACE'].read()
    
    dx = xt_a - xt_b
    dy = yt_a - yt_b
    
    print("\n--- Trace Comparison ---")
    print(f"  X-Trace: Mean Diff = {np.mean(dx):.6f}, RMS = {np.std(dx):.6f}, Max = {np.max(np.abs(dx)):.6f}")
    print(f"  Y-Trace: Mean Diff = {np.mean(dy):.6f}, RMS = {np.std(dy):.6f}, Max = {np.max(np.abs(dy)):.6f}")
    
    # 2. Compare PSF Parameters
    p_a = f_a['PSF'].read(); p_b = f_b['PSF'].read()
    names_a = [n.strip() for n in p_a['PARAM']]
    names_b = [n.strip() for n in p_b['PARAM']]
    
    common_names = sorted(list(set(names_a) & set(names_b)))
    print("\n--- PSF Parameter Comparison ---")
    print(f"{'Parameter':<12} | {'Mean Diff':<12} | {'RMS':<12} | {'Max Abs':<12}")
    print("-" * 55)
    
    for name in common_names:
        idx_a = names_a.index(name); idx_b = names_b.index(name)
        c_a = p_a['COEFF'][idx_a]; c_b = p_b['COEFF'][idx_b]
        
        # Specex stores coeffs as (500, 4) or similar.
        # We compare the actual values (not the Legendre coeffs directly, but let's start with coeffs)
        diff = c_a - c_b
        md = np.mean(diff); rms = np.std(diff); mx = np.max(np.abs(diff))
        
        # Only print major params for brevity
        if name in ['GHSIGX', 'GHSIGY', 'GH-0-0', 'TAILAMP', 'CONT']:
            print(f"{name:<12} | {md:>12.6e} | {rms:>12.6e} | {mx:>12.6e}")

if __name__ == "__main__":
    file_gpu = 'python-gpu-fit-z8-00344649.fits'
    file_cpp = '/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/exposures/20260401/00344649/fit-psf-z8-00344649.fits'
    compare_psfs(file_gpu, file_cpp)
