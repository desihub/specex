import os
import sys
import time
import numpy as np
import multiprocessing as mp

# Note: We do NOT import JAX here at the top level to ensure 
# CUDA_VISIBLE_DEVICES can be set cleanly in child processes.

from .io import load_python_psf, read_preproc, read_lamp_lines
from .fitter import PSF_Fitter, get_bundle_spots

def fit_bundle_task(bid, gpu_id, arc_file, in_psf_file, lamp_lines_file):
    """
    Isolated task for fitting a single bundle on a specific GPU.
    """
    # 1. Set GPU affinity before JAX initializes
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
    
    import jax
    import jax.numpy as jnp
    
    # 2. Setup (Each process needs its own local state)
    class Opts:
        def __init__(self):
            self.arc_image_filename = arc_file
            self.input_psf_filename = in_psf_file
    opts = Opts()
    
    ddata = read_preproc(opts)
    image = ddata['image'].T
    weight = ddata['ivar'].T
    
    psf = load_python_psf(in_psf_file, opts)
    lamp_lines = read_lamp_lines(lamp_lines_file)
    
    f_min, f_max = bid * 25, (bid + 1) * 25 - 1
    
    print(f"[Bundle {bid}] Starting fit on GPU {gpu_id}...", flush=True)
    
    spots = get_bundle_spots(psf, f_min, f_max, lamp_lines, 
                             image=image, weight=weight,
                             min_dist_angstrom=0.0, sn_threshold=3.0,
                             wave_min=psf.fiber_traces[f_min]['X_vs_W'].xmin,
                             wave_max=psf.fiber_traces[f_min]['X_vs_W'].xmax)
    
    fitter = PSF_Fitter(psf)
    # Run the high-performance JAX-GPU fit
    chi2, pc, tc, cont = fitter.fit(image, weight, spots, bid, fit_type='full', max_iter=50)
    
    return bid, {
        'chi2': chi2,
        'psf_coeffs': pc,
        'trace_coeffs': tc,
        'continuum': cont,
        'spots': spots # Useful for QA
    }

def fit_ccd_native(arc_file, in_psf_file, out_psf_file, lamp_lines_file, first_bundle=0, last_bundle=19):
    """
    Main driver using Python multiprocessing for Multi-GPU scaling.
    """
    print("--- SPECE-X Multi-GPU CCD Fit ---", flush=True)
    print(f"Arc Image: {arc_file}", flush=True)
    print(f"Input PSF: {in_psf_file}", flush=True)
    
    t_start = time.time()
    
    all_bundles = list(range(first_bundle, last_bundle + 1))
    n_gpus = 4 # A100 node
    
    # Use a pool to manage 4 concurrent fits
    bundle_results = {}
    with mp.Pool(processes=n_gpus) as pool:
        # Prepare arguments
        tasks = []
        for i, bid in enumerate(all_bundles):
            gpu_id = i % n_gpus
            tasks.append((bid, gpu_id, arc_file, in_psf_file, lamp_lines_file))
            
        # Run parallel
        chunk_results = pool.starmap(fit_bundle_task, tasks)
        
        for bid, res in chunk_results:
            bundle_results[bid] = res

    t_end = time.time()
    
    print("\n--- CCD Fit Summary ---", flush=True)
    for b in sorted(bundle_results.keys()):
        print(f"  Bundle {b:02d}: Chi2 = {bundle_results[b]['chi2']:.1f}", flush=True)
    
    print(f"Total CCD Fit Time: {t_end - t_start:.2f}s", flush=True)
    
    # Save Results
    if out_psf_file:
        from .io import write_python_psf
        print(f"Saving results to {out_psf_file}...", flush=True)
        write_python_psf(out_psf_file, bundle_results, in_psf_file)

if __name__ == "__main__":
    # Full 20-bundle CCD Stress Test
    fit_ccd_native(
        arc_file='/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/preproc/20260401/00344649/preproc-z8-00344649.fits.gz',
        in_psf_file='/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/exposures/20260401/00344649/shifted-input-psf-z8-00344649.fits',
        out_psf_file='python-gpu-fit-z8-00344649.fits',
        lamp_lines_file='py/specex/data/specex_linelist_desi.txt',
        first_bundle=0,
        last_bundle=19
    )
