import os
import sys
import time
import numpy as np
import multiprocessing as mp

from .io import load_python_psf, read_preproc, read_lamp_lines, write_python_psf
from .fitter import PSF_Fitter, get_bundle_spots

# --- Original C++ Wrapper (for baseline and legacy tools) ---

def run_specex(com):
    """
    Original C++ wrapper. This allows desi_psf_fit to run using the C++ core.
    """
    from ._libspecex import (PyOptions, PyIO, PyPrior, PyPSF, PyFitting, VectorString)
    from .io import read_psf, write_psf
    from .qa import specex_psf_qa
    import fitsio

    # instantiate specex C++ objects exposed to python        
    opts = PyOptions() 
    pyio = PyIO()      
    pypr = PyPrior()   
    pyps = PyPSF()     
    pyft = PyFitting() 
    
    spxargs = VectorString()
    for strs in com:
        spxargs.append(strs)

    # parse args
    retval = opts.parse(spxargs)
    if retval != 0: return retval

    # read psf
    read_psf(opts, pyps)

    pyio.set_inputpsf(opts,pyps)
    pypr.set_priors(opts)
    
    # We need read_preproc to return a C++ PyImage for the C++ fitter
    from .io import read_preproc_cpp
    pymg = read_preproc_cpp(opts) 
    
    retval = pyft.fit_psf(opts,pyio,pypr,pymg,pyps) 
    
    # write psf 
    write_psf(pyps,opts,pyio)        

    # do QA
    # retval += specex_psf_qa(opts)

    return retval

# --- New High-Performance Python/JAX Driver ---

def fit_bundle_task(bid, gpu_id, arc_file, in_psf_file, lamp_lines_file):
    """
    Isolated task for fitting a single bundle on a specific GPU.
    """
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
    
    import jax
    import jax.numpy as jnp
    
    class Opts:
        def __init__(self):
            self.arc_image_filename = arc_file
            self.input_psf_filename = in_psf_file
    opts = Opts()
    
    ddata = read_preproc(opts)
    image = ddata['image'].T
    weight = ddata['ivar'].T
    
    psf = load_python_psf(in_psf_file, opts)
    psf.h_size_y = 5
    lamp_lines = read_lamp_lines(lamp_lines_file)
    
    f_min, f_max = bid * 25, (bid + 1) * 25 - 1
    spots = get_bundle_spots(psf, f_min, f_max, lamp_lines, 
                             image=image, weight=weight,
                             min_dist_angstrom=0.0, sn_threshold=3.0)
    
    fitter = PSF_Fitter(psf)
    chi2, pc, tc, cc = fitter.fit(image, weight, spots, bid, max_iter=50)
    
    return bid, {'chi2': chi2, 'psf_coeffs': pc, 'trace_coeffs': tc, 'continuum': cc}

def fit_ccd_native(arc_file, in_psf_file, out_psf_file, lamp_lines_file, first_bundle=0, last_bundle=19):
    """
    Main driver using Python multiprocessing for Multi-GPU scaling.
    """
    print("--- SPECE-X Multi-GPU CCD Fit ---", flush=True)
    t_start = time.time()
    
    all_bundles = list(range(first_bundle, last_bundle + 1))
    n_gpus = 4 
    
    bundle_results = {}
    with mp.Pool(processes=n_gpus) as pool:
        tasks = []
        for i, bid in enumerate(all_bundles):
            gpu_id = i % n_gpus
            tasks.append((bid, gpu_id, arc_file, in_psf_file, lamp_lines_file))
        chunk_results = pool.starmap(fit_bundle_task, tasks)
        for bid, res in chunk_results:
            bundle_results[bid] = res

    t_end = time.time()
    print(f"Total CCD Fit Time: {t_end - t_start:.2f}s", flush=True)
    
    if out_psf_file:
        write_python_psf(out_psf_file, bundle_results, in_psf_file)

if __name__ == "__main__":
    fit_ccd_native(
        arc_file='/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/preproc/20260401/00344649/preproc-z8-00344649.fits.gz',
        in_psf_file='/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/exposures/20260401/00344649/shifted-input-psf-z8-00344649.fits',
        out_psf_file='python-gpu-fit-z8-00344649.fits',
        lamp_lines_file='py/specex/data/specex_linelist_desi.txt',
        first_bundle=0,
        last_bundle=19
    )
