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

def fit_bundle_task(bid, gpu_id, arc_file, in_psf_file, lamp_lines_file, broken_fibers=None, sn_threshold=3.0, h_size_y=None):
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
    if h_size_y is not None:
        psf.h_size_y = h_size_y
        
    lamp_lines = read_lamp_lines(lamp_lines_file)
    
    f_min, f_max = bid * 25, (bid + 1) * 25 - 1
    spots = get_bundle_spots(psf, f_min, f_max, lamp_lines, 
                             image=image, weight=weight,
                             min_dist_angstrom=0.0, sn_threshold=sn_threshold,
                             broken_fibers=broken_fibers)
    
    fitter = PSF_Fitter(psf)
    chi2, pc, tc, cc = fitter.fit(image, weight, spots, bid, max_iter=50)
    
    return bid, {'chi2': chi2, 'psf_coeffs': pc, 'trace_coeffs': tc, 'continuum': cc}

def fit_ccd_native(arc_file, in_psf_file, out_psf_file, lamp_lines_file, 
                   first_bundle=0, last_bundle=19, n_gpus=4, 
                   broken_fibers=None, sn_threshold=3.0, h_size_y=5):
    """
    Main driver using Python multiprocessing for Multi-GPU scaling.
    """
    print("--- SPECE-X Multi-GPU CCD Fit ---", flush=True)
    print(f"  Arc: {arc_file}")
    print(f"  In PSF: {in_psf_file}")
    print(f"  Out PSF: {out_psf_file}")
    print(f"  Broken Fibers: {broken_fibers}")
    t_start = time.time()
    
    all_bundles = list(range(first_bundle, last_bundle + 1))
    
    bundle_results = {}
    # Use 'spawn' for clean GPU initialization in child processes
    ctx = mp.get_context('spawn')
    with ctx.Pool(processes=n_gpus) as pool:
        tasks = []
        for i, bid in enumerate(all_bundles):
            gpu_id = i % n_gpus
            tasks.append((bid, gpu_id, arc_file, in_psf_file, lamp_lines_file, broken_fibers, sn_threshold, h_size_y))
        chunk_results = pool.starmap(fit_bundle_task, tasks)
        for bid, res in chunk_results:
            bundle_results[bid] = res

    t_end = time.time()
    print(f"Total CCD Fit Time: {t_end - t_start:.2f}s", flush=True)
    
    if out_psf_file:
        write_python_psf(out_psf_file, bundle_results, in_psf_file)

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Specex Python/JAX PSF Fitter")
    parser.add_argument("-a", "--arc", "--input-image", type=str, required=True, help="Input preproc arc image")
    parser.add_argument("--in-psf", "--input-psf", type=str, required=True, help="Input (shifted) PSF file")
    parser.add_argument("--out-psf", "--output-psf", type=str, required=True, help="Output PSF file")
    parser.add_argument("--lamp-lines", type=str, help="Lamp lines file")
    parser.add_argument("--first-bundle", type=int, default=0)
    parser.add_argument("--last-bundle", type=int, default=19)
    parser.add_argument("--first-fiber", type=int, help="First fiber to fit (used to derive bundle)")
    parser.add_argument("--last-fiber", type=int, help="Last fiber to fit")
    parser.add_argument("--legendre-deg-wave", type=int, default=3, help="Legendre degree for wavelength")
    parser.add_argument("--fit-continuum", action="store_true", default=True, help="Enable continuum fitting (always on in coupled mode)")
    parser.add_argument("--gpu", type=int, default=4, help="Number of GPUs to use")
    parser.add_argument("--backend", type=str, default="gpu", choices=["cpu", "gpu"])
    parser.add_argument("--broken-fibers", type=str, help="Comma-separated list of broken fibers")
    parser.add_argument("--sn-threshold", type=float, default=3.0, help="S/N threshold for spot selection")
    parser.add_argument("--h-size-y", type=int, default=5, help="Override PSF stamp half-size in Y")
    
    args = parser.parse_args()
    
    if not args.lamp_lines:
        # Try to find default lamp lines
        base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        args.lamp_lines = os.path.join(base, "specex/data/specex_linelist_desi.txt")

    os.environ["JAX_PLATFORM_NAME"] = args.backend
    
    fit_ccd_native(
        arc_file=args.arc,
        in_psf_file=args.in_psf,
        out_psf_file=args.out_psf,
        lamp_lines_file=args.lamp_lines,
        first_bundle=args.first_bundle,
        last_bundle=args.last_bundle,
        n_gpus=args.gpu,
        broken_fibers=args.broken_fibers,
        sn_threshold=args.sn_threshold,
        h_size_y=args.h_size_y
    )

if __name__ == "__main__":
    main()
