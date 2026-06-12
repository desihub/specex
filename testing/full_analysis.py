import os
import sys
import time
import re
import numpy as np
import fitsio

# Ensure we use the current workspace code
current_dir = os.getcwd()
sys.path.insert(0, os.path.join(current_dir, 'py'))
os.environ["PYTHONPATH"] = os.path.join(current_dir, 'py') + ":" + os.path.join(current_dir, 'build') + ":" + os.environ.get("PYTHONPATH", "")

from specex.io import read_lamp_lines, load_python_psf, read_image
from specex.fitter import get_bundle_spots, PSF_Fitter
from specex.specex import run_specex

def get_trace_rms(file_a, file_b):
    try:
        f_a = fitsio.FITS(file_a)
        f_b = fitsio.FITS(file_b)
        xt_a = f_a['XTRACE'].read(); xt_b = f_b['XTRACE'].read()
        yt_a = f_a['YTRACE'].read(); yt_b = f_b['YTRACE'].read()
        x_rms = np.std(xt_a - xt_b)
        y_rms = np.std(yt_a - yt_b)
        return x_rms, y_rms
    except Exception as e:
        print(f"Error comparing traces: {e}")
        return -1, -1

def run_comparison_suite(camera_list=None, night="20260401", expid="00344649"):
    bundle_id = 5
    sn_threshold = 3.0
    
    if camera_list is None:
        cameras = [f"z{i}" for i in range(10)]
    else:
        cameras = camera_list
    
    # Path templates
    arc_base = "/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/preproc"
    psf_base = "/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/exposures"
    lamp_lines = "py/specex/data/specex_linelist_desi.txt"
    
    print(f"--- Global Analysis: Night {night} Exp {expid} ---")
    print(f"{'Cam':<5} | {'Spots':<10} | {'Time (s)':<20} | {'Chi2':<20} | {'X-Trace RMS'}")
    print(f"{'':<5} | {'Py / C++':<10} | {'Py / C++':<20} | {'Py / C++':<20} |")
    print("-" * 85)
    
    for cam in cameras:
        arc_file = f"{arc_base}/{night}/{expid}/preproc-{cam}-{expid}.fits.gz"
        in_psf = f"{psf_base}/{night}/{expid}/shifted-input-psf-{cam}-{expid}.fits"
        
        # Check if files exist
        if not os.path.exists(arc_file):
            print(f"Skipping {cam}, {arc_file} not found.")
            continue

        broken = "367" if cam == "z0" else "473,474"
        
        # 1. Run Python GPU Fit
        # ... logic remains same ...
        # ... rest of logic remains same ...
        os.environ["JAX_PLATFORM_NAME"] = "gpu"
        py_out = f"validation_py_{cam}.fits"
        
        # We'll use a direct fit call to avoid mp overhead for measurement
        ddata = read_image(arc_file)
        image = ddata['image'].T; weight = ddata['ivar'].T
        
        # Dummy opts
        class Dummy: pass
        opts = Dummy(); opts.arc_image_filename = arc_file; opts.input_psf_filename = in_psf
        
        psf_py = load_python_psf(in_psf, opts)
        lines = read_lamp_lines(lamp_lines)
        
        fmin, fmax = bundle_id * 25, (bundle_id + 1) * 25 - 1
        spots = get_bundle_spots(psf_py, fmin, fmax, lines, image=image, weight=weight, sn_threshold=sn_threshold, broken_fibers=broken)
        n_spots_py = len(spots)
        
        fitter = PSF_Fitter(psf_py)
        t0 = time.time()
        # Warm-up (already done by first iteration in fit usually, but let's be explicit)
        chi2_py, pc, tc, cc = fitter.fit(image, weight, spots, bundle_id, max_iter=15)
        t1 = time.time()
        dt_py = t1 - t0
        
        # Save Python result
        from specex.io import write_python_psf
        write_python_psf(py_out, {bundle_id: {'chi2': chi2_py, 'psf_coeffs': pc, 'trace_coeffs': tc, 'continuum': cc}}, in_psf)

        # 2. Run C++ Fit (via Python wrapper using local build)
        cpp_out = f"validation_cpp_{cam}.fits"
        com = [
            "desi_psf_fit",
            "-a", arc_file,
            "--in-psf", in_psf,
            "--lamp-lines", lamp_lines,
            "--out-psf", cpp_out,
            "--first-bundle", str(bundle_id),
            "--last-bundle", str(bundle_id),
            "--legendre-deg-wave", "3",
            "--fit-continuum",
            "--broken-fibers", broken
        ]
        
        t0 = time.time()
        # Redirect stdout to capture spot count if needed, but for now just run
        # Actually run_specex writes to stdout
        run_specex(com)
        t1 = time.time()
        dt_cpp = t1 - t0
        
        # Extract C++ info from result file header
        f_cpp = fitsio.FITS(cpp_out)
        hdr = f_cpp['PSF'].read_header()
        chi2_cpp = hdr.get(f'B{bundle_id:02d}RCHI2', -1.0) * (120000.0) # Approx scale back to raw chi2
        # Actually, let's just get it from the fits comparison script if possible
        
        # 3. Compare
        x_rms, y_rms = get_trace_rms(py_out, cpp_out)
        
        # Quick hack to get spot count from debug if we could, but let's assume parity for now
        # or just read from the C++ log we just created.
        n_spots_cpp = "?" # would need to parse C++ output
        
        line = f"{cam:<5} | {n_spots_py:>4} / {n_spots_cpp:<3} | {dt_py:>7.1f} / {dt_cpp:<7.1f} | {chi2_py:>9.0f} / {chi2_cpp:<9.0f} | {x_rms:>10.6f}"
        print(line)

if __name__ == "__main__":
    arm = "z"
    if len(sys.argv) > 1:
        arm = sys.argv[1]
    
    cams = [f"{arm}{i}" for i in range(10)]
    run_comparison_suite(cams)
