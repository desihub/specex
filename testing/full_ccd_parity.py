import os
import sys
import time
import argparse
import subprocess
import re
import glob
import random
import numpy as np
import fitsio

# Ensure we use the current workspace code
current_dir = os.getcwd()
sys.path.insert(0, os.path.join(current_dir, 'py'))
from specex.math import Legendre1DPol

def get_predicted_pos(psf_file, fiber, wave):
    """
    Predicts the X, Y pixel position for a given fiber and wavelength.
    """
    f = fitsio.FITS(psf_file)
    if 'PSF' in f:
        hdr = f['PSF'].read_header()
    elif 'XTRACE' in f:
        hdr = f['XTRACE'].read_header()
    else:
        raise KeyError("FITS file missing required extensions for trace evaluation")
        
    wmin = hdr.get('WAVEMIN', hdr.get('WAVE_MIN', 3500.0))
    wmax = hdr.get('WAVEMAX', hdr.get('WAVE_MAX', 10000.0))
    
    xt = f['XTRACE'].read()[fiber]
    yt = f['YTRACE'].read()[fiber]
    
    poly_x = Legendre1DPol(deg=len(xt)-1, xmin=wmin, xmax=wmax, coeff=xt)
    poly_y = Legendre1DPol(deg=len(yt)-1, xmin=wmin, xmax=wmax, coeff=yt)
    
    return poly_x.value(wave), poly_y.value(wave)

def run_full_fit(mode, arc_file, psf_file, broken_fibers, camera, out_file, sn_threshold=3.0):
    """
    Runs a full CCD fit for the specified mode.
    """
    env = os.environ.copy()
    env["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
    env["XLA_PYTHON_CLIENT_MEM_FRACTION"] = ".70"
    
    lamp_lines_file = os.path.join(current_dir, 'py/specex/data/specex_linelist_desi.txt')
    
    if mode == "cpp":
        cmd = [
            "module", "load", "libfabric", "&&",
            "desi_compute_psf",
            "-a", arc_file,
            "--in-psf", psf_file,
            "--out-psf", out_file,
            "--lamp-lines", lamp_lines_file,
        ]
        if broken_fibers:
            cmd.extend(["--broken-fibers", broken_fibers])
        full_cmd = " ".join(cmd)
        result = subprocess.run(full_cmd, capture_output=True, text=True, env=env, shell=True, executable="/bin/bash")
        
        m_spots = re.findall(r'selected\s+(\d+)\s+spots', result.stdout + result.stderr)
        nspots = sum(map(int, m_spots)) if m_spots else 0
        return nspots

    elif mode == "py":
        env["JAX_PLATFORM_NAME"] = "gpu"
        env["PYTHONPATH"] = os.path.join(current_dir, 'py') + ":" + env.get("PYTHONPATH", "")
        
        script = f'''
import os
import sys
import numpy as np
import fitsio
sys.path.insert(0, '{os.path.join(current_dir, 'py')}')
from specex.io import read_preproc, load_python_psf, read_lamp_lines, write_python_psf
from specex.fitter import PSF_Fitter, get_bundle_spots

arc = '{arc_file}'
psf = '{psf_file}'
broken = '{broken_fibers}'
lines_file = '{lamp_lines_file}'
out_file = '{out_file}'
sn = {sn_threshold}

class Opts:
    def __init__(self):
        self.arc_image_filename = arc
        self.input_psf_filename = psf
opts = Opts()
ddata = read_preproc(opts)
image = ddata['image'].T
weight = ddata['ivar'].T
psf_py = load_python_psf(psf, opts)
psf_py.h_size_y = 5
lamp_lines = read_lamp_lines(lines_file)

all_results = {{}}
total_spots = 0
for bid in range(20):
    fmin, fmax = bid * 25, (bid + 1) * 25 - 1
    spots = get_bundle_spots(psf_py, fmin, fmax, lamp_lines, image=image, weight=weight, sn_threshold=sn, broken_fibers=broken)
    total_spots += len(spots)
    fitter = PSF_Fitter(psf_py)
    chi2, pc, tc, cc = fitter.fit(image, weight, spots, bid, max_iter=50)
    all_results[bid] = {{'chi2': chi2, 'psf_coeffs': pc, 'trace_coeffs': tc, 'continuum': cc}}

write_python_psf(out_file, all_results, psf)
print(f"NSPOTS_RESULT: {{total_spots}}")
'''
        result = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True, env=env)
        m_spots = re.search(r'NSPOTS_RESULT:\s*(\d+)', result.stdout)
        nspots = int(m_spots.group(1)) if m_spots else 0
        return nspots

def calculate_full_rms(cpp_file, py_file):
    """
    Calculates the RMS difference across all fibers and bundles.
    Since full fits use a global model, we evaluate the trace for all 500 fibers.
    """
    diffs_x, diffs_y = [], []
    for fib in range(500):
        # We use a representative grid of wavelengths for each fiber
        f = fitsio.FITS(py_file)
        hdr = f['PSF'].read_header() if 'PSF' in f else f['XTRACE'].read_header()
        wmin = hdr.get('WAVEMIN', hdr.get('WAVE_MIN', 3500.0))
        wmax = hdr.get('WAVEMAX', hdr.get('WAVE_MAX', 10000.0))
        wave_grid = np.linspace(wmin, wmax, 50)
        
        for w in wave_grid:
            try:
                x_cpp, y_cpp = get_predicted_pos(cpp_file, fib, w)
                x_py, y_py = get_predicted_pos(py_file, fib, w)
                diffs_x.append(x_cpp - x_py)
                diffs_y.append(y_cpp - y_py)
            except Exception:
                continue
                
    return np.sqrt(np.mean(np.array(diffs_x)**2)), np.sqrt(np.mean(np.array(diffs_y)**2))

def calculate_wavecorr_rms(cpp_file, py_file):
    """
    Compares the 4th extension (WAVECORR) of the FITS files.
    """
    try:
        f_cpp = fitsio.FITS(cpp_file)
        f_py = fitsio.FITS(py_file)
        
        # HDU 3 is the 4th extension
        data_cpp = f_cpp[3].read()
        data_py = f_py[3].read()
        
        # We only compare rows where wavelengths match
        waves_cpp = data_cpp['WAVE']
        waves_py = data_py['WAVE']
        
        # Find matching wavelengths
        diffs_dwave = []
        for i, w_py in enumerate(waves_py):
            # Find closest matching wave in cpp
            idx = np.argmin(np.abs(waves_cpp - w_py))
            if np.abs(waves_cpp[idx] - w_py) < 0.1:
                diffs_dwave.append(data_cpp['DWAVE'][idx] - data_py['DWAVE'][i])
        
        if not diffs_dwave:
            return -1.0
            
        return np.sqrt(np.mean(np.array(diffs_dwave)**2))
    except Exception as e:
        # print(f"WAVECORR error: {e}")
        return -1.0

def main():
    parser = argparse.ArgumentParser(description="Full CCD Parity Check")
    parser.add_argument("--n-cases", type=int, default=5, help="Number of full CCDs per camera to compare")
    parser.add_argument("--night", type=str, default="20260401")
    parser.add_argument("--output-dir", type=str, default="/pscratch/sd/c/cdwarner/specex/full_ccd_parity/")
    parser.add_argument("--sn", type=float, default=3.0)
    args = parser.parse_args()
    
    out_dir = args.output_dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
    from select_test_case import parse_log_line
    log_dir = f"/global/cfs/cdirs/desi/spectro/redux/matterhorn/run/scripts/night/{args.night}"
    log_pattern = os.path.join(log_dir, "arc*.log")
    logs = glob.glob(log_pattern)
    
    cases = []
    for log in logs:
        with open(log, 'r') as f:
            for line in f:
                if "desi_compute_psf" in line:
                    case = parse_log_line(line)
                    if case: cases.append(case)
    
    unique_cases = { (c['camera'], c['image']): c for c in cases }
    cases = list(unique_cases.values())
    
    results_file = os.path.join(out_dir, "full_ccd_parity_results.txt")
    with open(results_file, "w") as f:
        header = f"{'Night':<10} | {'ExpID':<10} | {'Cam':<5} | {'CPP Spots':<10} | {'Py Spots':<10} | {'DX RMS':<10} | {'DY RMS':<10} | {'DWAVE RMS':<10}\\n"
        f.write(header)
        f.write("-" * 95 + "\\n")
        print(header, end="")

        random.shuffle(cases)
        cameras = sorted(list({c['camera'] for c in cases}))
        for cam in cameras:
            cam_cases = [c for c in cases if c['camera'] == cam]
            sample_size = min(len(cam_cases), args.n_cases)
            selected_cases = random.sample(cam_cases, sample_size)
            
            for case in selected_cases:
                cam_name = case['camera']
                expid = case['expid']
                night = args.night
                
                cpp_out = os.path.join(out_dir, f"full_cpp_{cam_name}_{expid}.fits")
                py_out = os.path.join(out_dir, f"full_py_{cam_name}_{expid}.fits")
                
                print(f"Processing Full CCD {night} | {expid} | {cam_name}...")
                
                n_cpp = run_full_fit("cpp", case['image'], case['input_psf'], case['broken_fibers'], cam_name, cpp_out)
                n_py = run_full_fit("py", case['image'], case['input_psf'], case['broken_fibers'], cam_name, py_out)
                
                dx, dy = calculate_full_rms(cpp_out, py_out)
                dwave_rms = calculate_wavecorr_rms(cpp_out, py_out)
                
                line = f"{night:<10} | {expid:<10} | {cam_name:<5} | {n_cpp:<10} | {n_py:<10} | {dx:>10.4f} | {dy:>10.4f} | {dwave_rms:>10.4f}\\n"
                f.write(line)
                f.flush()
                print(line, end="")

if __name__ == "__main__":
    main()
