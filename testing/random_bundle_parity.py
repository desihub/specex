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
    # Get header for wavelength range
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

def run_fit(mode, arc_file, psf_file, broken_fibers, camera, bundle_id, out_file, sn_threshold=3.0):
    """
    Runs either the C++ or Python fit for a specific bundle.
    """
    env = os.environ.copy()
    env["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
    env["XLA_PYTHON_CLIENT_MEM_FRACTION"] = ".70"
    
    lamp_lines_file = os.path.join(current_dir, 'py/specex/data/specex_linelist_desi.txt')
    
    if mode == "cpp":
        cmd = [
            "module", "load", "libfabric", "&&",
            "desi_psf_fit",
            "-a", arc_file,
            "--in-psf", psf_file,
            "--lamp-lines", lamp_lines_file,
            "--out-psf", out_file,
            "--first-bundle", str(bundle_id),
            "--last-bundle", str(bundle_id),
            "--legendre-deg-wave", "3",
            "--fit-continuum"
        ]
        if broken_fibers:
            cmd.extend(["--broken-fibers", broken_fibers])
        full_cmd = " ".join(cmd)
        result = subprocess.run(full_cmd, capture_output=True, text=True, env=env, shell=True, executable="/bin/bash")
        
        m_spots = re.search(r'selected\s+(\d+)\s+spots', result.stdout + result.stderr)
        nspots = int(m_spots.group(1)) if m_spots else 0
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
from specex.specex import fit_ccd_native

arc = '{arc_file}'
psf = '{psf_file}'
bid = {bundle_id}
broken = '{broken_fibers}'
lines_file = '{lamp_lines_file}'
out_file = '{out_file}'
sn = {sn_threshold}

fit_ccd_native(
    arc_file=arc,
    in_psf_file=psf,
    out_psf_file=out_file,
    lamp_lines_file=lines_file,
    first_bundle=bid,
    last_bundle=bid,
    n_gpus=1,
    backend="gpu",
    broken_fibers=broken,
    sn_threshold=sn
)

# Extract spots from the resulting PSF to save for RMS calculation
# (The original script saved spots manually; fit_ccd_native does the fit, but not the .spots file)
# We need to re-run spot selection or read from the PSF if available.
# For parity check, we can just use the PSF output.
'''
        result = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True, env=env)
        # The original script expected NSPOTS_RESULT. fit_ccd_native doesn't print it in that format.
        # We'll search for the "Spot selection took" line.
        m_spots = re.search(r'Spot selection took .*\((\d+) spots\)', result.stdout)
        nspots = int(m_spots.group(1)) if m_spots else 0
        
        # Since fit_ccd_native doesn't write .spots, we manually write them here for calculate_rms
        # to maintain compatibility with the existing RMS logic.
        spots_script = f'''
import os
import sys
sys.path.insert(0, '{os.path.join(current_dir, 'py')}')
from specex.io import read_preproc, load_python_psf, read_lamp_lines
from specex.fitter import get_bundle_spots

arc = '{arc_file}'
psf_file = '{psf_file}'
bid = {bundle_id}
broken = '{broken_fibers}'
lines_file = '{lamp_lines_file}'
out_spots_file = '{out_file}.spots'
sn = {sn_threshold}

class Opts:
    def __init__(self):
        self.arc_image_filename = arc
        self.input_psf_filename = psf_file
opts = Opts()
ddata = read_preproc(opts)
image = ddata['image'].T
weight = ddata['ivar'].T
psf_py = load_python_psf(psf_file, opts)
psf_py.h_size_y = 5
lamp_lines = read_lamp_lines(lines_file)

fmin, fmax = bid * 25, (bid + 1) * 25 - 1
spots = get_bundle_spots(psf_py, fmin, fmax, lamp_lines, image=image, weight=weight, sn_threshold=sn, broken_fibers=broken)

with open(out_spots_file, 'w') as f:
    for s in spots:
        f.write(f"{{s['fiber']}},{{s['wave']}}\\n")
'''
        subprocess.run([sys.executable, "-c", spots_script], env=env)
        return nspots

def calculate_rms(cpp_file, py_file, spots_file):
    """
    Calculates the RMS difference between C++ and Py positions 
    at the wavelengths identified by the Python fit.
    """
    if not os.path.exists(spots_file):
        return -1.0, -1.0

    diffs_x, diffs_y = [], []
    with open(spots_file, 'r') as f:
        for line in f:
            try:
                parts = line.strip().split(',')
                fiber = int(parts[0])
                wave = float(parts[1])
                x_cpp, y_cpp = get_predicted_pos(cpp_file, fiber, wave)
                x_py, y_py = get_predicted_pos(py_file, fiber, wave)
                diffs_x.append(x_cpp - x_py)
                diffs_y.append(y_cpp - y_py)
            except Exception:
                continue
                
    if not diffs_x:
        return -1.0, -1.0
        
    return np.sqrt(np.mean(np.array(diffs_x)**2)), np.sqrt(np.mean(np.array(diffs_y)**2))

def main():
    parser = argparse.ArgumentParser(description="Random Bundle Parity Check")
    parser.add_argument("--n-bundles", type=int, default=1, help="Bundles per camera to compare")
    parser.add_argument("--night", type=str, default="20260401")
    parser.add_argument("--output-dir", type=str, default="/pscratch/sd/c/cdwarner/specex/random_bundles/")
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
    
    results_file = os.path.join(out_dir, "bundle_parity_results.txt")
    file_exists = os.path.exists(results_file)
    with open(results_file, "a") as f:
        if not file_exists:
            header = f"{'Night':<10} | {'ExpID':<10} | {'Cam':<5} | {'Bundle':<8} | {'CPP Spots':<10} | {'Py Spots':<10} | {'DX RMS':<10} | {'DY RMS':<10}\n"
            f.write(header)
            f.write("-" * 85 + "\n")
            print(header, end="")
        else:
            print("Appending to existing results file...")

        # Select random cases across the whole night
        random.shuffle(cases)
        
        # Group by camera "color" (first char of cam name, e.g. 'b', 'r', 'z')
        color_groups = {}
        for c in cases:
            color = c['camera'][0]
            if color not in color_groups:
                color_groups[color] = []
            color_groups[color].append(c)
        
        for color in sorted(color_groups.keys()):
            cam_cases = color_groups[color]
            sample_size = min(len(cam_cases), args.n_bundles)
            selected_cases = random.sample(cam_cases, sample_size)
            
            for case in selected_cases:
                bundle_id = random.randint(0, 19)
                cam_name = case['camera']
                expid = case['expid']
                night = args.night
                
                cpp_out = os.path.join(out_dir, f"fit_cpp_{cam_name}_{expid}_b{bundle_id}.fits")
                py_out = os.path.join(out_dir, f"fit_py_{cam_name}_{expid}_b{bundle_id}.fits")
                spots_file = py_out + '.spots'
                
                print(f"Processing {night} | {expid} | {cam_name} | Bundle {bundle_id}...")
                
                n_cpp = run_fit("cpp", case['image'], case['input_psf'], case['broken_fibers'], cam_name, bundle_id, cpp_out)
                n_py = run_fit("py", case['image'], case['input_psf'], case['broken_fibers'], cam_name, bundle_id, py_out)
                
                dx, dy = calculate_rms(cpp_out, py_out, spots_file)
                
                line = f"{night:<10} | {expid:<10} | {cam_name:<5} | {bundle_id:<8} | {n_cpp:<10} | {n_py:<10} | {dx:>10.4f} | {dy:>10.4f}\n"
                f.write(line)
                f.flush()
                print(line, end="")

if __name__ == "__main__":
    main()
