import os
import sys
import time
import argparse
import subprocess
import re
import glob
import numpy as np
import fitsio

# Ensure we use the current workspace code
current_dir = os.getcwd()
sys.path.insert(0, os.path.join(current_dir, 'py'))
from specex.math import Legendre1DPol

def evaluate_wavelength(psf_file, fiber):
    """
    Evaluates the X and Y positions for a fiber over its wavelength range.
    Returns (wave_grid, x_vals, y_vals).
    """
    f = fitsio.FITS(psf_file)
    if 'PSF' in f:
        hdr = f['PSF'].read_header()
    else:
        hdr = f['XTRACE'].read_header()
        
    xt = f['XTRACE'].read()[fiber]
    yt = f['YTRACE'].read()[fiber]
    
    wmin = hdr.get('WAVEMIN', hdr.get('WAVE_MIN', 3500.0))
    wmax = hdr.get('WAVEMAX', hdr.get('WAVE_MAX', 10000.0))
    wave_grid = np.linspace(float(wmin), float(wmax), 100)
    
    poly_x = Legendre1DPol(deg=len(xt)-1, xmin=wmin, xmax=wmax, coeff=xt)
    poly_y = Legendre1DPol(deg=len(yt)-1, xmin=wmin, xmax=wmax, coeff=yt)
    
    return wave_grid, np.array(poly_x.value(wave_grid)), np.array(poly_y.value(wave_grid))

def get_wavelength_diff(psf1_file, psf2_file, bundle_id):
    """
    Compares wavelength solution by evaluating positions on a grid.
    Returns (dx_rms, dy_rms) in pixels.
    """
    try:
        fmin, fmax = bundle_id * 25, (bundle_id + 1) * 25 - 1
        all_dx = []
        all_dy = []
        
        for fib in range(fmin, fmax + 1):
            w1, x1, y1 = evaluate_wavelength(psf1_file, fib)
            w2, x2, y2 = evaluate_wavelength(psf2_file, fib)
            all_dx.append(x1 - x2)
            all_dy.append(y1 - y2)
            
        dx_rms = np.sqrt(np.mean(np.array(all_dx)**2))
        dy_rms = np.sqrt(np.mean(np.array(all_dy)**2))
        
        return dx_rms, dy_rms
    except Exception as e:
        print(f"Error comparing wavelengths: {e}")
        return -1.0, -1.0

def run_subprocess_fit(mode, arc_file, psf_file, broken_fibers, camera, bundle_id, sn_threshold=5.0):
    """
    Runs a fit in a subprocess to ensure clean backend initialization.
    Returns dict of metrics.
    """
    env = os.environ.copy()
    results = {'chi2': -1.0, 'time': 0.0, 'nspots': 0, 'dx_rms': 0.0, 'dy_rms': 0.0, 'error': None}
    
    lamp_lines_file = os.path.join(current_dir, 'py/specex/data/specex_linelist_desi.txt')
    out_file = f"fit_{mode}_{camera}_b{bundle_id}.fits"
    if os.path.exists(out_file): os.remove(out_file)

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
        print(f"  [CPP] Running: {full_cmd}")
        t0 = time.time()
        result = subprocess.run(full_cmd, capture_output=True, text=True, env=env, shell=True, executable="/bin/bash")
        t1 = time.time()
        
        results['time'] = t1 - t0
        if result.returncode != 0:
            results['error'] = result.stderr
            return results
            
        m_chi2 = re.findall(r'chi2=\s*([\d\.]+)', result.stdout + result.stderr)
        m_spots = re.findall(r'selected\s+(\d+)\s+spots', result.stdout + result.stderr)
        
        results['chi2'] = float(m_chi2[-1]) if m_chi2 else -1.0
        results['nspots'] = int(m_spots[-1]) if m_spots else 0
        return results

    elif mode in ["py_cpu", "py_gpu"]:
        env["JAX_PLATFORM_NAME"] = "cpu" if mode == "py_cpu" else "gpu"
        env["PYTHONPATH"] = os.path.join(current_dir, 'py') + ":" + env.get("PYTHONPATH", "")
        
        script = f"""
import sys, os
sys.path.insert(0, '{os.path.join(current_dir, 'py')}')
from specex.io import read_preproc, load_python_psf, read_lamp_lines, write_python_psf
from specex.fitter import PSF_Fitter, get_bundle_spots
import numpy as np

arc = '{arc_file}'
psf = '{psf_file}'
bid = {bundle_id}
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

fmin, fmax = bid * 25, (bid + 1) * 25 - 1
spots = get_bundle_spots(psf_py, fmin, fmax, lamp_lines, image=image, weight=weight, sn_threshold=sn, broken_fibers=broken)
print(f"NSPOTS_RESULT: {{len(spots)}}")

fitter = PSF_Fitter(psf_py)
chi2, pc, tc, cc = fitter.fit(image, weight, spots, bid, max_iter=50)
print(f"CHI2_RESULT: {{chi2}}")

res_map = {{bid: {{'chi2': chi2, 'psf_coeffs': pc, 'trace_coeffs': tc, 'continuum': cc}}}}
write_python_psf(out_file, res_map, psf)
"""
        t0 = time.time()
        result = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True, env=env)
        t1 = time.time()
        
        results['time'] = t1 - t0
        if result.returncode != 0:
            results['error'] = result.stderr
            return results
            
        m_chi2 = re.search(r'CHI2_RESULT:\s*([\d\.]+)', result.stdout)
        m_spots = re.search(r'NSPOTS_RESULT:\s*(\d+)', result.stdout)
        
        results['chi2'] = float(m_chi2.group(1)) if m_chi2 else -1.0
        results['nspots'] = int(m_spots.group(1)) if m_spots else 0
        return results

def main():
    parser = argparse.ArgumentParser(description="Full Instrumentation Analysis for Specex.")
    parser.add_argument("--night", type=str, default="20260401")
    parser.add_argument("--expid", type=str, default="00344649")
    parser.add_argument("--cameras", type=str, help="Comma-separated list (e.g. b0,r3,z8)")
    parser.add_argument("--bundle", type=int, default=5)
    parser.add_argument("--sn", type=float, default=5.0)
    parser.add_argument("--output", type=str, default="instrumentation_analysis.txt")
    
    args = parser.parse_args()
    
    # Check for select_test_case logic
    from select_test_case import parse_log_line
    log_dir = f"/global/cfs/cdirs/desi/spectro/redux/matterhorn/run/scripts/night/{args.night}"
    log_pattern = os.path.join(log_dir, "arc*.log")
    logs = glob.glob(log_pattern)
    
    case_map = {}
    for log in logs:
        with open(log, 'r') as f:
            for line in f:
                if "desi_compute_psf" in line and args.expid in line:
                    case = parse_log_line(line)
                    if case: case_map[case['camera']] = case
                    
    selected_cams = args.cameras.split(',') if args.cameras else sorted(case_map.keys())
    
    with open(args.output, "w") as f:
        header = f"{'Cam':<5} | {'Mode':<10} | {'Time(s)':<8} | {'Chi2':<12} | {'Spots':<6} | {'XT RMS':<8} | {'YT RMS':<8}\n"
        f.write(header)
        f.write("-" * 85 + "\n")
        print(header, end="")
        
        for cam in selected_cams:
            if cam not in case_map: continue
            case = case_map[cam]
            print(f"\n--- Testing Camera {cam} ---")
            
            modes_results = {}
            for mode in ["cpp", "py_gpu"]:
                res = run_subprocess_fit(mode, case['image'], case['input_psf'], case['broken_fibers'], cam, args.bundle, sn_threshold=args.sn)
                modes_results[mode] = res
                
                xt_rms, yt_rms = 0.0, 0.0
                if mode == "py_gpu" and modes_results.get("cpp") and modes_results['cpp']['chi2'] > 0:
                    cpp_out = f"fit_cpp_{cam}_b{args.bundle}.fits"
                    py_out = f"fit_py_gpu_{cam}_b{args.bundle}.fits"
                    if os.path.exists(cpp_out) and os.path.exists(py_out):
                        xt_rms, yt_rms = get_wavelength_diff(cpp_out, py_out, args.bundle)
                
                line = f"{cam:<5} | {mode:<10} | {res['time']:>8.2f} | {res['chi2']:>12.1f} | {res['nspots']:>6} | {xt_rms:>8.4f} | {yt_rms:>8.4f}\n"
                f.write(line)
                f.flush()
                print(line, end="")
                if res['error']:
                    print(f"  ERROR: {res['error'][:500]}")

if __name__ == "__main__":
    main()
