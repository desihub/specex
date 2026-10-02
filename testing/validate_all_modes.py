import os
import sys
import time
import argparse
import subprocess
import re
import glob

# Ensure we use the current workspace code
current_dir = os.getcwd()

def run_subprocess_fit(mode, arc_file, psf_file, broken_fibers, camera, bundle_id):
    """
    Runs a fit in a subprocess to ensure clean backend initialization.
    Returns (chi2, time).
    """
    env = os.environ.copy()
    
    if mode == "cpp":
        # Don't add local py to path for CPP, use production environment
        # But we need the lamp lines file from our local data if possible, 
        # or use production one. Let's try production one if available.
        lamp_lines = os.path.join(current_dir, 'py/specex/data/specex_linelist_desi.txt')
        out_file = f"baseline_cpp_{camera}_b{bundle_id}.fits"
        
        cmd = [
            "module", "load", "libfabric", "&&",
            "desi_psf_fit",
            "-a", arc_file,
            "--in-psf", psf_file,
            "--lamp-lines", lamp_lines,
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
        
        if result.returncode != 0:
            print(f"  [CPP] Failed: {result.stderr}")
            return None, 0.0
            
        # Parse Chi2 from output
        match = re.search(r'chi2=\s*([\d\.]+)', result.stdout + result.stderr)
        chi2 = float(match.group(1)) if match else -1.0
        return chi2, t1 - t0

    elif mode in ["py_cpu", "py_gpu"]:
        env["JAX_PLATFORM_NAME"] = "cpu" if mode == "py_cpu" else "gpu"
        env["PYTHONPATH"] = os.path.join(current_dir, 'py') + ":" + env.get("PYTHONPATH", "")
        
        script = f"""
import sys, os
sys.path.insert(0, '{os.path.join(current_dir, 'py')}')
from specex.io import read_preproc, load_python_psf, read_lamp_lines
from specex.fitter import PSF_Fitter, get_bundle_spots
import numpy as np

arc = '{arc_file}'
psf = '{psf_file}'
bid = {bundle_id}
broken = '{broken_fibers}'
lines_file = '{os.path.join(current_dir, 'py/specex/data/specex_linelist_desi.txt')}'

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
spots = get_bundle_spots(psf_py, fmin, fmax, lamp_lines, image=image, weight=weight, sn_threshold=3.0, broken_fibers=broken)
fitter = PSF_Fitter(psf_py)
chi2, pc, tc, cc = fitter.fit(image, weight, spots, bid, max_iter=50)
print(f"CHI2_RESULT: {{chi2}}")
"""
        t0 = time.time()
        result = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True, env=env)
        t1 = time.time()
        
        if result.returncode != 0:
            print(f"  [{mode}] Failed: {result.stderr}")
            return None, 0.0
            
        match = re.search(r'CHI2_RESULT:\s*([\d\.]+)', result.stdout)
        chi2 = float(match.group(1)) if match else -1.0
        return chi2, t1 - t0

def main():
    parser = argparse.ArgumentParser(description="Full 3-mode validator for Specex.")
    parser.add_argument("--night", type=str, default="20260401")
    parser.add_argument("--expid", type=str, default="00344649")
    parser.add_argument("--cameras", type=str, help="Comma-separated list (e.g. b0,r3,z8)")
    parser.add_argument("--bundle", type=int, default=5)
    parser.add_argument("--output", type=str, default="full_comparison_results.txt")
    
    args = parser.parse_args()
    
    # 1. Scrape logs for parameters
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
        header = f"{'Cam':<5} | {'Mode':<12} | {'Time(s)':<10} | {'Chi2':<12} | {'Speedup'}\n"
        f.write(header)
        f.write("-" * 60 + "\n")
        print(header, end="")
        
        for cam in selected_cams:
            if cam not in case_map: continue
            case = case_map[cam]
            
            results = {}
            for mode in ["cpp", "py_cpu", "py_gpu"]:
                print(f"Running {cam} in {mode} mode...")
                chi2, dt = run_subprocess_fit(mode, case['image'], case['input_psf'], case['broken_fibers'], cam, args.bundle)
                results[mode] = (chi2, dt)
                
                speedup = results['cpp'][1] / dt if results.get('cpp') and results['cpp'][1] > 0 and dt > 0 else 1.0
                line = f"{cam:<5} | {mode:<12} | {dt:>10.2f} | {chi2 if chi2 is not None else -1.0:>12.1f} | {speedup:>7.2f}x\n"
                f.write(line)
                f.flush()
                print(line, end="")

if __name__ == "__main__":
    main()
