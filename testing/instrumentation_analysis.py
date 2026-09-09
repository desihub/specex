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

# Per-band settings the real desispec `desi_compute_psf` wrapper passes to
# the C++ `desi_psf_fit` binary (desispec/scripts/specex.py:224-228) --
# z-band gets a higher wavelength-basis degree plus a continuum fit, b/r
# don't. The Python side auto-detects the same split from the CAMERA header
# (specex.specex.fit_ccd_native's own docstring) so it's left unset below
# and doesn't need duplicating here.
def cpp_band_flags(camera):
    band = camera[0]
    if band == "z":
        return ["--legendre-deg-wave", "3", "--fit-continuum"]
    return ["--legendre-deg-wave", "1"]


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

    # Try multiple common keys
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
        # print(f"Error comparing wavelengths: {e}")
        return -1.0, -1.0

def run_subprocess_fit(mode, arc_file, psf_file, broken_fibers, camera, bundle_id, sn_threshold=3.0, outdir="."):
    """
    Runs a fit in a subprocess to ensure clean backend initialization.
    Returns dict of metrics.
    """
    results = {'chi2': -1.0, 'time': 0.0, 'nspots': 0, 'dx_rms': 0.0, 'dy_rms': 0.0, 'error': None}

    lamp_lines_file = os.path.join(current_dir, 'py/specex/data/specex_linelist_desi.txt')
    out_file = os.path.join(outdir, f"fit_{mode}_{camera}_b{bundle_id}.fits")
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
        ] + cpp_band_flags(camera)
        if broken_fibers:
            cmd.extend(["--broken-fibers", broken_fibers])

        full_cmd = " ".join(cmd)
        print(f"  [CPP] Running: {full_cmd}")
        t0 = time.time()
        result = subprocess.run(full_cmd, capture_output=True, text=True, shell=True, executable="/bin/bash")
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
        # Goes through the real `python -m specex.specex` CLI (same entry
        # point documented in docs/python-port/how-to-run.md) rather than calling
        # PSF_Fitter.fit() directly -- that way this script automatically
        # tracks whatever specex.py's current production defaults are
        # (legendre-deg-wave/fit-continuum auto-detected per band,
        # trace-per-fiber-deg=6 + the ndead-gated trace prior as of
        # 2026-08-05, etc.) instead of silently drifting out of date every
        # time those defaults change, which is what happened to the
        # previous version of this script (it called PSF_Fitter.fit() with
        # none of those kwargs set, so it always tested pre-2026-08-05
        # shared-trace-basis behavior with z-band's wdeg/continuum
        # hardcoded for every band).
        backend = "cpu" if mode == "py_cpu" else "gpu"
        cmd = [
            sys.executable, "-m", "specex.specex",
            "-a", arc_file,
            "--in-psf", psf_file,
            "--out-psf", out_file,
            "--first-bundle", str(bundle_id),
            "--last-bundle", str(bundle_id),
            "--backend", backend,
            "--sn-threshold", str(sn_threshold),
            "--gpu", "1",
        ]
        if broken_fibers:
            cmd.extend(["--broken-fibers", broken_fibers])

        env = os.environ.copy()
        env["PYTHONPATH"] = os.path.join(current_dir, 'py') + ":" + env.get("PYTHONPATH", "")
        # Memory-sharing courtesy only (doesn't affect correctness) --
        # --backend already handles CUDA_VISIBLE_DEVICES/JAX_PLATFORMS
        # internally, see specex.py's fit_bundle_task.
        env["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
        env["XLA_PYTHON_CLIENT_MEM_FRACTION"] = ".70"

        print(f"  [{mode.upper()}] Running: {' '.join(cmd)}")
        t0 = time.time()
        result = subprocess.run(cmd, capture_output=True, text=True, env=env, cwd=current_dir)
        t1 = time.time()

        results['time'] = t1 - t0
        if result.returncode != 0:
            results['error'] = result.stderr
            return results

        out = result.stdout + result.stderr
        m_chi2 = re.findall(r'Iter \d+: chi2 = ([\d\.]+)', out)
        m_spots = re.findall(r'(?:Iterative spot selection|Spot selection) took [\d\.]+s \((\d+) spots\)', out)

        results['chi2'] = float(m_chi2[-1]) if m_chi2 else -1.0
        results['nspots'] = int(m_spots[-1]) if m_spots else 0
        return results

def main():
    parser = argparse.ArgumentParser(description="Full Instrumentation Analysis for Specex.")
    parser.add_argument("--night", type=str, default="20260401")
    parser.add_argument("--expid", type=str, default="00344649")
    parser.add_argument("--cameras", type=str, help="Comma-separated list (e.g. b0,r3,z8)")
    parser.add_argument("--bundle", type=int, default=5)
    parser.add_argument("--sn", type=float, default=3.0, help="S/N threshold for spot selection, passed as --sn-threshold to the Python CLI (default matches specex.py's own default)")
    parser.add_argument("--outdir", type=str, default=".", help="Directory for intermediate per-mode fit output FITS files")
    parser.add_argument("--output", type=str, default="instrumentation_analysis.txt")

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

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
                res = run_subprocess_fit(mode, case['image'], case['input_psf'], case['broken_fibers'], cam, args.bundle, sn_threshold=args.sn, outdir=args.outdir)
                modes_results[mode] = res

                xt_rms, yt_rms = 0.0, 0.0
                if mode == "py_gpu" and modes_results.get("cpp") and modes_results['cpp']['chi2'] > 0:
                    cpp_out = os.path.join(args.outdir, f"fit_cpp_{cam}_b{args.bundle}.fits")
                    py_out = os.path.join(args.outdir, f"fit_py_gpu_{cam}_b{args.bundle}.fits")
                    if os.path.exists(cpp_out) and os.path.exists(py_out):
                        xt_rms, yt_rms = get_wavelength_diff(cpp_out, py_out, args.bundle)

                line = f"{cam:<5} | {mode:<10} | {res['time']:>8.2f} | {res['chi2']:>12.1f} | {res['nspots']:>6} | {xt_rms:>8.4f} | {yt_rms:>8.4f}\n"
                f.write(line)
                f.flush()
                print(line, end="")
                if res['error']:
                    print(f"  ERROR: {res['error']}")

if __name__ == "__main__":
    main()
