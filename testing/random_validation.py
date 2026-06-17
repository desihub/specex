import os
import sys
import time
import random
import numpy as np
import fitsio
import subprocess
import argparse
import glob
from concurrent.futures import ProcessPoolExecutor, as_completed

# Ensure we use the current workspace code
current_dir = os.getcwd()
sys.path.insert(0, os.path.join(current_dir, 'py'))
os.environ["PYTHONPATH"] = os.path.join(current_dir, 'py') + ":" + os.path.join(current_dir, 'build') + ":" + os.environ.get("PYTHONPATH", "")

# Fix for libfabric on Perlmutter
libfabric_path = "/opt/cray/libfabric/1.22.0/lib64"
if os.path.exists(libfabric_path):
    os.environ["LD_LIBRARY_PATH"] = libfabric_path + ":" + os.environ.get("LD_LIBRARY_PATH", "")

def get_trace_rms(file_a, file_b):
    """Calculates the X and Y RMS difference between two PSF files."""
    if not os.path.exists(file_a) or not os.path.exists(file_b):
        return -1, -1
    try:
        f_a = fitsio.FITS(file_a)
        f_b = fitsio.FITS(file_b)
        
        # X-Trace
        xt_a = f_a['XTRACE'].read()
        xt_b = f_b['XTRACE'].read()
        mask_a = np.sum(np.abs(xt_a), axis=1) > 0
        mask_b = np.sum(np.abs(xt_b), axis=1) > 0
        mask = mask_a & mask_b
        x_rms = np.std(xt_a[mask] - xt_b[mask]) if np.any(mask) else -1
        
        # Y-Trace
        yt_a = f_a['YTRACE'].read()
        yt_b = f_b['YTRACE'].read()
        y_rms = np.std(yt_a[mask] - yt_b[mask]) if np.any(mask) else -1
        
        return x_rms, y_rms
    except:
        return -1, -1

def write_summary(outdir, data):
    """Appends a result row to the summary.txt file."""
    summary_file = os.path.join(outdir, "validation_summary.txt")
    header = "# Mode Night ExpID Cam Bundle Time Spots X-RMS Y-RMS\n"
    write_header = not os.path.exists(summary_file)
    
    with open(summary_file, "a") as f:
        if write_header:
            f.write(header)
        line = f"{data['mode']:10} {data['night']} {data['expid']:8} {data['cam']:5} {data['bundle']:6} {data['time']:8.2f} {data['spots']:6} {data['x_rms']:10.6f} {data['y_rms']:10.6f}\n"
        f.write(line)

def get_random_night_exp():
    """Finds a random night."""
    preproc_base = "/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/preproc"
    if not os.path.exists(preproc_base): return None
    nights = sorted([n for n in os.listdir(preproc_base) if n.isdigit() and int(n) > 20210101])
    if not nights: return None
    return random.choice(nights)

def run_bundle_test(args):
    night, expid, cam, bundle_id, gpu_id, outdir = args
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
    
    arc_file = f"/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/preproc/{night}/{expid}/preproc-{cam}-{expid}.fits.gz"
    in_psf = f"/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/exposures/{night}/{expid}/shifted-input-psf-{cam}-{expid}.fits"
    lamp_lines = "py/specex/data/specex_linelist_desi.txt"
    py_out = os.path.join(outdir, f"rand_py_{night}_{expid}_{cam}_b{bundle_id:02d}.fits")
    cpp_out = os.path.join(outdir, f"rand_cpp_{night}_{expid}_{cam}_b{bundle_id:02d}.fits")
    broken = "367" if cam == "z0" else "473,474"
    
    print(f"[Worker] Starting Bundle Fit: {night}/{expid} {cam} B{bundle_id}")

    # 1. Python Fit
    n_spots = 0
    dt_py = -1
    try:
        from specex.io import read_lamp_lines, load_python_psf, read_image, write_python_psf
        from specex.fitter import get_bundle_spots, PSF_Fitter
        ddata = read_image(arc_file)
        class Dummy: pass
        opts = Dummy(); opts.arc_image_filename = arc_file; opts.input_psf_filename = in_psf
        psf_py = load_python_psf(in_psf, opts)
        lines = read_lamp_lines(lamp_lines)
        fmin, fmax = bundle_id * 25, (bundle_id + 1) * 25 - 1
        spots = get_bundle_spots(psf_py, fmin, fmax, lines, image=ddata['image'].T, weight=ddata['ivar'].T, sn_threshold=3.0, broken_fibers=broken)
        n_spots = len(spots)
        fitter = PSF_Fitter(psf_py)
        t0 = time.time()
        chi2_py, pc, tc, cc = fitter.fit(ddata['image'].T, ddata['ivar'].T, spots, bundle_id, max_iter=15)
        dt_py = time.time() - t0
        write_python_psf(py_out, {bundle_id: {'chi2': chi2_py, 'psf_coeffs': pc, 'trace_coeffs': tc, 'continuum': cc}}, in_psf)
    except Exception as e: print(f"Python Bundle Failed for {cam} B{bundle_id}: {e}")

    # 2. C++ Fit
    dt_cpp = -1
    try:
        env_cpp = os.environ.copy()
        env_cpp["OMP_NUM_THREADS"] = "16"
        com_cpp = ["desi_psf_fit", "-a", arc_file, "--in-psf", in_psf, "--lamp-lines", lamp_lines, "--out-psf", cpp_out, "--first-bundle", str(bundle_id), "--last-bundle", str(bundle_id), "--fit-continuum", "--legendre-deg-wave", "3", "--broken-fibers", broken]
        t0 = time.time()
        res = subprocess.run(com_cpp, env=env_cpp, check=True, capture_output=True, timeout=300)
        dt_cpp = time.time() - t0
    except Exception as e: 
        print(f"C++ Bundle Failed for {cam} B{bundle_id}: {e}")
        if hasattr(e, 'stderr'): print(f"Stderr: {e.stderr.decode()}")

    # Calculate Parity
    x_rms, y_rms = get_trace_rms(py_out, cpp_out)
    
    if dt_py > 0:
        write_summary(outdir, {'mode': 'python', 'night': night, 'expid': expid, 'cam': cam, 'bundle': str(bundle_id), 'time': dt_py, 'spots': n_spots, 'x_rms': x_rms, 'y_rms': y_rms})
    if dt_cpp > 0:
        write_summary(outdir, {'mode': 'cpp', 'night': night, 'expid': expid, 'cam': cam, 'bundle': str(bundle_id), 'time': dt_cpp, 'spots': n_spots, 'x_rms': 0, 'y_rms': 0})
        
    return {"night": night, "expid": expid, "cam": cam, "bundle": bundle_id, "dt_py": dt_py, "dt_cpp": dt_cpp, "x_rms": x_rms, "y_rms": y_rms}

def run_ccd_test(args):
    night, expid, cam, gpu_id, outdir = args
    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    broken = "367" if cam == "z0" else "473,474"
    arc_file = f"/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/preproc/{night}/{expid}/preproc-{cam}-{expid}.fits.gz"
    in_psf = f"/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/exposures/{night}/{expid}/shifted-input-psf-{cam}-{expid}.fits"
    py_out = os.path.join(outdir, f"rand_py_{night}_{expid}_{cam}_full.fits")
    cpp_out = os.path.join(outdir, f"rand_cpp_{night}_{expid}_{cam}_full.fits")
    
    print(f"[Worker] Starting Full CCD Test: {night}/{expid} {cam}")

    # 1. Python CCD
    dt_py = -1
    try:
        t0 = time.time()
        # Ensure we run from the project root so -m specex works
        res = subprocess.run(["python", "-m", "specex.specex", "-a", arc_file, "--in-psf", in_psf, "--out-psf", py_out, "--gpu", "1", "--broken-fibers", broken], env=env, check=True, capture_output=True, timeout=1800)
        dt_py = time.time() - t0
    except Exception as e: 
        print(f"Python CCD Failed for {cam}: {e}")
        if hasattr(e, 'stderr'): print(f"Stderr: {e.stderr.decode()}")
    
    # 2. C++ CCD
    dt_cpp = -1
    try:
        t0 = time.time()
        if "SLURM_JOB_ID" in os.environ:
            com_cpp = ["srun", "-n", "20", "--cpu-bind", "none", "desi_compute_psf", "--mpi", "--input-image", arc_file, "--input-psf", in_psf, "--output-psf", cpp_out, "--broken-fibers", broken]
        else:
            env["OMP_NUM_THREADS"] = "32"
            com_cpp = ["desi_compute_psf", "--input-image", arc_file, "--input-psf", in_psf, "--output-psf", cpp_out, "--broken-fibers", broken]
        res = subprocess.run(com_cpp, env=env, check=True, capture_output=True, timeout=3600)
        dt_cpp = time.time() - t0
    except Exception as e: 
        print(f"C++ CCD Failed for {cam}: {e}")
        if hasattr(e, 'stderr'): print(f"Stderr: {e.stderr.decode()}")
    
    x_rms, y_rms = get_trace_rms(py_out, cpp_out)
    
    if dt_py > 0:
        write_summary(outdir, {'mode': 'python-ccd', 'night': night, 'expid': expid, 'cam': cam, 'bundle': 'ALL', 'time': dt_py, 'spots': 0, 'x_rms': x_rms, 'y_rms': y_rms})
    if dt_cpp > 0:
        write_summary(outdir, {'mode': 'cpp-ccd', 'night': night, 'expid': expid, 'cam': cam, 'bundle': 'ALL', 'time': dt_cpp, 'spots': 0, 'x_rms': 0, 'y_rms': 0})

    return {"night": night, "expid": expid, "cam": cam, "dt_py": dt_py, "dt_cpp": dt_cpp, "x_rms": x_rms, "y_rms": y_rms}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--night", type=str, default=None)
    parser.add_argument("--expid", type=str, default=None)
    parser.add_argument("--bundles-per-arm", type=int, default=0)
    parser.add_argument("--ccds-per-arm", type=int, default=0)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--outdir", type=str, default=".")
    args = parser.parse_args()
    
    if not os.path.exists(args.outdir): os.makedirs(args.outdir)
    random.seed(int(time.time()))

    preproc_base = "/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/preproc"
    psf_base = "/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/exposures"

    # Prepare Tasks
    bundle_tasks = []
    if args.bundles_per_arm > 0:
        for band in ['b', 'r', 'z']:
            count = 0
            while count < args.bundles_per_arm:
                night = args.night if args.night else get_random_night_exp()
                if not night: break
                night_dir = os.path.join(preproc_base, night)
                if not os.path.isdir(night_dir): continue
                exps = sorted([e for e in os.listdir(night_dir) if e.isdigit()])
                if not exps: continue
                expid = args.expid if args.expid else random.choice(exps)
                cam = f"{band}{random.randint(0,9)}"
                
                arc_file = os.path.join(preproc_base, night, expid, f"preproc-{cam}-{expid}.fits.gz")
                in_psf = os.path.join(psf_base, night, expid, f"shifted-input-psf-{cam}-{expid}.fits")
                
                if os.path.exists(arc_file) and os.path.exists(in_psf):
                    bundle_tasks.append((night, expid, cam, random.randint(0,19), len(bundle_tasks) % args.workers, args.outdir))
                    count += 1

    ccd_tasks = []
    if args.ccds_per_arm > 0:
        for band in ['b', 'r', 'z']:
            count = 0
            while count < args.ccds_per_arm:
                night = args.night if args.night else get_random_night_exp()
                if not night: break
                night_dir = os.path.join(preproc_base, night)
                if not os.path.isdir(night_dir): continue
                exps = sorted([e for e in os.listdir(night_dir) if e.isdigit()])
                if not exps: continue
                expid = args.expid if args.expid else random.choice(exps)
                cam = f"{band}{random.randint(0,9)}"
                
                arc_file = os.path.join(preproc_base, night, expid, f"preproc-{cam}-{expid}.fits.gz")
                in_psf = os.path.join(psf_base, night, expid, f"shifted-input-psf-{cam}-{expid}.fits")
                
                if os.path.exists(arc_file) and os.path.exists(in_psf):
                    ccd_tasks.append((night, expid, cam, len(ccd_tasks) % args.workers, args.outdir))
                    count += 1

    # Run Bundle Tests
    results_bundle = []
    if bundle_tasks:
        print(f"\n>>> Running {len(bundle_tasks)} Bundle Tests...")
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = {executor.submit(run_bundle_test, task): task for task in bundle_tasks}
            for f in as_completed(futures):
                res = f.result()
                if res:
                    results_bundle.append(res)
                    print(f"CHECKPOINT BUNDLE: {res['cam']} B{res['bundle']} RMS={res['x_rms']:.6f}")

    # Run CCD Tests
    results_ccd = []
    if ccd_tasks:
        print(f"\n>>> Running {len(ccd_tasks)} Full CCD Tests...")
        with ProcessPoolExecutor(max_workers=min(args.workers, len(ccd_tasks))) as executor:
            futures = {executor.submit(run_ccd_test, task): task for task in ccd_tasks}
            for f in as_completed(futures):
                res = f.result()
                if res:
                    results_ccd.append(res)
                    print(f"CHECKPOINT CCD: {res['cam']} RMS={res['x_rms']:.6f}")

    print(f"\n>>> Summary written to {os.path.join(args.outdir, 'validation_summary.txt')}")

if __name__ == "__main__":
    main()
