#!/usr/bin/env python3
"""Single entry point to fit all cameras of a night/expid, switchable
between the C++ production pipeline and the Python/JAX GPU port via one
flag (or the SPECEX_BACKEND env var) -- see how-to-run.md Section 4.3.

    python testing/run_night.py --night 20260401 --expid 00344649 --backend python
    python testing/run_night.py --night 20260401 --expid 00344649 --backend cpp
    SPECEX_BACKEND=python python testing/run_night.py --night 20260401 --expid 00344649

Scope: this fits an ALREADY-PREPROCESSED exposure (preproc-*.fits.gz +
shifted-input-psf-*.fits must already exist -- true for any real matterhorn
production night/expid). It does not run raw-data preprocessing itself.

  --backend cpp    : runs the real production driver, `desi_proc --mpi`,
                      which does its own preprocessing (idempotent -- skips
                      it if outputs already exist) then calls the C++
                      desi_psf_fit binary per camera via MPI ranks. Single
                      command, one srun call, scales via --nodes.
  --backend python : runs this project's `python -m specex.specex` once per
                      camera, pinned to a dedicated GPU (no GPU sharing
                      across cameras), using a per-node dynamic work queue
                      so however many GPUs you have stay busy -- the scheme
                      validated (correctness + timing) against real C++
                      production output across two independent nights as of
                      2026-08-10 (13.2min/30cam on 1 node/4 GPU; 5.92-
                      6.7min/30cam on 2 nodes/8 GPU). Assumes preprocessing
                      already done (see Scope above) -- this backend does
                      NOT call desi_proc at all.
  --backend cpp-direct : runs the real C++ `desi_compute_psf --mpi` binary
                      once per camera (--cpp-ranks MPI ranks each, default
                      20, matching full_ccd_campaign.py's proven per-camera
                      invocation), reading the same already-preprocessed
                      inputs as --backend python and writing the same
                      fit-psf-<cam>-<expid>.fits/.log naming. Does NOT call
                      desi_proc -- no idempotent-preprocessing pass, no
                      whole-job MPI collective for 101+ ranks to hang on if
                      one camera's input is missing (see how-to-run.md
                      Section 7). Sequential by default (--cpp-concurrency);
                      CPU-only, safe to run alongside --backend python on
                      the same node (disjoint resources -- see
                      full_ccd_campaign.py's run_cpp_and_py_concurrent) but
                      intended for a dedicated CPU node for clean timing.

Both backends default to fitting the standard 30 cameras (b0-9, r0-9,
z0-9); pass --cameras to restrict.
"""
import os
import sys
import glob
import re
import json
import time
import argparse
import subprocess
import threading
import queue
import multiprocessing as mp
import traceback

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "testing"))
from select_test_case import parse_log_line

SCRIPTS_DIR = "/global/cfs/cdirs/desi/spectro/redux/matterhorn/run/scripts/night"
PRODUCTION_REDUX_ROOT = "/global/cfs/cdirs/desi/spectro/redux/matterhorn"
ALL_CAMERAS = [f"{b}{s}" for b in "brz" for s in range(10)]
LAMP_LINES_FILE = os.path.join(REPO, "py", "specex", "data", "specex_linelist_desi.txt")

# Validated per-band worker packing. Originally tuned lower (b10/r7/z4,
# testing/07Aug2026-30ccd-campaign/pinned30 -- subprocess-mode, one Pool
# spawned fresh per camera) for z's larger per-fiber design matrix
# (trace-per-fiber-deg=6, the production default) OOMing above that, and r
# silently dropping bundles (RESOURCE_EXHAUSTED) above 7. Bumped to these
# values after a 4-night persistent-mode campaign (2026-09-02, --worker-mode
# persistent + pool-reuse, single 40GB-A100 node) found 0 bundle failures /
# 0 OOM-retries at b12/r8/z5 across all 4 nights, for an additional ~11.5%
# mean wall-time win on top of pool-reuse's own ~6.7% -- reflects headroom
# now available with the persistent worker's post-pool MEM_FRACTION=0.05
# cap. Untested above these values and untested on subprocess mode -- pass
# --workers-per-gpu-{b,r,z} to override (e.g. back to 10/7/4) if a future
# OOM/regression surfaces here.
DEFAULT_WORKERS_PER_GPU = {"b": 12, "r": 8, "z": 5}


def find_cases(night, expid, cameras):
    """Single pass over the night's arc*.log files -- looking each camera up
    individually via full_ccd_campaign.find_case() is O(cameras x logs) and
    was measured to take >120s for 30 cameras; this is one pass total."""
    log_dir = os.path.join(SCRIPTS_DIR, night)
    wanted = set(cameras)
    found = {}
    for log in glob.glob(os.path.join(log_dir, "arc*.log")):
        with open(log) as f:
            for line in f:
                if "desi_compute_psf" not in line or expid not in line:
                    continue
                case = parse_log_line(line)
                if case and case.get("camera") in wanted:
                    found[case["camera"]] = case
        if len(found) == len(wanted):
            break
    missing = wanted - set(found)
    if missing:
        print(f"WARNING: no arc log entry found for {sorted(missing)} on {night}/{expid} -- skipping", file=sys.stderr)
    return found


def slurm_nodes():
    """Hostnames of this job's allocation, in a stable order. Falls back to
    [None] (meaning "run locally, no srun -w pinning") outside a SLURM job."""
    nodelist = os.environ.get("SLURM_JOB_NODELIST") or os.environ.get("SLURM_NODELIST")
    if not nodelist:
        return [None]
    out = subprocess.run(["scontrol", "show", "hostnames", nodelist], capture_output=True, text=True)
    hosts = [h for h in out.stdout.split() if h]
    return hosts or [None]


def resolve_python_run_dir(args):
    """Where --backend python or cpp-direct writes fit-psf-<cam>-<expid>.fits
    + per-camera logs for this night/expid. Three tiers, in order:
      1. --outdir, if given -- used exactly as-is (a flat directory, for
         ad-hoc/scratch runs where production-style nesting isn't wanted).
      2. $DESI_SPECTRO_REDUX (+ $SPECPROD, default $USER) if set --
         mirrors --backend cpp's/desi_proc's own exposures/<night>/<expid>/
         layout, so `export DESI_SPECTRO_REDUX=...` before running either
         backend now lands both in the same, directly comparable location
         instead of being silently ignored by --backend python.
      3. $SCRATCH/specex/run_night_<night>_<expid>/ (previous default).
    Refuses to resolve inside the real production tree -- never write
    fit-psf output there by accident, explicit override or not.
    """
    if args.outdir is not None:
        run_dir = args.outdir
    elif os.environ.get("DESI_SPECTRO_REDUX"):
        specprod = os.environ.get("SPECPROD") or os.environ.get("USER", "specex")
        run_dir = os.path.join(os.environ["DESI_SPECTRO_REDUX"], specprod,
                                "exposures", args.night, args.expid)
    else:
        scratch = os.environ.get("SCRATCH", "/tmp")
        run_dir = os.path.join(scratch, "specex", f"run_night_{args.night}_{args.expid}")

    real = os.path.realpath(run_dir)
    if real == PRODUCTION_REDUX_ROOT or real.startswith(PRODUCTION_REDUX_ROOT + os.sep):
        print(f"ERROR: resolved output dir {run_dir} is inside the real production "
              f"redux tree ({PRODUCTION_REDUX_ROOT}) -- refusing to write there. "
              f"Pass --outdir or point $DESI_SPECTRO_REDUX at a private tree.", file=sys.stderr)
        sys.exit(1)
    return run_dir


def detect_gpus_per_node(default=4):
    try:
        out = subprocess.run(["nvidia-smi", "-L"], capture_output=True, text=True, timeout=15)
        n = len([l for l in out.stdout.splitlines() if l.strip()])
        return n if n > 0 else default
    except Exception:
        return default


# ---------------------------------------------------------------------------
# --backend python
# ---------------------------------------------------------------------------

# Coarse per-band relative duration hint, used only to ORDER a node's own
# queue when no real --lpt-profile is available -- not meant to be an
# accurate absolute estimate, just enough to stop z (the heaviest, most
# variable band) from being queued dead last. Roughly matches typical
# wpg=10 b/7 r/4 z per-camera times (see porting-notes.md).
DEFAULT_DURATION_HINT = {"b": 90.0, "r": 100.0, "z": 150.0}


def split_cameras_for_nodes(cameras, n_nodes, lpt_profile):
    """LPT (longest-processing-time-first) balanced split across n_nodes if
    a per-camera timing profile is given (JSON: {"b0": 69.2, ...}, from a
    prior run's own measured wall times -- see how-to-run.md Section 4.3);
    otherwise a naive alternating split (band-diverse but not load-balanced,
    since we have no timing prior for an arbitrary fresh night/expid).

    Whichever way cameras get ASSIGNED to a node, each node's own returned
    list is separately sorted by descending expected duration before being
    handed to run_node_python's FIFO queue -- assignment and per-node queue
    ORDER are different problems. `cameras` comes in alphabetically sorted
    (b0..b9, r0..r9, z0..z9); left as-is, that queues every b before every
    r before every z, so by the time a node's 4 GPU slots reach z -- the
    heaviest, most variable band -- there's no lighter b/r work left to
    overlap a straggling z camera against, and the whole node's wall time
    is gated by z's own tail. This is the exact single-node instance of the
    cross-node tail-starvation pattern porting-notes.md's 2026-08-10 LPT
    session found and fixed at the multi-node level; fixing it here too
    (uses the real profile if given, else DEFAULT_DURATION_HINT) closes it
    for the single-node case as well, where it was previously untouched --
    n_nodes==1 didn't even look at --lpt-profile before this."""
    if n_nodes == 1:
        bins = [list(cameras)]
    elif lpt_profile:
        durations = {c: lpt_profile.get(c, lpt_profile.get("__default__", DEFAULT_DURATION_HINT[c[0]])) for c in cameras}
        order = sorted(cameras, key=lambda c: -durations[c])
        bins = [[] for _ in range(n_nodes)]
        load = [0.0] * n_nodes
        for cam in order:
            i = min(range(n_nodes), key=lambda k: load[k])
            bins[i].append(cam)
            load[i] += durations[cam]
    else:
        bins = [cameras[i::n_nodes] for i in range(n_nodes)]

    durations = {c: (lpt_profile.get(c, lpt_profile.get("__default__", DEFAULT_DURATION_HINT[c[0]]))
                      if lpt_profile else DEFAULT_DURATION_HINT[c[0]]) for c in cameras}
    return [sorted(b, key=lambda c: -durations[c]) for b in bins]


def tail_error(log_path, n=6):
    """Last few non-blank lines of a camera's log, for surfacing *why* a
    failure happened in the run summary -- rather than the C++ path's
    failure mode observed this session, where 2 of 101 MPI ranks failing
    fast (missing preproc input) left the whole job hung with zero further
    log output for 36+ minutes and no summary at all. Each camera here is
    an independent subprocess, so one failing can never block the others,
    but a bare 'rc=1' with no context is still not a *clean report* --
    this closes that gap."""
    try:
        with open(log_path) as f:
            lines = [l.rstrip() for l in f if l.strip()]
        return lines[-n:]
    except OSError:
        return []


def run_camera_python(cam, case, gpu_id, outdir, wpg_override, dry_run, footprint_margin=None):
    band = cam[0]
    wpg = wpg_override.get(band, DEFAULT_WORKERS_PER_GPU[band])
    # fit-psf-<cam>-<expid>.fits matches desi_proc's/desi_compute_psf's own
    # real output naming (see run_backend_cpp's printed output path below)
    # -- was bare "{cam}.fits", which didn't line up with anything C++
    # produces and made the two backends' outputs hard to tell apart or
    # script against generically.
    out_fits = os.path.join(outdir, f"fit-psf-{cam}-{case['expid']}.fits")
    # Same stem as out_fits, .log instead of .fits -- was "py-<cam>.log"
    # (no expid, so it collided across different expids of the same camera
    # sharing one outdir, and didn't visually pair with its own fits file).
    log = os.path.join(outdir, f"fit-psf-{cam}-{case['expid']}.log")

    # Check inputs up front rather than letting the subprocess fail deep
    # inside fitsio/io.py with a traceback that reads like a code bug --
    # a missing preproc/shifted-input-psf file is an input-data problem,
    # not a specex problem, and should be reported as one immediately.
    missing = [p for p in (case["image"], case["input_psf"]) if not os.path.exists(p)]
    if missing:
        with open(log, "w") as f:
            f.write(f"SKIPPED: missing input file(s): {missing}\n")
        return cam, 0.0, "SKIPPED", 0, [f"missing input file(s): {missing}"]

    cmd = [sys.executable, "-m", "specex.specex",
           "-a", case["image"], "--in-psf", case["input_psf"],
           "--out-psf", out_fits, "--gpu", "1", "--workers-per-gpu", str(wpg)]
    if case.get("broken_fibers"):
        cmd += ["--broken-fibers", case["broken_fibers"]]
    if footprint_margin is not None:
        cmd += ["--footprint-margin", str(footprint_margin)]
    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    if dry_run:
        print(f"  [DRY RUN] CUDA_VISIBLE_DEVICES={gpu_id} {' '.join(cmd)}")
        return cam, 0.0, 0, 0, []
    t0 = time.time()
    with open(log, "w") as f:
        rc = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, env=env).returncode
    dt = time.time() - t0
    with open(log) as f:
        n_bundle_fail = sum(1 for line in f if "WARNING: Bundle" in line)
    err_tail = tail_error(log) if rc != 0 else []
    return cam, dt, rc, n_bundle_fail, err_tail


def run_node_python(cameras, cases, n_gpus, outdir, wpg_override, dry_run, results, results_lock, footprint_margin=None):
    work_q = queue.Queue()
    for cam in cameras:
        work_q.put(cam)

    def slot_worker(gpu_id):
        while True:
            try:
                cam = work_q.get_nowait()
            except queue.Empty:
                return
            cam, dt, rc, n_bf, err_tail = run_camera_python(cam, cases[cam], gpu_id, outdir, wpg_override, dry_run, footprint_margin)
            with results_lock:
                results.append((cam, dt, rc, n_bf, err_tail))
            flag = f" *** {n_bf} BUNDLE FAILURES ***" if n_bf else ""
            if rc == "SKIPPED":
                print(f"[{time.strftime('%H:%M:%S')}] {os.uname().nodename} gpu{gpu_id}: {cam} SKIPPED -- {err_tail[0]}", flush=True)
            else:
                print(f"[{time.strftime('%H:%M:%S')}] {os.uname().nodename} gpu{gpu_id}: {cam} done in {dt:.1f}s rc={rc}{flag}", flush=True)
                if rc != 0:
                    for line in err_tail:
                        print(f"    | {line}", flush=True)

    threads = [threading.Thread(target=slot_worker, args=(g,)) for g in range(n_gpus)]
    [t.start() for t in threads]
    [t.join() for t in threads]


def _fit_one_camera_inprocess(cam, case, out_fits, log_path, wpg, footprint_margin, bundle_pool=None):
    """Fit one camera by calling fit_ccd_native() directly, instead of
    shelling out to a fresh `python -m specex.specex` subprocess -- the
    persistent-worker path's actual payload. Args mirror specex.py's own
    main()/argparse defaults exactly (trace_per_fiber_deg=6, trace_prior_deg=1,
    trace_prior_weight=1e5, trace_prior_ndead_threshold=500,
    masked_amp_ndead_threshold=8000, max_number_of_lines=200 -- note this is
    main()'s --max-lines default, NOT fit_ccd_native's own function-signature
    default of 100 -- footprint_margin from the caller, everything else auto)
    so a persistent-worker run is bit-for-bit equivalent to the existing
    subprocess path, not just "close enough". Import of specex.specex is
    deferred to inside this call (not module level) so CUDA_VISIBLE_DEVICES,
    set by the caller before the worker's first camera, is honored -- jax
    latches its visible-device set at first import.

    Returns: (rc, n_bundle_fail) -- rc=0 success, rc=1 failure/exception,
    rc="SKIPPED" for missing inputs (matching run_camera_python's contract).
    """
    missing = [p for p in (case["image"], case["input_psf"]) if not os.path.exists(p)]
    if missing:
        with open(log_path, "w") as f:
            f.write(f"SKIPPED: missing input file(s): {missing}\n")
        return "SKIPPED", 0

    # Redirect at the OS file-descriptor level (not just sys.stdout) so
    # fit_ccd_native's own internal spawn-context bundle-worker pool --
    # separate processes that inherit fds at spawn time -- also lands in
    # this camera's log file, matching what subprocess.run(stdout=f) gave
    # the old per-camera-subprocess path for free.
    from specex.specex import fit_ccd_native
    old_out, old_err = os.dup(1), os.dup(2)
    rc = 0
    with open(log_path, "w") as f:
        os.dup2(f.fileno(), 1)
        os.dup2(f.fileno(), 2)
        try:
            failed_bundles = fit_ccd_native(
                arc_file=case["image"], in_psf_file=case["input_psf"], out_psf_file=out_fits,
                lamp_lines_file=LAMP_LINES_FILE, n_gpus=1, backend="gpu",
                broken_fibers=case.get("broken_fibers"), sn_threshold=3.0,
                max_number_of_lines=200, h_size_y=5, workers_per_gpu=wpg,
                trace_per_fiber_deg=6, trace_prior_deg=1, trace_prior_weight=1e5,
                trace_prior_ndead_threshold=500, masked_amp_ndead_threshold=8000,
                footprint_margin=footprint_margin, bundle_pool=bundle_pool,
                bundle_log_path=(log_path if bundle_pool is not None else None),
            )
            if failed_bundles:
                rc = 1
        except Exception:
            traceback.print_exc()
            rc = 1
        finally:
            sys.stdout.flush(); sys.stderr.flush()
            os.dup2(old_out, 1); os.dup2(old_err, 2)
            os.close(old_out); os.close(old_err)

    with open(log_path) as f:
        n_bundle_fail = sum(1 for line in f if "WARNING: Bundle" in line)
    return rc, n_bundle_fail


def _gpu_persistent_worker(gpu_id, work_q, results_q, outdir, wpg_override, footprint_margin, pool_reuse=True):
    """Process target for the persistent-worker path: one long-lived process
    per GPU that imports jax/specex ONCE (paying the interpreter-startup +
    module-import cost a single time for the whole night) then pulls
    cameras off the shared queue until it sees the None sentinel -- as
    opposed to run_camera_python's design, which pays that cost fresh for
    every camera via a brand-new subprocess. Measured on 20260401/00344649:
    ~7-10s of pure launch overhead per camera regardless of band (outer
    subprocess wall time minus fit_ccd_native's own internally-timed "Total
    CCD Fit Time"), serialized ~7.5x per GPU lane on a 30-camera/4-GPU
    night -- this is the thing that overhead was going to.

    CUDA_VISIBLE_DEVICES must be set before the first `import jax` in this
    process (here, transitively via the first _fit_one_camera_inprocess
    call's `from specex.specex import fit_ccd_native`) -- set eagerly at
    worker startup rather than relying on that deferred import.

    XLA_PYTHON_CLIENT_MEM_FRACTION is capped low for the SAME reason: this
    persistent process's own JAX usage isn't limited to the per-camera
    bundle-worker child pool (which fit_bundle_task already isolates with
    its own explicit env overrides, see specex.py) -- fit_ccd_native's
    post-pool merge/write_python_psf() step runs real jax.numpy ops
    directly in THIS process, on THIS same physical GPU. Confirmed by a
    first real test on 20260401/00344649 (8 cameras, 2/GPU): every worker's
    FIRST camera succeeded cleanly but its SECOND hit partial
    RESOURCE_EXHAUSTED bundle failures (5-10 of 20 bundles) -- consistent
    with JAX's BFC allocator growing to ~75% of *whatever's still free* on
    this process's first GPU touch and never shrinking back, permanently
    starving every subsequent camera's freshly-spawned bundle pool on the
    same GPU. Capped at 5% here; fit_bundle_task's own children explicitly
    re-pin to 75% for their own process regardless of what they inherit
    from this parent's environment, so real per-bundle compute is
    unaffected."""
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
    os.environ["XLA_PYTHON_CLIENT_MEM_FRACTION"] = "0.05"

    # Reuse one bundle-worker pool across cameras instead of paying a fresh
    # spawn+JAX-import cost (measured far larger than this outer process's
    # own per-camera overhead -- fit_ccd_native's internal _run_batch
    # otherwise creates and tears down a brand-new Pool of `wpg` processes
    # on EVERY camera call, even in persistent-worker mode) -- only recreated
    # on a band transition, since workers_per_gpu differs by band (see
    # DEFAULT_WORKERS_PER_GPU above, tuned to each band's GPU-memory
    # footprint) and the shared work queue is filled band-sorted (b's, then
    # r's, then z's), so a single GPU worker sees at most ~2 real
    # transitions per night, not one per camera. Confirmed itself worth
    # ~6.7% mean wall time on a 4-night persistent-mode campaign
    # (2026-09-02) vs a --no-pool-reuse control on the same nights/node,
    # stacking with (not substituting for) the DEFAULT_WORKERS_PER_GPU bump
    # above (~11.5% more) for ~17.5% combined. --no-pool-reuse still exists
    # to revert to a fresh Pool per camera for future A/B comparisons.
    # OOM retries inside fit_ccd_native still fall back to
    # their own temporary reduced-packing pool (see specex.py's _run_batch)
    # regardless of this pool -- unaffected, still correctness-preserving.
    ctx = mp.get_context('spawn')
    pool = None
    pool_size = None
    try:
        while True:
            item = work_q.get()
            if item is None:
                return
            cam, case = item
            band = cam[0]
            wpg = wpg_override.get(band, DEFAULT_WORKERS_PER_GPU[band])
            if not pool_reuse:
                # A/B-test path: tear down and recreate every camera,
                # matching pre-pool-reuse behavior exactly (still goes
                # through the bundle_pool= plumbing/log-redirect, just
                # never actually reused across calls).
                if pool is not None:
                    pool.close()
                    pool.join()
                pool = ctx.Pool(processes=wpg)
                pool_size = wpg
            elif pool is None or pool_size != wpg:
                if pool is not None:
                    pool.close()
                    pool.join()
                pool = ctx.Pool(processes=wpg)
                pool_size = wpg
            out_fits = os.path.join(outdir, f"fit-psf-{cam}-{case['expid']}.fits")
            log_path = os.path.join(outdir, f"fit-psf-{cam}-{case['expid']}.log")
            t0 = time.time()
            rc, n_bf = _fit_one_camera_inprocess(cam, case, out_fits, log_path, wpg, footprint_margin, bundle_pool=pool)
            dt = time.time() - t0
            err_tail = tail_error(log_path) if rc != 0 else []
            results_q.put((cam, dt, rc, n_bf, err_tail, gpu_id))
    finally:
        if pool is not None:
            pool.close()
            pool.join()


def run_node_python_persistent(cameras, cases, n_gpus, outdir, wpg_override, dry_run, results, results_lock, footprint_margin=None, pool_reuse=True):
    """Persistent-worker counterpart to run_node_python: same shared-queue
    dynamic dispatch (whichever GPU frees up next claims the next camera),
    but n_gpus long-lived worker PROCESSES instead of n_gpus threads each
    shelling out to a fresh subprocess per camera. A worker process dying
    outright (e.g. a CUDA-driver-level crash, not a catchable Python
    exception -- those are already caught inside _fit_one_camera_inprocess)
    only loses whatever single camera it was actively running: every other
    camera stays in the shared queue for a still-alive worker to pick up,
    so the blast radius is one camera, not the whole night -- preserving
    the crash-isolation property the subprocess-per-camera design had,
    just at finer granularity than "the rest of this worker's queue"."""
    if dry_run:
        for cam in cameras:
            case = cases[cam]
            out_fits = os.path.join(outdir, f"fit-psf-{cam}-{case['expid']}.fits")
            fm = footprint_margin if footprint_margin is not None else 7
            print(f"  [DRY RUN] (persistent) {cam}: fit_ccd_native(arc={case['image']}, "
                  f"in_psf={case['input_psf']}, out_psf={out_fits}, footprint_margin={fm}) [in-process]")
            results.append((cam, 0.0, 0, 0, []))
        return

    ctx = mp.get_context("spawn")
    work_q, results_q = ctx.Queue(), ctx.Queue()
    for cam in cameras:
        work_q.put((cam, cases[cam]))
    for _ in range(n_gpus):
        work_q.put(None)

    workers = [ctx.Process(target=_gpu_persistent_worker, args=(g, work_q, results_q, outdir, wpg_override, footprint_margin, pool_reuse))
               for g in range(n_gpus)]
    [w.start() for w in workers]

    # A caught Python exception inside _fit_one_camera_inprocess still
    # reports normally via results_q (rc=1). But a hard crash (segfault,
    # CUDA driver abort) kills the worker process without ever putting a
    # result -- polling with a timeout, rather than a plain blocking get()
    # for exactly len(cameras) results, means that loses only the one
    # in-flight camera instead of hanging this whole night forever waiting
    # for a result that will never arrive.
    n_done = 0
    seen_cams = set()
    while n_done < len(cameras):
        try:
            cam, dt, rc, n_bf, err_tail, gpu_id = results_q.get(timeout=10)
        except queue.Empty:
            if not any(w.is_alive() for w in workers):
                missing = len(cameras) - n_done
                print(f"WARNING: all persistent workers exited but {missing} camera(s) never reported a result "
                      f"(worker crash) -- treating as failed and stopping this node's collection.", flush=True)
                break
            continue
        seen_cams.add(cam)
        n_done += 1
        with results_lock:
            results.append((cam, dt, rc, n_bf, err_tail))
        flag = f" *** {n_bf} BUNDLE FAILURES ***" if n_bf else ""
        if rc == "SKIPPED":
            print(f"[{time.strftime('%H:%M:%S')}] {os.uname().nodename} gpu{gpu_id}: {cam} SKIPPED -- {err_tail[0]}", flush=True)
        else:
            print(f"[{time.strftime('%H:%M:%S')}] {os.uname().nodename} gpu{gpu_id}: {cam} done in {dt:.1f}s rc={rc}{flag}", flush=True)
            if rc != 0:
                for line in err_tail:
                    print(f"    | {line}", flush=True)
    for cam in cameras:
        if cam not in seen_cams:
            with results_lock:
                results.append((cam, 0.0, 1, 0, ["worker crashed before reporting a result"]))

    [w.join() for w in workers]


def run_backend_python(args, cases):
    cameras = sorted(cases)
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    hosts = slurm_nodes()
    n_nodes = args.nodes or len(hosts)
    if n_nodes > len(hosts):
        print(f"WARNING: --nodes {n_nodes} exceeds this job's {len(hosts)}-node allocation; clamping.", file=sys.stderr)
        n_nodes = len(hosts)
    n_gpus = args.gpus_per_node or detect_gpus_per_node()

    lpt_profile = None
    if args.lpt_profile:
        with open(args.lpt_profile) as f:
            lpt_profile = json.load(f)

    wpg_override = {}
    if args.workers_per_gpu_b is not None: wpg_override["b"] = args.workers_per_gpu_b
    if args.workers_per_gpu_r is not None: wpg_override["r"] = args.workers_per_gpu_r
    if args.workers_per_gpu_z is not None: wpg_override["z"] = args.workers_per_gpu_z

    node_splits = split_cameras_for_nodes(cameras, n_nodes, lpt_profile)
    print(f"=== backend=python nodes={n_nodes} gpus/node={n_gpus} cameras={len(cameras)} "
          f"split={[len(s) for s in node_splits]}{' (LPT-balanced)' if lpt_profile else ' (naive)'} ===", flush=True)

    results = []
    results_lock = threading.Lock()
    t0 = time.time()

    if n_nodes == 1 and hosts[0] is None:
        # not in a SLURM job (or single node with no need to srun -w) -- run
        # directly in this process
        if args.worker_mode == "persistent":
            run_node_python_persistent(node_splits[0], cases, n_gpus, outdir, wpg_override, args.dry_run, results, results_lock, args.footprint_margin, pool_reuse=not args.no_pool_reuse)
        else:
            run_node_python(node_splits[0], cases, n_gpus, outdir, wpg_override, args.dry_run, results, results_lock, args.footprint_margin)
    else:
        procs = []
        for i in range(n_nodes):
            node_cams = node_splits[i]
            if not node_cams:
                continue
            inner = [sys.executable, __file__, "--_node-worker",
                     "--night", args.night, "--expid", args.expid,
                     "--outdir", outdir, "--gpus-per-node", str(n_gpus),
                     "--cameras", ",".join(node_cams)]
            for band, wpg in wpg_override.items():
                inner += [f"--workers-per-gpu-{band}", str(wpg)]
            if args.footprint_margin is not None:
                inner += ["--footprint-margin", str(args.footprint_margin)]
            if args.worker_mode == "persistent":
                inner += ["--worker-mode", "persistent"]
                if args.no_pool_reuse:
                    inner += ["--no-pool-reuse"]
            if args.dry_run:
                inner += ["--dry-run"]
            cmd = ["srun", "-N1", "-n1", "-w", hosts[i]] + inner if hosts[i] else inner
            print(f"  launching node {i} ({hosts[i] or 'local'}): {len(node_cams)} cameras", flush=True)
            procs.append(subprocess.Popen(cmd))
        for p in procs:
            p.wait()
        # _node-worker mode writes its own per-camera logs under outdir;
        # reparse them here for the summary rather than trying to pipe
        # results back across process boundaries.
        for cam in cameras:
            log = os.path.join(outdir, f"fit-psf-{cam}-{args.expid}.log")
            if not os.path.exists(log):
                continue
            with open(log) as f:
                text = f.read()
            n_bf = text.count("WARNING: Bundle")
            ok = "Total CCD Fit Time" in text
            skipped = text.startswith("SKIPPED:")
            rc = "SKIPPED" if skipped else (0 if ok else -1)
            err_tail = [] if (ok and not skipped) else tail_error(log)
            results.append((cam, None, rc, n_bf, err_tail))

    total = time.time() - t0
    n_fail = sum(1 for r in results if r[2] not in (0, None))
    n_bf_cams = sum(1 for r in results if r[3])
    print(f"\n=== SUMMARY backend=python: {len(results)}/{len(cameras)} cameras, "
          f"{n_fail} rc!=0/SKIPPED, {n_bf_cams} cameras with bundle failures, "
          f"TOTAL WALL TIME: {total:.1f}s ({total/60:.1f} min) ===", flush=True)
    print(f"  output: {outdir}/", flush=True)
    for cam, dt, rc, n_bf, err_tail in results:
        if rc not in (0, None) or n_bf:
            print(f"  PROBLEM: {cam} rc={rc} bundle_failures={n_bf}")
            for line in err_tail:
                print(f"    | {line}")


# ---------------------------------------------------------------------------
# --backend cpp
# ---------------------------------------------------------------------------

def run_backend_cpp(args):
    """Shells out to the real production driver, desi_proc --mpi, into a
    private SPECPROD (never the real 'matterhorn' production tree). Rank
    count follows the validated single/multi-node formula from
    porting-notes.md (-N1<->n101 ~11min, -N3<->n301 ~7min -- ranks =
    100*nodes + 1; 600 total bundles / (ranks-1) workers divides evenly for
    any node count on this formula)."""
    n_nodes = args.nodes or 1
    ranks = 100 * n_nodes + 1
    redux_dir = args.redux_dir or os.path.join(args.outdir, "redux")
    os.makedirs(redux_dir, exist_ok=True)

    env = os.environ.copy()
    env["DESI_SPECTRO_REDUX"] = redux_dir
    env["SPECPROD"] = args.specprod

    cmd = ["srun", "-N", str(n_nodes), "-n", str(ranks), "-c", "2", "--cpu-bind=cores",
           "desi_proc", "-n", args.night, "-e", args.expid, "--mpi"]
    if args.cameras:
        cmd += ["--cameras", args.cameras]

    print(f"=== backend=cpp nodes={n_nodes} ranks={ranks} redux={redux_dir} SPECPROD={args.specprod} ===", flush=True)
    print(f"  {' '.join(cmd)}", flush=True)
    if args.dry_run:
        return

    log_path = os.path.join(args.outdir, f"cpp_{args.night}_{args.expid}.log")
    os.makedirs(args.outdir, exist_ok=True)
    t0 = time.time()
    with open(log_path, "w") as f:
        rc = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, env=env).returncode
    total = time.time() - t0
    print(f"\n=== SUMMARY backend=cpp: rc={rc} TOTAL WALL TIME: {total:.1f}s ({total/60:.1f} min) "
          f"(includes idempotent preprocessing -- see log for a specex-only breakdown) ===", flush=True)
    print(f"  log: {log_path}", flush=True)
    print(f"  output: {redux_dir}/{args.specprod}/exposures/{args.night}/{args.expid}/", flush=True)


# ---------------------------------------------------------------------------
# --backend cpp-direct
# ---------------------------------------------------------------------------

def run_camera_cpp(cam, case, outdir, ranks, dry_run):
    """One `desi_compute_psf --mpi` call for a single camera, same pattern
    as full_ccd_campaign.py's start_cpp_full but writing this project's
    fit-psf-<cam>-<expid> naming so cpp-direct and python outputs sit in
    the same outdir, directly comparable by filename alone."""
    tag = f"fit-psf-{cam}-{case['expid']}"
    out_fits = os.path.join(outdir, f"{tag}.fits")
    log = os.path.join(outdir, f"{tag}.log")

    missing = [p for p in (case["image"], case["input_psf"]) if not os.path.exists(p)]
    if missing:
        with open(log, "w") as f:
            f.write(f"SKIPPED: missing input file(s): {missing}\n")
        return cam, 0.0, "SKIPPED", [f"missing input file(s): {missing}"]

    cmd = ["srun", "-n", str(ranks), "--cpu-bind=cores", "desi_compute_psf", "--mpi",
           "--input-image", case["image"], "--input-psf", case["input_psf"],
           "-o", out_fits]
    if case.get("broken_fibers"):
        cmd += ["--broken-fibers", case["broken_fibers"]]
    if dry_run:
        print(f"  [DRY RUN] {' '.join(cmd)}")
        return cam, 0.0, 0, []
    t0 = time.time()
    with open(log, "w") as f:
        rc = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT).returncode
    dt = time.time() - t0
    err_tail = tail_error(log) if rc != 0 else []
    return cam, dt, rc, err_tail


def run_backend_cpp_direct(args, cases):
    cameras = sorted(cases)
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    concurrency = args.cpp_concurrency

    print(f"=== backend=cpp-direct cameras={len(cameras)} ranks/camera={args.cpp_ranks} "
          f"concurrency={concurrency} (no desi_proc) ===", flush=True)

    work_q = queue.Queue()
    for cam in cameras:
        work_q.put(cam)
    results = []
    results_lock = threading.Lock()
    t0 = time.time()

    def worker():
        while True:
            try:
                cam = work_q.get_nowait()
            except queue.Empty:
                return
            cam, dt, rc, err_tail = run_camera_cpp(cam, cases[cam], outdir, args.cpp_ranks, args.dry_run)
            with results_lock:
                results.append((cam, dt, rc, err_tail))
            if rc == "SKIPPED":
                print(f"[{time.strftime('%H:%M:%S')}] {cam} SKIPPED -- {err_tail[0]}", flush=True)
            else:
                print(f"[{time.strftime('%H:%M:%S')}] {cam} done in {dt:.1f}s rc={rc}", flush=True)
                if rc != 0:
                    for line in err_tail:
                        print(f"    | {line}", flush=True)

    threads = [threading.Thread(target=worker) for _ in range(concurrency)]
    [t.start() for t in threads]
    [t.join() for t in threads]

    total = time.time() - t0
    n_fail = sum(1 for r in results if r[2] not in (0, None))
    print(f"\n=== SUMMARY backend=cpp-direct: {len(results)}/{len(cameras)} cameras, "
          f"{n_fail} rc!=0/SKIPPED, TOTAL WALL TIME: {total:.1f}s ({total/60:.1f} min) ===", flush=True)
    print(f"  output: {outdir}/", flush=True)
    for cam, dt, rc, err_tail in results:
        if rc not in (0, None):
            print(f"  PROBLEM: {cam} rc={rc}")
            for line in err_tail:
                print(f"    | {line}")


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--night", required=True)
    ap.add_argument("--expid", required=True)
    ap.add_argument("--backend", choices=["cpp", "cpp-direct", "python"], default=os.environ.get("SPECEX_BACKEND", "python"),
                     help="Default: $SPECEX_BACKEND env var, or 'python' if unset.")
    ap.add_argument("--cameras", help="Comma-separated camera list (default: all 30 standard b0-9/r0-9/z0-9)")
    ap.add_argument("--outdir", default=None, help="Default (--backend python or cpp-direct): $DESI_SPECTRO_REDUX/$SPECPROD/exposures/<night>/<expid> if $DESI_SPECTRO_REDUX is set (matches desi_proc's own layout), else $SCRATCH/specex/run_night_<night>_<expid>. Default (--backend cpp): always $SCRATCH/specex/run_night_<night>_<expid> (this backend manages its own private redux tree under it -- see --redux-dir).")
    ap.add_argument("--nodes", type=int, default=None, help="Default: this SLURM job's full node allocation")
    ap.add_argument("--dry-run", action="store_true", help="Print planned commands without executing them")
    # --backend python only
    ap.add_argument("--gpus-per-node", type=int, default=None, help="Default: auto-detect via nvidia-smi")
    ap.add_argument("--lpt-profile", help="JSON file of {camera: seconds} from a prior run, for LPT-balanced multi-node splitting (see how-to-run.md). Default: naive alternating split.")
    ap.add_argument("--workers-per-gpu-b", type=int, default=None)
    ap.add_argument("--workers-per-gpu-r", type=int, default=None)
    ap.add_argument("--workers-per-gpu-z", type=int, default=None)
    ap.add_argument("--footprint-margin", type=int, default=None, help="Passed through to `python -m specex.specex --footprint-margin` for every camera. Default: None (specex.specex's own default, currently 7). Pass 0 to reproduce pre-fix zero-margin behavior for comparison reruns.")
    ap.add_argument("--worker-mode", choices=["subprocess", "persistent"], default="subprocess", help="'subprocess' (default): fresh `python -m specex.specex` process per camera, matching every prior campaign's methodology exactly. 'persistent': one long-lived worker process per GPU calling fit_ccd_native() in-process for a stream of cameras, avoiding ~7-10s of per-camera interpreter/JAX-import overhead measured on 20260401/00344649 (~65s/night on a 4-GPU node) -- experimental, not yet validated at the same scale as 'subprocess'.")
    ap.add_argument("--no-pool-reuse", action="store_true", help="Persistent mode only: disable bundle-worker Pool reuse across cameras, reverting to a fresh Pool per camera (pre-pool-reuse behavior) -- for A/B timing comparisons only, no correctness effect either way.")
    # internal, used to re-invoke this script once per node via srun
    ap.add_argument("--_node-worker", action="store_true", help=argparse.SUPPRESS)
    # --backend cpp only
    ap.add_argument("--redux-dir", default=None, help="Default: <outdir>/redux (a private tree, never the real matterhorn production output)")
    ap.add_argument("--specprod", default=os.environ.get("USER", "specex"))
    # --backend cpp-direct only
    ap.add_argument("--cpp-ranks", type=int, default=20, help="MPI ranks per desi_compute_psf call (default: 20, matches full_ccd_campaign.py)")
    ap.add_argument("--cpp-concurrency", type=int, default=1, help="How many cameras' desi_compute_psf calls to run at once (default: 1, sequential -- for clean per-camera timing; raise only if you have rank budget to spare)")

    args = ap.parse_args()
    cameras = args.cameras.split(",") if args.cameras else ALL_CAMERAS

    if args._node_worker:
        # invoked by run_backend_python via srun, one call per node --
        # the parent always passes --outdir explicitly (it resolved it
        # once, below, before launching any node), so no re-resolution
        # needed/wanted here.
        wpg_override = {}
        if args.workers_per_gpu_b is not None: wpg_override["b"] = args.workers_per_gpu_b
        if args.workers_per_gpu_r is not None: wpg_override["r"] = args.workers_per_gpu_r
        if args.workers_per_gpu_z is not None: wpg_override["z"] = args.workers_per_gpu_z
        cases = find_cases(args.night, args.expid, cameras)
        results, lock = [], threading.Lock()
        if args.worker_mode == "persistent":
            run_node_python_persistent(cameras, cases, args.gpus_per_node or detect_gpus_per_node(),
                   args.outdir, wpg_override, args.dry_run, results, lock, args.footprint_margin, pool_reuse=not args.no_pool_reuse)
        else:
            run_node_python(cameras, cases, args.gpus_per_node or detect_gpus_per_node(),
                   args.outdir, wpg_override, args.dry_run, results, lock, args.footprint_margin)
        return

    if args.backend == "cpp":
        # Deliberately does NOT honor a pre-set $DESI_SPECTRO_REDUX for its
        # own outdir default -- run_backend_cpp always manages its own
        # private redux tree (see --redux-dir) so a real production
        # DESI_SPECTRO_REDUX left set in the shell can never get written to.
        if args.outdir is None:
            scratch = os.environ.get("SCRATCH", "/tmp")
            args.outdir = os.path.join(scratch, "specex", f"run_night_{args.night}_{args.expid}")
        run_backend_cpp(args)
        return

    args.outdir = resolve_python_run_dir(args)
    cases = find_cases(args.night, args.expid, cameras)
    if not cases:
        print("No cases found -- nothing to do.", file=sys.stderr)
        sys.exit(1)
    os.makedirs(args.outdir, exist_ok=True)
    if args.backend == "cpp-direct":
        run_backend_cpp_direct(args, cases)
    else:
        run_backend_python(args, cases)


if __name__ == "__main__":
    main()
