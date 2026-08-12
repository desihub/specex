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

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "testing"))
from select_test_case import parse_log_line

SCRIPTS_DIR = "/global/cfs/cdirs/desi/spectro/redux/matterhorn/run/scripts/night"
PRODUCTION_REDUX_ROOT = "/global/cfs/cdirs/desi/spectro/redux/matterhorn"
ALL_CAMERAS = [f"{b}{s}" for b in "brz" for s in range(10)]

# Validated per-band worker packing (testing/07Aug2026-30ccd-campaign/pinned30
# and this session's 2-node follow-up) -- z's larger per-fiber design matrix
# (trace-per-fiber-deg=6, the production default) OOMs above this; r silently
# drops bundles (RESOURCE_EXHAUSTED) above 7.
DEFAULT_WORKERS_PER_GPU = {"b": 10, "r": 7, "z": 4}


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


def run_camera_python(cam, case, gpu_id, outdir, wpg_override, dry_run):
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


def run_node_python(cameras, cases, n_gpus, outdir, wpg_override, dry_run, results, results_lock):
    work_q = queue.Queue()
    for cam in cameras:
        work_q.put(cam)

    def slot_worker(gpu_id):
        while True:
            try:
                cam = work_q.get_nowait()
            except queue.Empty:
                return
            cam, dt, rc, n_bf, err_tail = run_camera_python(cam, cases[cam], gpu_id, outdir, wpg_override, dry_run)
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
        run_node_python(node_splits[0], cases, n_gpus, outdir, wpg_override, args.dry_run, results, results_lock)
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
        run_node_python(cameras, cases, args.gpus_per_node or detect_gpus_per_node(),
                         args.outdir, wpg_override, args.dry_run, results, lock)
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
