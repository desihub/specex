import os
import sys
import time
import numpy as np
import multiprocessing as mp

from .io import load_python_psf, read_preproc, read_lamp_lines, write_python_psf
from .fitter import PSF_Fitter, get_bundle_spots, select_bundle_spots_iterative

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
def fit_bundle_task(bid, gpu_id, arc_file, in_psf_file, out_psf_file, lamp_lines_file, backend="gpu", broken_fibers=None, sn_threshold=3.0, h_size_y=None, stagger_s=0.0, force_spots_path=None, max_number_of_lines=100, raw_spots_path=None, wdeg=3, fit_continuum=True, double_precision=False, trace_wdeg=None, trace_wdeg_x=None, trace_wdeg_y=None, trace_per_fiber_deg=None, cpu_threads_per_worker=None, line_search='grid', trace_prior_deg=None, trace_prior_weight=None, trace_prior_ndead_threshold=None, debug_spots=False, masked_amp_ndead_threshold=8000):
    """
    Isolated task for fitting a single bundle.
    """
    t_entry = time.time()
    # trace_wdeg defaults to wdeg (old behavior) -- see fitter.py's
    # PSF_Fitter.fit() docstring/comment and porting-notes.md's
    # r2@20250109 investigation for why this is a separate knob rather
    # than always reusing wdeg (the PSF-shape correction's degree).
    # trace_wdeg_x/trace_wdeg_y further override it per axis -- added
    # after finding X didn't need (and was mildly destabilized by) the
    # same extra curvature that fixed Y.
    trace_wdeg = wdeg if trace_wdeg is None else trace_wdeg
    trace_wdeg_x = trace_wdeg if trace_wdeg_x is None else trace_wdeg_x
    trace_wdeg_y = trace_wdeg if trace_wdeg_y is None else trace_wdeg_y
    if stagger_s > 0:
        time.sleep(stagger_s)
        
    # STRICT ISOLATION: Set before ANY JAX imports in this process.
    # If the parent already restricted CUDA_VISIBLE_DEVICES, map gpu_id
    # within that restriction instead of clobbering it, so multiple driver
    # instances can be pinned to disjoint GPUs from the outside.
    if backend == "gpu":
        existing = os.environ.get("CUDA_VISIBLE_DEVICES")
        if existing:
            devs = existing.split(',')
            os.environ["CUDA_VISIBLE_DEVICES"] = devs[gpu_id % len(devs)]
        else:
            os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
        os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
    else:
        # Hard-exclude CUDA from non-GPU workers. Without this, JAX's CUDA
        # plugin still probes/initializes on whatever GPUs are inherited as
        # visible even with JAX_PLATFORM_NAME=cpu, so concurrent CPU-backend
        # workers collide with each other (and with real GPU workers) over
        # GPU memory and crash with CUDA_ERROR_OUT_OF_MEMORY.
        os.environ["CUDA_VISIBLE_DEVICES"] = ""
        # With no devices visible, the CUDA plugin's own version-check probe
        # (cuInit via cuda_device_count()) fails with CUDA_ERROR_NO_DEVICE --
        # harmless (JAX falls back to CPU regardless) but logged as an ERROR
        # with a full traceback. This is JAX's own supported switch to skip
        # that probe outright instead of just hiding its output.
        os.environ["JAX_SKIP_CUDA_CONSTRAINTS_CHECK"] = "1"
        # JAX_PLATFORM_NAME (set below) is deprecated and, on this JAX
        # version, no longer consulted by backend selection at all -- only
        # JAX_PLATFORMS (plural) is. Without restricting to it, skipping the
        # constraints check above lets the CUDA plugin register successfully,
        # and JAX then genuinely tries to initialize 'cuda' as a real backend
        # (it's highest-priority), which fails hard with a real GPU present
        # but hidden by CUDA_VISIBLE_DEVICES="" -- turning the harmless log
        # noise into a fatal crash instead. Set via jax.config below, not
        # os.environ here -- see the jax_compilation_cache_dir note below for
        # why an os.environ write in this function body is already too late.
        #
        # Thread-count limiting -- must also happen before any JAX import.
        # Without this, JAX's CPU/XLA backend (and whatever BLAS it
        # delegates to) defaults to claiming *all* available hardware
        # threads per process. With N concurrent CPU workers each
        # independently trying to grab every thread, the result is severe
        # intra-node oversubscription/thrashing among the workers
        # themselves -- confirmed directly as the root cause of a real
        # Perlmutter hybrid-allocation pilot going ~3x below its own
        # achievable per-worker throughput while also slowing a concurrent
        # GPU job by 1.56x (see porting-notes.md, 2026-07-21). Scale each
        # worker's thread budget to roughly (available cores / worker
        # count) so N workers collectively stay within the node's real
        # core count instead of each claiming all of them.
        if cpu_threads_per_worker is not None:
            n_threads_str = str(max(1, int(cpu_threads_per_worker)))
            for _env_key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
                os.environ[_env_key] = n_threads_str
            os.environ["XLA_FLAGS"] = (
                os.environ.get("XLA_FLAGS", "")
                + f" --xla_cpu_multi_thread_eigen=true intra_op_parallelism_threads={n_threads_str}"
            ).strip()
    os.environ["JAX_PLATFORM_NAME"] = backend
    # Persistent JAX/XLA compilation cache (analogous to CuPy's .cubin disk
    # cache) -- this driver spawns a fresh process per bundle, so without
    # this every worker pays a full JIT-compile cost from scratch even when
    # an identical (function, array-shape) pair was already compiled by an
    # earlier bundle/camera in the same campaign. Validated in porting-notes.md
    # "JAX persistent compilation cache" -- ~43% faster end-to-end wall time
    # on a warm cache in isolated single-bundle testing, with the previously
    # documented "final joint fit is slower in Python than C++" finding
    # reversing once warm (Python becomes ~3.2x faster on that phase).
    # Respects an operator-set JAX_COMPILATION_CACHE_DIR if present.
    default_cache_dir = os.path.join(os.path.expanduser("~"), ".cache", "specex", "jax_compilation_cache")
    cache_dir = os.environ.get("JAX_COMPILATION_CACHE_DIR", default_cache_dir)
    os.makedirs(cache_dir, exist_ok=True)
    # Mixed precision (float32 Jacobian in the joint-fit accumulate step,
    # float64 everywhere else) is the default -- see porting-notes.md
    # "Mixed precision, tested exactly as directed" for validation (single-
    # bundle chi2 relative error 2.4e-6, full-CCD wavelength RMS matches the
    # float64 pipeline to 4 decimal places, 71% GPU memory cut). Pass
    # --double-precision on the CLI to force full float64 if ever needed.
    os.environ["SPECEX_MIXED_PRECISION"] = "0" if double_precision else "1"
    
    try:
        import jax
        import jax.numpy as jnp
        # NOTE: env-var-based cache config (JAX_COMPILATION_CACHE_DIR etc.)
        # is too late here -- .psf imports jax.numpy at module level, which
        # this multiprocessing 'spawn' worker triggers while resolving this
        # very function *before* its body runs, so JAX's own config is
        # already locked in by the time any os.environ write below would
        # take effect. Use the jax.config API directly instead, which is
        # read fresh at the point of the call.
        jax.config.update("jax_compilation_cache_dir", cache_dir)
        jax.config.update("jax_persistent_cache_min_compile_time_secs", 0)
        jax.config.update("jax_persistent_cache_min_entry_size_bytes", 0)
        if backend != "gpu":
            jax.config.update("jax_platforms", "cpu")
        t_jax_import = time.time()
        print(f"PHASE_TIMING bundle={bid} jax_import={t_jax_import - t_entry:.2f}s", flush=True)

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

        # Add output path for spot/debug-checkpoint writing. Suffixed with
        # the bundle id (matching C++'s own `_NN` per-bundle checkpoint
        # naming, e.g. cppspots_pass4) -- NOT used for the actual merged
        # FITS output (that's written once, separately, by fit_ccd_native
        # via write_python_psf(out_psf_file, ...) after all bundle workers
        # return). Without this suffix, every bundle worker in a multi-
        # bundle/full-CCD run derives the exact same checkpoint filenames
        # from the shared out_psf_file and clobbers every other bundle's
        # spot list -- confirmed 2026-08-04: a 20-bundle full-CCD run's
        # final .pyspots.txt contained only one bundle's ~25 fibers, not
        # all 500, silently (no error, just whichever bundle's worker
        # finished writing last "wins").
        psf.output_psf_path = out_psf_file.replace('.fits', f'_bundle{bid:02d}.fits')
        psf.debug_spots = debug_spots

        lamp_lines = read_lamp_lines(lamp_lines_file)
        t_io = time.time()
        print(f"PHASE_TIMING bundle={bid} image_psf_io={t_io - t_jax_import:.2f}s", flush=True)

        f_min, f_max = bid * 25, (bid + 1) * 25 - 1

        # Detect fibers with essentially no real data anywhere along their
        # trace (a masked/dead CCD amp, not the milder single-bad-column
        # case the trace-prior ndead gate already handles) -- see
        # find_masked_amp_fibers' docstring and porting-notes.md's
        # 2026-08-14 writeup. These are excluded from candidate generation
        # below exactly like an explicitly-broken fiber, and get the same
        # "propagate the input starting-guess PSF, flag STATUS=-1" write-
        # time treatment (io.py's write_python_psf).
        from .fitter import find_masked_amp_fibers
        masked_amp_fibers, _ = find_masked_amp_fibers(psf, f_min, f_max, weight, ndead_threshold=masked_amp_ndead_threshold)
        if masked_amp_fibers:
            print(f"  MASKED-AMP DETECTION: {len(masked_amp_fibers)} fiber(s) in bundle {bid} flagged as "
                  f"no-data (ndead>{masked_amp_ndead_threshold}, contiguous run): {sorted(masked_amp_fibers)}", flush=True)
        explicit_broken_set = set()
        if broken_fibers:
            explicit_broken_set = {int(f) for f in str(broken_fibers).split(",") if f.strip()}
        candidate_exclude_fibers = explicit_broken_set | masked_amp_fibers

        # If EVERY fiber in this bundle is excluded (a bundle fully inside
        # a masked amp, or an all-broken bundle), there's nothing left to
        # fit -- select_bundle_spots_iterative would return an empty spot
        # list and hit the "no spots" error path below, wrongly reporting
        # a genuine no-op (whole bundle correctly propagated from input) as
        # a bundle failure. Short-circuit cleanly instead: write_python_psf
        # already does the right thing for a bundle_results entry with no
        # pc/tc/chi2 at all, as long as 'skip_bundle' tells it to skip the
        # normal per-bundle correction-write block entirely and fall
        # straight through to the explicitly_broken_fibers/masked_amp_fibers
        # pass-through + STATUS=-1 restoration (which doesn't need pc/tc).
        if set(range(f_min, f_max + 1)) <= candidate_exclude_fibers:
            print(f"  Bundle {bid}: all {f_max - f_min + 1} fibers excluded (broken/masked-amp) -- "
                  f"no fit performed, propagating input PSF for the whole bundle.", flush=True)
            return bid, {
                'skip_bundle': True,
                'masked_amp_fibers': sorted(masked_amp_fibers),
                'explicitly_broken_fibers': sorted(explicit_broken_set),
            }

        if force_spots_path:
            # Load spots from file: fiber,wave,xc,yc
            spots = []
            with open(force_spots_path, 'r') as f:
                for line in f:
                    parts = line.strip().split(',')
                    if len(parts) < 4: continue
                    s = {'fiber': int(parts[0]), 'wave': float(parts[1]), 'xc_init': float(parts[2]), 'yc_init': float(parts[3])}
                    # Pre-calculate stamp boundaries to avoid KeyError in fitter.get_bundle_footprint
                    s['stamp_imin'] = int(np.floor(s['xc_init'] + 0.5)) - psf.h_size_x
                    s['stamp_imax'] = int(np.floor(s['xc_init'] + 0.5)) + psf.h_size_x + 1
                    s['stamp_jmin'] = int(np.floor(s['yc_init'] + 0.5)) - psf.h_size_y
                    s['stamp_jmax'] = int(np.floor(s['yc_init'] + 0.5)) + psf.h_size_y + 1
                    spots.append(s)
            print(f"  Forcing {len(spots)} spots from {force_spots_path}", flush=True)
            
            # C++ spots files don't have flux; we must estimate it to initialize the fitter
            # We do this by running the internal spot stats tool once.
            import jax.numpy as jnp
            from .fitter import _get_spot_stats_jax
            c_xc = jnp.array([s['xc_init'] for s in spots])
            c_yc = jnp.array([s['yc_init'] for s in spots])
            gh_all = jnp.array([psf.gh_params(s['fiber'], s['wave']) for s in spots])
            housekeeping_hsize_x = min(3, psf.h_size_x); housekeeping_hsize_y = min(3, psf.h_size_y)
            fluxes, snrs, chi2s, efluxes = _get_spot_stats_jax(jnp.array(image), jnp.array(weight), c_xc, c_yc, gh_all, psf.gh_psf.degree, housekeeping_hsize_x, housekeeping_hsize_y)
            fluxes = np.array(fluxes)
            for i in range(len(spots)):
                spots[i]['flux'] = float(fluxes[i])
        else:
            # Mirrors C++ FitEverything's housekeeping/selection phase: multi-pass
            # reselection from the full candidate list with a trace warm-up loop.
            spots = select_bundle_spots_iterative(psf, f_min, f_max, lamp_lines,
                                                   image, weight, bid,
                                                   broken_fibers=candidate_exclude_fibers,
                                                   max_number_of_lines=max_number_of_lines,
                                                   wdeg=wdeg, fit_continuum=fit_continuum)

        t_select = time.time()
        print(f"PHASE_TIMING bundle={bid} selection={t_select - t_io:.2f}s", flush=True)

        if not spots:
            return bid, {"error": "No spots found for bundle"}

        # A fiber with (near-)zero surviving spots in the final selection is
        # genuinely unconstrained -- C++ detects exactly this condition
        # (specex_psf_fitter.cc:2594-2596, "No selected spot for fiber",
        # trace.mask=3) and excludes that fiber from the trace fit
        # entirely, which (specex_psf_proc.cc:49,58 -- coeff2d starts
        # zero-initialized and the fiber's now-empty coeff array never
        # gets copied in) leaves its XTRACE/YTRACE output row as literal
        # zero. Matched here rather than left to silently interpolate a
        # plausible-looking but never-actually-validated position from
        # neighboring fibers (see porting-notes.md's b3@20241208 bundle-2
        # fiber-65 writeup) -- downstream consumers presumably rely on
        # this all-zero convention to recognize an untrustworthy fiber.
        # Threshold is <2, not strictly 0: C++'s own selection runs an
        # earlier, stricter pass specifically for trace-fitting (not
        # replicated here, which only has one, broader, final pass), so a
        # fiber can have literally 0 spots by C++'s count while Python's
        # broader pass finds exactly 1 -- seen directly on z3@20260401
        # bundle 14/fiber 368 (a single spot at the very top wavelength
        # edge). A single spot can't independently validate any
        # wavelength-dependent trace behavior regardless of which pass
        # found it -- it's the same "too thin to trust" case C++ zeros,
        # just not always caught by an exact spot-count match to C++'s own
        # (unreplicated) selection stages.
        # Explicitly-listed --broken-fibers are a *different* case from a
        # dynamically-discovered zero-spot fiber, and C++ treats them
        # differently too: a broken fiber is simply excluded from the fit
        # entirely and its XTRACE/YTRACE row is left completely untouched
        # at the *input* template's value (confirmed bit-for-bit identical
        # to the input PSF, both for z8@20260401 fibers 473/474 and
        # z3@20260401 fiber 368 -- see porting-notes.md) -- it is NOT
        # zeroed the way a dynamically-discovered zero-spot fiber is
        # (that's the mask=3/resize(0) case above). Excluding a fiber from
        # candidate generation naturally drops its spot count to 0, which
        # would otherwise wrongly pull it into the zeroing set below; skip
        # it here so it's left alone and inherits the input value exactly
        # like C++ does.
        explicitly_broken = explicit_broken_set
        spot_fiber_counts = {}
        for s in spots:
            spot_fiber_counts[s['fiber']] = spot_fiber_counts.get(s['fiber'], 0) + 1
        # masked_amp_fibers excluded here too -- like explicitly_broken,
        # they were kept out of candidate generation entirely (see
        # candidate_exclude_fibers above), so they'd otherwise be wrongly
        # swept into the zero-spot/literal-zero-trace convention below
        # instead of their own pass-through/STATUS=-1 treatment
        # (write_python_psf).
        zero_spot_fibers = [fib for fib in range(f_min, f_max + 1)
                             if spot_fiber_counts.get(fib, 0) < 2
                             and fib not in explicitly_broken
                             and fib not in masked_amp_fibers]

        # trace_prior_weight/trace_prior_ndead_threshold are read by fit()
        # via os.environ (see SPECEX_TRACE_PRIOR_WEIGHT/NDEAD_THRESHOLD in
        # fitter.py) rather than as direct parameters -- set them here, in
        # this already-spawned worker process, right before the call, so a
        # CLI-supplied value takes precedence over the env-only default.
        if trace_prior_weight is not None:
            os.environ["SPECEX_TRACE_PRIOR_WEIGHT"] = str(trace_prior_weight)
        if trace_prior_ndead_threshold is not None:
            os.environ["SPECEX_TRACE_PRIOR_NDEAD_THRESHOLD"] = str(trace_prior_ndead_threshold)

        fitter = PSF_Fitter(psf)
        chi2, pc, tc, cc, final_flux, xc_final, yc_final, spots = fitter.fit(image, weight, spots, bid, max_iter=50, wdeg=wdeg, fit_continuum=fit_continuum, trace_wdeg=trace_wdeg, trace_wdeg_x=trace_wdeg_x, trace_wdeg_y=trace_wdeg_y, trace_per_fiber_deg=trace_per_fiber_deg, line_search=line_search, trace_prior_deg=trace_prior_deg)
        # Reassigning `spots` here (not just capturing it under a new name)
        # is deliberate: everything below this line -- x_orig/y_orig,
        # trace_monomials_abs, the `spots[i]['xc_init'] = ...` refresh loop,
        # pyspots.txt writing, final_selected -- assumes 1:1 correspondence
        # with xc_final/yc_final by index/length. fit() may return a SHORTER
        # list than it was given (SPECEX_MATCH_CPP_DEAD_COLUMN can drop
        # spots), so every downstream use must see that same list, not the
        # original pre-fit one -- see fitter.fit()'s return-statement
        # comment for the real crash this fixes.
        t_finalfit = time.time()
        print(f"PHASE_TIMING bundle={bid} final_joint_fit={t_finalfit - t_select:.2f}s", flush=True)

        # --- Recompute trace_coeffs as the true correction relative to the
        # *original input* trace, fit directly to the final absolute
        # positions (xc_final/yc_final). fitter.fit()'s own `tc` only
        # captures the residual relative to xc_init (the anchor held fixed
        # during the joint fit) -- if xc_init already deviates from the
        # input trace (e.g. after selection's trace-warmup snapping, or
        # when using --force-spots with externally supplied positions),
        # that deviation would otherwise be silently dropped when
        # write_python_psf() adds `tc` onto the input trace to build the
        # output XTRACE/YTRACE.
        from .fitter import get_bundle_monomials_jnp, get_sparse_nz, get_bundle_block_diagonal_trace_monomials
        x_orig = np.array([psf.x_ccd(s['fiber'], s['wave']) for s in spots])
        y_orig = np.array([psf.y_ccd(s['fiber'], s['wave']) for s in spots])
        res_x = np.array(xc_final) - x_orig
        res_y = np.array(yc_final) - y_orig
        if trace_per_fiber_deg is not None:
            # Block-diagonal by fiber (stage 1 of the full per-fiber
            # redesign) -- both axes already share the exact same
            # full-width basis (no freeze-masking, unlike the shared-basis
            # path below), so this is a plain lstsq per axis, no
            # prefix/padding bookkeeping needed.
            trace_monomials_abs = np.array(get_bundle_block_diagonal_trace_monomials(psf, bid, spots, trace_per_fiber_deg))
            tc_x_abs, _, _, _ = np.linalg.lstsq(trace_monomials_abs, res_x, rcond=None)
            tc_y_abs, _, _, _ = np.linalg.lstsq(trace_monomials_abs, res_y, rcond=None)
        else:
            # trace_wdeg_x/trace_wdeg_y (not wdeg) -- this recompute rebuilds tc
            # from scratch against the absolute final positions, so it must use
            # the same basis PSF_Fitter.fit() actually optimized trace_coeffs
            # in, not the (possibly different) PSF-shape basis. x and y can
            # have different degrees (see fitter.py's freeze-mask comment for
            # why): build one shared design matrix at the higher of the two
            # degrees, lstsq each axis against only its own leading-column
            # prefix (get_sparse_nz(1, d) is a strict prefix of any higher
            # degree's), then zero-pad the shorter one back up to the shared
            # width so tc stays a plain (2, N) array like every other caller
            # (write_python_psf included) expects.
            trace_wdeg_shared_abs = max(trace_wdeg_x, trace_wdeg_y)
            trace_monomials_abs = np.array(get_bundle_monomials_jnp(psf, bid, spots, wdeg=trace_wdeg_shared_abs))
            npoly_x_abs = len(get_sparse_nz(1, trace_wdeg_x)); npoly_y_abs = len(get_sparse_nz(1, trace_wdeg_y))
            tc_x_fit, _, _, _ = np.linalg.lstsq(trace_monomials_abs[:, :npoly_x_abs], res_x, rcond=None)
            tc_y_fit, _, _, _ = np.linalg.lstsq(trace_monomials_abs[:, :npoly_y_abs], res_y, rcond=None)
            npoly_shared_abs = trace_monomials_abs.shape[1]
            tc_x_abs = np.zeros(npoly_shared_abs); tc_x_abs[:npoly_x_abs] = tc_x_fit
            tc_y_abs = np.zeros(npoly_shared_abs); tc_y_abs[:npoly_y_abs] = tc_y_fit
        tc = np.stack([tc_x_abs, tc_y_abs], axis=0)

        # --- C++ Parity: Update spots with refined model centroids ---
        # The C++ implementation snaps final spots to the model PSF position
        # We need to preserve the original raw values for the debug file
        raw_centroids = [(s['xc_init'], s['yc_init']) for s in spots]
        
        # Use the xc_final/yc_final already computed by the fitter
        for i in range(len(spots)):
            spots[i]['xc_init'] = float(xc_final[i])
            spots[i]['yc_init'] = float(yc_final[i])

        # CRITICAL: Overwrite pyspots.txt with refined centroids.
        # get_bundle_spots wrote the raw selection; we must update it with fit results.
        if debug_spots and hasattr(psf, 'output_psf_path') and psf.output_psf_path:
            spots_path = psf.output_psf_path.replace('.fits', '.pyspots.txt')
            with open(spots_path, 'w') as f:
                for s in spots:
                    f.write(f"{s['fiber']},{s['wave']:.15f},{s['xc_init']:.15f},{s['yc_init']:.15f}\n")

        # Debug export: verify that xc_final/yc_final differ from initial raw values
        # (same per-bundle clobbering issue as psf.output_psf_path above -- fixed the same way)
        if debug_spots:
            debug_path = out_psf_file.replace('.fits', f'_bundle{bid:02d}.refined_centroids_debug.txt')
            with open(debug_path, 'w') as f:
                for i in range(len(spots)):
                    f.write(f"spot {i}: raw({raw_centroids[i][0]:.15f}, {raw_centroids[i][1]:.15f}) -> refined({float(xc_final[i]):.15f}, {float(yc_final[i]):.15f})\n")

        # We no longer re-run selection with the refined centroids as a second pass.
        # This matches the iterative snapping logic now implemented inside the fitter,
        # and prevents the "Selected is subset of Raw: True" behavior.
        # We simply update the 'spots' metadata for the final result.
        final_selected = spots
        for i in range(len(final_selected)):
            final_selected[i]['flux'] = float(final_flux[i])
            final_selected[i]['xc_init'] = float(xc_final[i])
            final_selected[i]['yc_init'] = float(yc_final[i])
        
        t_postproc = time.time()
        print(f"PHASE_TIMING bundle={bid} postproc={t_postproc - t_finalfit:.2f}s total={t_postproc - t_entry:.2f}s", flush=True)

        # Keep the cross-process payload minimal to avoid multiprocessing
        # pipe/pickling overhead.
        return bid, {

            'psf_coeffs': np.array(pc),
            'trace_coeffs': np.array(tc),
            'continuum_coeffs': np.array(cc),
            'wdeg': wdeg,
            # The actual width tc was built at (max(trace_wdeg_x,
            # trace_wdeg_y), not the shared-default trace_wdeg convenience
            # value above) -- write_python_psf indexes tc with this, and
            # the shorter axis's unused trailing columns are exact zeros,
            # so using the wider of the two here is required for correct
            # reconstruction, not just for the wider axis's own sake.
            # Meaningless (ignored by write_python_psf) when
            # trace_per_fiber_deg is set instead.
            'trace_wdeg': None if trace_per_fiber_deg is not None else trace_wdeg_shared_abs,
            # Stage 1 of the full per-fiber redesign (see porting-notes.md):
            # tells write_python_psf to use the per-fiber (not
            # broadcast-across-fibers) write-back path.
            'trace_per_fiber_deg': trace_per_fiber_deg,
            'zero_spot_fibers': zero_spot_fibers,
            'explicitly_broken_fibers': sorted(explicitly_broken),
            'masked_amp_fibers': sorted(masked_amp_fibers),
            'chi2': float(chi2),
            's_fiber': np.array([s['fiber'] for s in final_selected]),
            's_wave': np.array([s['wave'] for s in final_selected]),
            's_flux': np.array([s['flux'] for s in final_selected])
        }
    except Exception as e:
        import traceback
        err_msg = traceback.format_exc()
        # JAX/XLA's OOM exception (jaxlib.xla_extension.XlaRuntimeError)
        # always includes this token in its message -- distinguishing it
        # from a genuine bug lets fit_ccd_native retry just these bundles
        # at a lower packing instead of dropping them silently.
        is_oom = "RESOURCE_EXHAUSTED" in err_msg
        print(f"FAILED Bundle {bid} on {backend.upper()} {gpu_id}:\n{err_msg}", flush=True)
        return bid, {"error": str(e), "traceback": err_msg, "is_oom": is_oom}


def fit_ccd_native(arc_file, in_psf_file, out_psf_file, lamp_lines_file,
                     first_bundle=0, last_bundle=19, n_gpus=4, backend="gpu",
                     broken_fibers=None, sn_threshold=3.0, h_size_y=5, force_spots_path=None, max_number_of_lines=100,
                     workers_per_gpu=None, cpu_workers=None, gpu_worker_threads=None, legendre_deg_wave=None, fit_continuum=None, double_precision=False,
                     trace_legendre_deg_wave=None, trace_legendre_deg_wave_x=None, trace_legendre_deg_wave_y=None,
                     trace_per_fiber_deg=6, trace_prior_deg=1, trace_prior_weight=None, trace_prior_ndead_threshold=None,
                     line_search='grid', debug_spots=False, masked_amp_ndead_threshold=8000):
    """
    Fits a full CCD (20 bundles) using parallel processes.

    For backend="gpu", multiple worker processes are packed onto each
    physical GPU (workers_per_gpu) since a single bundle fit doesn't
    saturate an A100 -- benchmarked on bundle 5 (z8/00344649): 1/GPU
    takes ~54s/bundle, 4/GPU takes ~75s/bundle (~1.4x slower) but yields
    ~2.7x more aggregate throughput. 5/GPU reliably hits
    RESOURCE_EXHAUSTED (peak ~8.6GB/worker x 5 exceeds the 40GB A100),
    so 4/GPU is the validated safe ceiling. The pool size need not equal
    the bundle count -- Pool.starmap dynamically queues remaining tasks
    onto whichever worker frees up first.

    legendre_deg_wave/fit_continuum default to None, which auto-selects
    real C++ production defaults (desispec/scripts/specex.py:224-228) by
    detecting the band from the input image's CAMERA header: degree 3 +
    continuum for z-band, degree 1 + no continuum otherwise. Pass either
    explicitly to override (matching desi_psf_fit's own CLI override
    behavior) for controlled A/B comparisons.

    trace_legendre_deg_wave_x/trace_legendre_deg_wave_y independently set
    the wavelength degree of the trace-position correction per axis,
    leaving legendre_deg_wave's value governing just the PSF-shape
    (Gauss-Hermite) correction. trace_legendre_deg_wave (no suffix), if
    given, sets *both* axes at once as a convenience override; otherwise
    each axis defaults independently: X to legendre_deg_wave's own value
    (1 for b/r, 3 for z -- i.e. unchanged from the pre-decoupling
    behavior), Y to 2 for b/r bands and legendre_deg_wave's value (3) for
    z-band. See porting-notes.md's r2@20250109 investigation: raising the
    *shared* wdeg to give the trace fit more wavelength curvature also
    handed the PSF-shape fit the same extra freedom, opening a
    trace-position/PSF-asymmetry degeneracy that made xrms worse even as
    it fixed yrms; decoupling trace from PSF-shape and validating
    trace_wdeg=2 against the real C++ engine (run_specex()) on 6 bundles
    across both flagged exposures plus 3 on a clean control exposure found
    a clean win on Y (yrms cut ~68%) but a smaller, real xrms cost on 2 of
    those 6 bundles. Root-caused to X sharing the same raised degree as Y
    even though X's own true residual (checked against C++) is well
    described by degree 1 already -- decoupling X and Y within the trace
    correction (this parameter split) isolates the extra freedom to Y
    only. z-band was not part of that validation, so its trace correction
    stays coupled to its own wdeg (3) on both axes unless overridden.

    trace_per_fiber_deg (default 6, as of 2026-08-05) replaces the shared
    low-degree trace_legendre_deg_wave_x/y basis entirely with a
    block-diagonal-by-fiber one -- each of the bundle's 25 fibers gets
    its own independent (trace_per_fiber_deg+1)-term wavelength basis
    with zero cross-fiber sharing, matching C++'s per-fiber trace
    parameter count. Pass 0/None to fall back to the old shared basis.
    Paired with trace_prior_deg (default 1), an ndead-gated cross-fiber
    regularization on the per-fiber coefficients at that Legendre degree
    and above (ported from C++'s own trace_prior_deg mechanism,
    specex_psf_fitter.cc -- off by default in C++ itself, but a real,
    validated fix here for the handful of fibers with severe local dead-
    pixel contamination that per-fiber independence alone handles badly;
    see porting-notes.md). Pass a negative trace_prior_deg to disable just
    the prior while keeping per-fiber trace on. **Validated on a
    definitive 30-CCD isolated-JAX-cache campaign, 2026-08-05**: 30/30
    cases improved on both xrms (mean -37.5%) and yrms (mean -60.7%)
    vs. the old shared-basis default, for a ~9%% aggregate timing cost
    (b +6.5%, r -0.3%, z +20.7% -- still 3-4x+ faster than C++ overall).
    Two known, small, already-understood residual limitations remain
    (not blocking): a handful of fibers with extreme dead-pixel counts
    (ndead>>threshold) only partially respond even with the prior active
    (a genuine data floor, not a bug), and bundle-boundary fibers improve
    less than interior fibers under full per-fiber independence (no
    cross-fiber sharing to lean on at the edge) -- see porting-notes.md's
    2026-08-05 entries for both.

    line_search (default 'grid') selects the final joint fit's per-
    iteration step-size search: 'grid' is the long-standing coarse
    3-point [0.2,0.5,1.0] search; 'brent' is a continuous but NOT
    C++-faithful search (wrong bracket/tolerance, kept for reference,
    tested negative); 'cpp' is a faithful replica of C++'s actual
    algorithm (specex_psf_fitter.cc/specex_brent.cc -- mode-dependent
    skip logic plus a direct Numerical Recipes brent() port). Both
    'brent' and 'cpp' were tested on 2 hard + 2 normal bundles and found
    to produce no meaningful xrms/yrms change vs 'grid' -- see
    porting-notes.md. Experimental/opt-in, not the default.
    """
    t_start = time.time()
    all_bundles = range(first_bundle, last_bundle + 1)
    bundle_results = {}

    # trace_per_fiber_deg<=0 means "off" (fall back to the shared trace
    # basis); trace_prior_deg<0 means "on but with the cross-fiber prior
    # disabled" -- both CLI-level escape hatches from the 2026-08-05
    # promoted defaults (6 / 1), see --trace-per-fiber-deg/--trace-prior-deg
    # help text.
    if trace_per_fiber_deg is not None and trace_per_fiber_deg <= 0:
        trace_per_fiber_deg = None
    if trace_prior_deg is not None and trace_prior_deg < 0:
        trace_prior_deg = None

    if legendre_deg_wave is None or fit_continuum is None or trace_legendre_deg_wave_x is None or trace_legendre_deg_wave_y is None or workers_per_gpu is None:
        import fitsio
        cam = fitsio.read_header(arc_file, ext=0)['CAMERA'].strip().lower()
        band = cam[0]
        if workers_per_gpu is None:
            # z-band's larger per-fiber design matrix (~350 params/bundle)
            # combined with its higher spot density hits GPU
            # RESOURCE_EXHAUSTED at 5 workers/GPU specifically -- see
            # porting-notes.md's 2026-08-05 OOM investigation. 3/GPU fully
            # avoids it (validated OOM-free across all 10 z-band cases in
            # the definitive 30-CCD campaign) at a modest packing cost.
            workers_per_gpu = 3 if (band == 'z' and trace_per_fiber_deg is not None) else 5
        if legendre_deg_wave is None:
            legendre_deg_wave = 3 if band == 'z' else 1
        if fit_continuum is None:
            fit_continuum = (band == 'z')
        if trace_legendre_deg_wave_x is None:
            trace_legendre_deg_wave_x = trace_legendre_deg_wave if trace_legendre_deg_wave is not None else legendre_deg_wave
        if trace_legendre_deg_wave_y is None:
            # b/r default: 2, not 1 -- validated against the real C++ engine
            # (run_specex(), now working locally -- see porting-notes.md)
            # across 6 bundles on both flagged exposures (r2@20250109,
            # r2@20241208): mean yrms 0.2182px -> 0.0686px (68% cut, into
            # normal-case territory). z-band kept coupled to its own wdeg
            # (3) -- not part of this validation.
            trace_legendre_deg_wave_y = trace_legendre_deg_wave if trace_legendre_deg_wave is not None else (2 if band != 'z' else legendre_deg_wave)

    print(f"--- SPECE-X Multi-Process CCD Fit ({backend.upper()}) ---")
    print(f"  Arc: {arc_file}")
    print(f"  In PSF: {in_psf_file}")
    print(f"  Out PSF: {out_psf_file}")
    print(f"  legendre-deg-wave: {legendre_deg_wave}  trace-legendre-deg-wave: x={trace_legendre_deg_wave_x} y={trace_legendre_deg_wave_y}  fit-continuum: {fit_continuum}")
    if trace_per_fiber_deg is not None:
        print(f"  trace-per-fiber-deg: {trace_per_fiber_deg}  trace-prior-deg: {trace_prior_deg}  workers-per-gpu: {workers_per_gpu}")
    if broken_fibers:
        print(f"  Broken Fibers: {broken_fibers}")

    available_cores = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)

    if backend == "gpu":
        packing = workers_per_gpu  # bundle-fit workers packed onto each GPU
        # Normally unset -- GPU-backend workers' host-side NumPy/BLAS calls
        # (the selection/housekeeping phase, ~70% of a bundle's wall time,
        # see porting-notes.md) default to unconstrained thread counts,
        # unlike --backend cpu workers which always get a computed budget
        # below. --gpu-worker-threads lets this be forced explicitly, to
        # test whether that default threading is (a) load-bearing for
        # per-worker speed or (b) pure oversubscription noise that's
        # crowding out any concurrent CPU-backend work -- see the 2026-08-09
        # profiling session in porting-notes.md.
        gpu_thread_cap = gpu_worker_threads
        if gpu_thread_cap is not None:
            print(f"  GPU-worker thread cap: {gpu_thread_cap} threads/worker (forced via --gpu-worker-threads)", flush=True)
    else:
        packing = cpu_workers or n_gpus  # concurrent CPU-backend worker processes

    def _cpu_threads_for(n_workers_this):
        # See fit_bundle_task's thread-limiting comment for why this exists
        # at all: without it, every CPU worker independently claims all
        # available cores, and N concurrent workers thrash each other.
        # Scaling to (real core count / worker count) keeps the pool's
        # aggregate thread demand within the node's actual budget. Floor
        # of 1 thread/worker (an oversubscribed-but-not-zero fallback) if
        # there happen to be more workers than cores. Recomputed per batch
        # (not just once up front) so an OOM-retry batch, which runs with
        # fewer workers, correctly gets a bigger per-worker thread budget.
        if backend == "gpu":
            return gpu_thread_cap
        t = max(1, available_cores // n_workers_this)
        print(f"  CPU thread budget: {available_cores} cores / {n_workers_this} workers = {t} threads/worker", flush=True)
        return t

    def _run_batch(bundle_ids, n_workers_this, cpu_threads_this):
        ctx = mp.get_context('spawn')
        with ctx.Pool(processes=n_workers_this) as pool:
            tasks = []
            for i, bid in enumerate(bundle_ids):
                gpu_id = i % n_gpus
                # Use 2s stagger to prevent JIT compilation contention on CPU
                stagger_s = i * 2.0 if backend == "cpu" else 0.0
                tasks.append((bid, gpu_id, arc_file, in_psf_file, out_psf_file, lamp_lines_file, backend, broken_fibers, sn_threshold, h_size_y, stagger_s, force_spots_path, max_number_of_lines, None, legendre_deg_wave, fit_continuum, double_precision, None, trace_legendre_deg_wave_x, trace_legendre_deg_wave_y, trace_per_fiber_deg, cpu_threads_this, line_search, trace_prior_deg, trace_prior_weight, trace_prior_ndead_threshold, debug_spots, masked_amp_ndead_threshold))
            return pool.starmap(fit_bundle_task, tasks)

    # Bundles are independent tasks (Pool.starmap queues them dynamically),
    # so a bundle that fails with a GPU OOM can simply be resubmitted in a
    # smaller follow-up batch at reduced packing -- no need to restart the
    # whole CCD. Non-OOM failures are NOT retried (retrying a real bug just
    # wastes GPU time and reproduces the same failure); only OOM gets this
    # treatment. Capped at 2 retry rounds, halving packing each time (floor
    # of 1), so a bundle that's simply too large to ever fit still fails
    # fast rather than looping.
    pending = list(all_bundles)
    failed_bundles = {}
    attempt = 0
    max_oom_retries = 2
    while pending:
        n_workers_this = min(len(pending), n_gpus * packing) if backend == "gpu" else min(len(pending), packing)
        cpu_threads_this = _cpu_threads_for(n_workers_this)
        label = "initial" if attempt == 0 else f"OOM-retry {attempt}"
        print(f"Launching {len(pending)} bundles across {n_workers_this} workers ({label}, packing={packing})...", flush=True)
        chunk_results = _run_batch(pending, n_workers_this, cpu_threads_this)

        oom_bids = []
        pending = []
        for bid, res in chunk_results:
            if "error" in res:
                if res.get("is_oom"):
                    oom_bids.append(bid)
                else:
                    print(f"WARNING: Bundle {bid} failed: {res['error']}")
                    failed_bundles[bid] = res["error"]
            else:
                bundle_results[bid] = res

        if not oom_bids:
            break
        if attempt >= max_oom_retries or packing <= 1:
            for bid in oom_bids:
                print(f"WARNING: Bundle {bid} failed: GPU OOM persisted after {attempt} retry round(s) down to packing={packing}")
                failed_bundles[bid] = "GPU OOM persisted after retries"
            break

        packing = max(1, packing // 2)
        print(f"OOM RETRY: {len(oom_bids)} bundle(s) {sorted(oom_bids)} hit GPU OOM (RESOURCE_EXHAUSTED); retrying at packing={packing}", flush=True)
        pending = oom_bids
        attempt += 1

    print(f"Total CCD Fit Time: {time.time() - t_start:.2f}s")

    # Merge: C++'s desi_compute_psf --mpi fits each bundle independently
    # and merge_psf() (desispec/scripts/specex.py) does a straight
    # per-fiber copy of each bundle's XTRACE/YTRACE/PSF coefficients into
    # the shared output arrays -- no cross-bundle smoothing or refit.
    # write_python_psf already does this per-bundle slice-copy correctly;
    # a prior "Phase 2: Global Wavelength Refinement" here fabricated a
    # smoothed CCD-wide re-fit (with hardcoded, non-z-band wavelengths and
    # unexplained scale factors) that overwrote every bundle's real fitted
    # trace_coeffs and wrote a bogus WAVECORR table. That table's real-world
    # counterpart (WAVE/DWAVE/DWAVE_ERR as 'EXTOFF') comes from desispec's
    # trace_shifts.py, a separate downstream pipeline stage that runs after
    # desi_compute_psf on sky/arc reference lines -- desi_compute_psf itself
    # never produces it, so there is nothing to replicate here.
    if bundle_results and out_psf_file:
        write_python_psf(out_psf_file, bundle_results, in_psf_file)

    n_total = len(all_bundles)
    if failed_bundles:
        print(f"SPECEX_RESULT: FAILED {len(failed_bundles)}/{n_total} bundles: {sorted(failed_bundles)}", flush=True)
    else:
        print(f"SPECEX_RESULT: OK {n_total}/{n_total} bundles", flush=True)
    return failed_bundles

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Specex Python/JAX PSF Fitter")
    parser.add_argument("-a", "--arc", "--input-image", type=str, required=True, help="Input preproc arc image")
    parser.add_argument("--in-psf", "--input-psf", type=str, required=True, help="Input (shifted) PSF file")
    parser.add_argument("--out-psf", "--output-psf", type=str, required=True, help="Output PSF file")
    parser.add_argument("--lamp-lines", type=str, help="Lamp lines file")
    parser.add_argument("--first-bundle", type=int, default=0)
    parser.add_argument("--last-bundle", type=int, default=19)
    parser.add_argument("--first-fiber", type=int, help="First fiber to fit")
    parser.add_argument("--last-fiber", type=int, help="Last fiber to fit")
    parser.add_argument("--legendre-deg-wave", type=int, default=None, help="Legendre degree for the joint fit's PSF-shape wavelength basis (default: auto, matching real C++ production -- 3 for z-band, 1 otherwise, detected from the input image's CAMERA header).")
    parser.add_argument("--trace-legendre-deg-wave", type=int, default=None, help="Legendre degree for the joint fit's trace-position wavelength basis, both axes at once (independent of --legendre-deg-wave's PSF-shape degree). Overridden per-axis by --trace-legendre-deg-wave-x/-y if either is also given. Default: auto per axis -- see those flags' help.")
    parser.add_argument("--trace-legendre-deg-wave-x", type=int, default=None, help="Legendre degree for the trace-position X basis only (default: auto -- same as --legendre-deg-wave, i.e. unchanged from the pre-decoupling behavior; X was found not to need the extra curvature Y does -- see porting-notes.md's r2@20250109 investigation)")
    parser.add_argument("--trace-legendre-deg-wave-y", type=int, default=None, help="Legendre degree for the trace-position Y basis only (default: auto -- 2 for b/r bands, same as --legendre-deg-wave for z-band; validated against the real C++ engine -- see porting-notes.md's r2@20250109 investigation)")
    parser.add_argument("--trace-per-fiber-deg", type=int, default=6, help="Replaces the shared trace basis with a block-diagonal-by-fiber one at this wavelength degree (6 matches the input PSF's own native trace degree, and C++'s per-fiber parameter count), paired with an ndead-gated cross-fiber trace-coefficient prior (see --trace-prior-* below). DEFAULT AS OF 2026-08-05: on at degree 6 -- validated on a definitive 30-CCD isolated-cache campaign (porting-notes.md), 30/30 cases improved on both xrms (-37.5%% mean) and yrms (-60.7%% mean) vs the old shared-basis default, for a ~9%% timing cost. Pass 0 to fall back to the old shared trace_legendre_deg_wave_x/y basis.")
    parser.add_argument("--trace-prior-deg", type=int, default=1, help="Legendre degree at/above which --trace-per-fiber-deg's per-fiber coefficients are pulled toward the bundle's cross-fiber consensus (C++'s trace prior, specex_psf_fitter.cc, ported and ndead-gated -- see porting-notes.md). Only active when --trace-per-fiber-deg is on. Default 1 (each fiber's own physical position, degree 0, stays fully independent; only higher-order shape terms are regularized). Pass a negative value to disable the prior entirely while keeping per-fiber trace on.")
    parser.add_argument("--trace-prior-weight", type=float, default=1e5, help="Trace-prior penalty weight (C++'s own hardcoded value, 1e8, was found to measurably harm healthy bundles when applied blanket-style -- see porting-notes.md's weight sweep). Only matters for fibers flagged by --trace-prior-ndead-threshold.")
    parser.add_argument("--trace-prior-ndead-threshold", type=int, default=500, help="A fiber's C++-style dead-pixel count (ndead) above this triggers the trace prior for that fiber only; fibers below it are completely unaffected (bit-identical to no-prior). 500 comfortably separates normal fibers (ndead ~20-120) from the known bad cases (ndead ~2300-17500).")
    parser.add_argument("--masked-amp-ndead-threshold", type=int, default=8000, help="A fiber's ndead above this, PLUS a contiguous run of >=3 such fibers, marks it as overlapping a masked/dead CCD amp (no real data at all, not the milder single-bad-column case --trace-prior-ndead-threshold handles) -- the fiber is excluded from the fit entirely, its input starting-guess PSF is propagated unchanged, and its STATUS is set to -1, matching real C++'s own observed behavior. FIRST-PASS HEURISTIC: calibrated against one real case (r8@20211028/00106399's amp-A mask, see porting-notes.md's 2026-08-14 writeup) -- treat as tunable, not load-bearing precision.")
    parser.add_argument("--fit-continuum", action=argparse.BooleanOptionalAction, default=None, help="Fit a per-bundle continuum background (default: auto, matching real C++ production -- on for z-band, off otherwise)")
    parser.add_argument("--gpu", type=int, default=4, help="Number of GPUs to use")
    parser.add_argument("--workers-per-gpu", type=int, default=None, help="Concurrent bundle-fit worker processes packed onto each GPU. Default: auto -- 5, except 3 for z-band when --trace-per-fiber-deg is active (its larger per-fiber design matrix hits GPU RESOURCE_EXHAUSTED at 5/GPU on z-band specifically -- see porting-notes.md's OOM investigation). Pass explicitly to override.")
    parser.add_argument("--cpu-workers", type=int, help="Concurrent worker processes for --backend cpu (default: --gpu count)")
    parser.add_argument("--gpu-worker-threads", type=int, default=None, help="Force an OMP/BLAS/XLA thread cap on each --backend gpu worker's host-side (CPU) computation, mirroring --backend cpu's own auto-computed budget. Default: unconstrained (each worker's BLAS calls may claim all visible threads). Diagnostic flag for probing whether GPU-worker host threading is load-bearing or pure oversubscription -- see porting-notes.md.")
    parser.add_argument("--backend", type=str, default="gpu", choices=["cpu", "gpu"])
    parser.add_argument("--broken-fibers", type=str, help="Comma-separated list of broken fibers")
    parser.add_argument("--sn-threshold", type=float, default=3.0, help="S/N threshold for spot selection")
    parser.add_argument("--max-lines", type=int, default=200, help="Maximum number of lines to keep per bundle")
    parser.add_argument("--h-size-y", type=int, default=5, help="Override PSF stamp half-size in Y")
    parser.add_argument("--force-spots", type=str, help="Path to a file containing spots to fit (fiber,wave,xc,yc)")
    parser.add_argument("--double-precision", action="store_true", help="Force full float64 precision for the joint-fit Jacobian (default: mixed float32/float64 -- see porting-notes.md; validated equivalent accuracy, ~71%% less GPU memory/worker)")
    parser.add_argument("--line-search", type=str, default="grid", choices=["grid", "brent", "cpp"], help="EXPERIMENTAL: the final joint fit's per-iteration step-size search. 'grid' (default): the long-standing coarse 3-point [0.2,0.5,1.0] search. 'brent': a continuous but NOT C++-faithful search, kept for reference. 'cpp': a faithful replica of C++'s actual algorithm (mode-dependent skip logic + a direct Numerical Recipes brent() port, see specex_psf_fitter.cc/specex_brent.cc). Both 'brent' and 'cpp' tested negative (no xrms/yrms change on 2 hard + 2 normal bundles) -- see porting-notes.md.")
    parser.add_argument("--debug-spots", action="store_true", help="Write per-pass spot-selection debug dump files (.pyrawspots.txt, .pyspots_pass*.txt, .pyrawspots_final.txt, .pyspots.txt, .refined_centroids_debug.txt), the direct Python analog of C++'s --debug-spots. Off by default -- adds I/O overhead (one set of files per bundle worker) with no effect on the fitted output.")

    args = parser.parse_args()
    
    if not args.lamp_lines:
        base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        args.lamp_lines = os.path.join(base, "specex/data/specex_linelist_desi.txt")

    os.environ["JAX_PLATFORM_NAME"] = args.backend  # deprecated/ineffective on this JAX version, see below

    if args.backend == "gpu":
        # Fail fast with an actionable message instead of either (a) JAX's
        # own opaque "Unknown backend: 'gpu' requested... Platforms present
        # are: cpu" RuntimeError, which gives no hint at the actual cause, or
        # (b) silently falling back to CPU and running ~10x slower with no
        # indication anything is wrong. --gpu is the default and this is a
        # perf-critical batch pipeline, so a silent fallback would be worse
        # than a loud failure -- see how-to-run.md Section 0 for the real
        # fix (confirmed root cause once: an invalid pip extra name, e.g.
        # "jax[cuda13_pip]", is not a hard error -- pip only *warns* and
        # silently installs a CPU-only jaxlib).
        import jax
        if not any(d.platform == "gpu" for d in jax.devices()):
            raise RuntimeError(
                "--backend gpu (the default) was requested, but no GPU-capable JAX "
                "platform is available -- jax.devices() found only "
                f"{sorted(set(d.platform for d in jax.devices()))}. This almost always means "
                "jaxlib was installed without CUDA support (see how-to-run.md Section 0). "
                "Reinstall with `pip install --upgrade \"jax[cuda13]\"` and confirm "
                "`python -c \"import jax; print(jax.devices())\"` reports a CudaDevice, "
                "or pass --backend cpu to run on CPU deliberately."
            )

    if args.backend != "gpu":
        # STRICT ISOLATION for the main process itself, mirroring
        # fit_bundle_task's per-worker isolation above. Without this, the
        # main process's own post-pool JAX usage (write_python_psf ->
        # legendre_pol_jnp, io.py) has no backend restriction at all and
        # defaults to JAX's normal CUDA-first device selection regardless
        # of --backend cpu (the JAX_PLATFORM_NAME line above doesn't help --
        # it's the deprecated singular name; only JAX_PLATFORMS, plural, is
        # consulted, and even that needs to be set via jax.config, not
        # os.environ, once jax is already imported -- see the matching
        # comment in fit_bundle_task). Harmless when no other GPU job is
        # running (CUDA init just succeeds or falls back cleanly), but under
        # real concurrent CPU+GPU production use this touches CUDA while
        # concurrent GPU-backend jobs have already exhausted GPU memory --
        # confirmed directly to crash with CUDA_ERROR_OUT_OF_MEMORY on all
        # visible devices (see porting-notes.md), and the prime suspect for
        # an earlier session's CPU+GPU hybrid deadlock (same code path, a
        # hang instead of a crash is plausible under different CUDA-driver
        # contention timing).
        os.environ["CUDA_VISIBLE_DEVICES"] = ""
        os.environ["JAX_SKIP_CUDA_CONSTRAINTS_CHECK"] = "1"
        import jax
        jax.config.update("jax_platforms", "cpu")

    failed_bundles = fit_ccd_native(
        arc_file=args.arc,
        in_psf_file=args.in_psf,
        out_psf_file=args.out_psf,
        lamp_lines_file=args.lamp_lines,
        first_bundle=args.first_bundle,
        last_bundle=args.last_bundle,
        n_gpus=args.gpu,
        backend=args.backend,
        broken_fibers=args.broken_fibers,
        sn_threshold=args.sn_threshold,
        max_number_of_lines=args.max_lines,
        h_size_y=args.h_size_y,
        force_spots_path=args.force_spots,
        workers_per_gpu=args.workers_per_gpu,
        cpu_workers=args.cpu_workers,
        gpu_worker_threads=args.gpu_worker_threads,
        legendre_deg_wave=args.legendre_deg_wave,
        trace_legendre_deg_wave=args.trace_legendre_deg_wave,
        trace_legendre_deg_wave_x=args.trace_legendre_deg_wave_x,
        trace_legendre_deg_wave_y=args.trace_legendre_deg_wave_y,
        trace_per_fiber_deg=args.trace_per_fiber_deg,
        trace_prior_deg=args.trace_prior_deg,
        trace_prior_weight=args.trace_prior_weight,
        trace_prior_ndead_threshold=args.trace_prior_ndead_threshold,
        masked_amp_ndead_threshold=args.masked_amp_ndead_threshold,
        fit_continuum=args.fit_continuum,
        double_precision=args.double_precision,
        line_search=args.line_search,
        debug_spots=args.debug_spots
    )

    if failed_bundles:
        # Previously this always exited 0 even when bundles were silently
        # dropped from the output -- rc==0 alone was never sufficient to
        # confirm a real success (see porting-notes.md's OOM investigation).
        # A non-zero exit here lets callers (run_night.py, desi_proc) tell
        # a genuine failure apart from success without grepping logs.
        sys.exit(1)

if __name__ == "__main__":
    main()
