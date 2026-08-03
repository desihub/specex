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
def fit_bundle_task(bid, gpu_id, arc_file, in_psf_file, out_psf_file, lamp_lines_file, backend="gpu", broken_fibers=None, sn_threshold=3.0, h_size_y=None, stagger_s=0.0, force_spots_path=None, max_number_of_lines=100, raw_spots_path=None, wdeg=3, fit_continuum=True, double_precision=False, trace_wdeg=None, trace_wdeg_x=None, trace_wdeg_y=None, trace_per_fiber_deg=None, cpu_threads_per_worker=None, line_search='grid'):
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

        # Add output path for spot writing
        psf.output_psf_path = out_psf_file

        lamp_lines = read_lamp_lines(lamp_lines_file)
        t_io = time.time()
        print(f"PHASE_TIMING bundle={bid} image_psf_io={t_io - t_jax_import:.2f}s", flush=True)

        f_min, f_max = bid * 25, (bid + 1) * 25 - 1
        
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
                                                   broken_fibers=broken_fibers,
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
        explicitly_broken = set()
        if broken_fibers:
            explicitly_broken = {int(f) for f in str(broken_fibers).split(",") if f.strip()}
        spot_fiber_counts = {}
        for s in spots:
            spot_fiber_counts[s['fiber']] = spot_fiber_counts.get(s['fiber'], 0) + 1
        zero_spot_fibers = [fib for fib in range(f_min, f_max + 1)
                             if spot_fiber_counts.get(fib, 0) < 2 and fib not in explicitly_broken]

        fitter = PSF_Fitter(psf)
        chi2, pc, tc, cc, final_flux, xc_final, yc_final = fitter.fit(image, weight, spots, bid, max_iter=50, wdeg=wdeg, fit_continuum=fit_continuum, trace_wdeg=trace_wdeg, trace_wdeg_x=trace_wdeg_x, trace_wdeg_y=trace_wdeg_y, trace_per_fiber_deg=trace_per_fiber_deg, line_search=line_search)
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
        if hasattr(psf, 'output_psf_path') and psf.output_psf_path:
            spots_path = psf.output_psf_path.replace('.fits', '.pyspots.txt')
            with open(spots_path, 'w') as f:
                for s in spots:
                    f.write(f"{s['fiber']},{s['wave']:.15f},{s['xc_init']:.15f},{s['yc_init']:.15f}\n")

        # Debug export: verify that xc_final/yc_final differ from initial raw values
        debug_path = out_psf_file.replace('.fits', '.refined_centroids_debug.txt')
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
            'chi2': float(chi2),
            's_fiber': np.array([s['fiber'] for s in final_selected]),
            's_wave': np.array([s['wave'] for s in final_selected]),
            's_flux': np.array([s['flux'] for s in final_selected])
        }
    except Exception as e:
        import traceback
        err_msg = traceback.format_exc()
        print(f"FAILED Bundle {bid} on {backend.upper()} {gpu_id}:\n{err_msg}", flush=True)
        return bid, {"error": str(e), "traceback": err_msg}


def fit_ccd_native(arc_file, in_psf_file, out_psf_file, lamp_lines_file,
                     first_bundle=0, last_bundle=19, n_gpus=4, backend="gpu",
                     broken_fibers=None, sn_threshold=3.0, h_size_y=5, force_spots_path=None, max_number_of_lines=100,
                     workers_per_gpu=4, cpu_workers=None, legendre_deg_wave=None, fit_continuum=None, double_precision=False,
                     trace_legendre_deg_wave=None, trace_legendre_deg_wave_x=None, trace_legendre_deg_wave_y=None,
                     trace_per_fiber_deg=None, line_search='grid'):
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

    trace_per_fiber_deg (default None = off) is stage 1 of the full
    per-fiber-independent trace redesign (see porting-notes.md): when
    set to an integer degree (6, matching the input PSF's own native
    trace degree, is the natural first thing to try), it replaces the
    shared low-degree trace_legendre_deg_wave_x/y basis entirely with a
    block-diagonal-by-fiber one -- each of the bundle's 25 fibers gets
    its own independent (trace_per_fiber_deg+1)-term wavelength basis
    with zero cross-fiber sharing, matching C++'s per-fiber trace
    parameter count. Experimental/opt-in -- not validated at production
    scale yet, and real GPU memory/wall-time cost has not been measured
    beyond the single bundle-0 forced-spots test in porting-notes.md.

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

    if legendre_deg_wave is None or fit_continuum is None or trace_legendre_deg_wave_x is None or trace_legendre_deg_wave_y is None:
        import fitsio
        cam = fitsio.read_header(arc_file, ext=0)['CAMERA'].strip().lower()
        band = cam[0]
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
    if broken_fibers:
        print(f"  Broken Fibers: {broken_fibers}")

    if backend == "gpu":
        n_workers = min(len(all_bundles), n_gpus * workers_per_gpu)
        cpu_threads_per_worker = None
    else:
        n_workers = min(len(all_bundles), cpu_workers or n_gpus)
        # See fit_bundle_task's thread-limiting comment for why this exists
        # at all: without it, every CPU worker independently claims all
        # available cores, and N concurrent workers thrash each other.
        # Scaling to (real core count / worker count) keeps the pool's
        # aggregate thread demand within the node's actual budget. Floor
        # of 1 thread/worker (an oversubscribed-but-not-zero fallback) if
        # there happen to be more workers than cores.
        available_cores = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
        cpu_threads_per_worker = max(1, available_cores // n_workers)
        print(f"  CPU thread budget: {available_cores} cores / {n_workers} workers = {cpu_threads_per_worker} threads/worker", flush=True)

    ctx = mp.get_context('spawn')
    with ctx.Pool(processes=n_workers) as pool:
        tasks = []
        for i, bid in enumerate(all_bundles):
            gpu_id = i % n_gpus
            # Use 2s stagger to prevent JIT compilation contention on CPU
            stagger_s = i * 2.0 if backend == "cpu" else 0.0
            tasks.append((bid, gpu_id, arc_file, in_psf_file, out_psf_file, lamp_lines_file, backend, broken_fibers, sn_threshold, h_size_y, stagger_s, force_spots_path, max_number_of_lines, None, legendre_deg_wave, fit_continuum, double_precision, None, trace_legendre_deg_wave_x, trace_legendre_deg_wave_y, trace_per_fiber_deg, cpu_threads_per_worker, line_search))

        print(f"Launching {len(tasks)} bundles across {n_workers} workers...", flush=True)
        chunk_results = pool.starmap(fit_bundle_task, tasks)
        for bid, res in chunk_results:
            if "error" in res:
                print(f"WARNING: Bundle {bid} failed: {res['error']}")
            else:
                bundle_results[bid] = res

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
    parser.add_argument("--trace-per-fiber-deg", type=int, default=None, help="EXPERIMENTAL (stage 1 of the full per-fiber trace redesign, see porting-notes.md): replaces the shared trace basis with a block-diagonal-by-fiber one at this wavelength degree (6 matches the input PSF's own native trace degree, and C++'s per-fiber parameter count). Overrides --trace-legendre-deg-wave-x/-y entirely when set. Default: off (None).")
    parser.add_argument("--fit-continuum", action=argparse.BooleanOptionalAction, default=None, help="Fit a per-bundle continuum background (default: auto, matching real C++ production -- on for z-band, off otherwise)")
    parser.add_argument("--gpu", type=int, default=4, help="Number of GPUs to use")
    parser.add_argument("--workers-per-gpu", type=int, default=4, help="Concurrent bundle-fit worker processes packed onto each GPU (validated safe ceiling: 4)")
    parser.add_argument("--cpu-workers", type=int, help="Concurrent worker processes for --backend cpu (default: --gpu count)")
    parser.add_argument("--backend", type=str, default="gpu", choices=["cpu", "gpu"])
    parser.add_argument("--broken-fibers", type=str, help="Comma-separated list of broken fibers")
    parser.add_argument("--sn-threshold", type=float, default=3.0, help="S/N threshold for spot selection")
    parser.add_argument("--max-lines", type=int, default=200, help="Maximum number of lines to keep per bundle")
    parser.add_argument("--h-size-y", type=int, default=5, help="Override PSF stamp half-size in Y")
    parser.add_argument("--force-spots", type=str, help="Path to a file containing spots to fit (fiber,wave,xc,yc)")
    parser.add_argument("--double-precision", action="store_true", help="Force full float64 precision for the joint-fit Jacobian (default: mixed float32/float64 -- see porting-notes.md; validated equivalent accuracy, ~71%% less GPU memory/worker)")
    parser.add_argument("--line-search", type=str, default="grid", choices=["grid", "brent", "cpp"], help="EXPERIMENTAL: the final joint fit's per-iteration step-size search. 'grid' (default): the long-standing coarse 3-point [0.2,0.5,1.0] search. 'brent': a continuous but NOT C++-faithful search, kept for reference. 'cpp': a faithful replica of C++'s actual algorithm (mode-dependent skip logic + a direct Numerical Recipes brent() port, see specex_psf_fitter.cc/specex_brent.cc). Both 'brent' and 'cpp' tested negative (no xrms/yrms change on 2 hard + 2 normal bundles) -- see porting-notes.md.")

    args = parser.parse_args()
    
    if not args.lamp_lines:
        base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        args.lamp_lines = os.path.join(base, "specex/data/specex_linelist_desi.txt")

    os.environ["JAX_PLATFORM_NAME"] = args.backend
    
    fit_ccd_native(
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
        legendre_deg_wave=args.legendre_deg_wave,
        trace_legendre_deg_wave=args.trace_legendre_deg_wave,
        trace_legendre_deg_wave_x=args.trace_legendre_deg_wave_x,
        trace_legendre_deg_wave_y=args.trace_legendre_deg_wave_y,
        trace_per_fiber_deg=args.trace_per_fiber_deg,
        fit_continuum=args.fit_continuum,
        double_precision=args.double_precision,
        line_search=args.line_search
    )

if __name__ == "__main__":
    main()
