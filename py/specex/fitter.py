import os
import time
import numpy as np
from functools import partial

from .psf import GaussHermitePSF

# --- Helpers ---

def _cpp_brent(f, ax, bx, cx, tol, itmax=100):
    """Direct line-by-line port of specex_brent.cc's brent() (itself a
    near-verbatim copy of Numerical Recipes' brent()), used by
    SPECEX_CPP_LINESEARCH for a genuinely faithful replica of C++'s line
    search. Not scipy's minimize_scalar(method='brent') -- that raises on a
    "loose" (non-strictly-bracketing) triple, which C++'s raw NR
    implementation tolerates by just falling back to golden-section
    stepping; this port preserves that tolerance instead of failing.
    Returns (x_min, f_min)."""
    CGOLD = 0.3819660
    ZEPS = 1e-60
    a = min(ax, cx); b = max(ax, cx)
    x = w = v = bx
    fx = fw = fv = f(x)
    d = 0.0; e = 0.0
    for _ in range(itmax):
        xm = 0.5 * (a + b)
        tol1 = tol * abs(x) + ZEPS
        tol2 = 2.0 * tol1
        if abs(x - xm) <= (tol2 - 0.5 * (b - a)):
            return x, fx
        if abs(e) > tol1:
            r = (x - w) * (fx - fv)
            q = (x - v) * (fx - fw)
            p = (x - v) * q - (x - w) * r
            q = 2.0 * (q - r)
            if q > 0.0: p = -p
            q = abs(q)
            etemp = e; e = d
            if abs(p) >= abs(0.5 * q * etemp) or p <= q * (a - x) or p >= q * (b - x):
                e = (a - x) if x >= xm else (b - x); d = CGOLD * e
            else:
                d = p / q; u = x + d
                if (u - a) < tol2 or (b - u) < tol2:
                    d = tol1 if (xm - x) > 0 else -tol1
        else:
            e = (a - x) if x >= xm else (b - x); d = CGOLD * e
        u = x + d if abs(d) >= tol1 else x + (tol1 if d > 0 else -tol1)
        fu = f(u)
        if fu <= fx:
            if u >= x: a = x
            else: b = x
            v, w, x = w, x, u; fv, fw, fx = fw, fx, fu
        else:
            if u < x: a = u
            else: b = u
            if fu <= fw or w == x:
                v, w = w, u; fv, fw = fw, fu
            elif fu <= fv or v == x or v == w:
                v = u; fv = fu
    return x, fx

def next_pow2_bucket(n, min_bucket=256):
    """Smallest power of 2 >= n (floored at min_bucket), for padding
    variable-length JAX inputs to a small, campaign-stable set of shapes so
    the persistent JIT-compilation cache (see porting-notes.md "JAX
    persistent compilation cache") gets reused across bundles/cameras
    instead of triggering a fresh XLA compile per distinct array length.
    """
    if n <= min_bucket:
        return min_bucket
    return 1 << (n - 1).bit_length()

def get_sparse_nz(xdeg, ydeg):
    nz = []
    for j in range(ydeg + 1):
        for i in range(xdeg + 1):
            if i == 0: nz.append(i + j*(xdeg + 1))
            elif i == 1 and j < 2: nz.append(i + j*(xdeg + 1))
            elif i > 1 and j == 0: nz.append(i + j*(xdeg + 1))
    return nz

def build_warm_start_pc(psf, bundle_id, spots, gh_deg, monomials):
    """
    Mirrors C++'s default behavior (specex_pyio.cc: use_input_specex_psf=True
    unless psf/trace degrees are overridden on the CLI) of warm-starting the
    PSF shape fit from the input PSF's own already-fit Gauss-Hermite
    coefficients, rather than cold-starting from a flat default. Projects the
    input PSF's per-fiber Legendre-in-wave GH coefficients into the same
    sparse 2D (fiber, wave) Legendre basis used by the joint fit, via
    least-squares.
    """
    param_mapping = ['GHSIGX', 'GHSIGY']
    for j_gh in range(gh_deg + 1):
        for i_gh in range(gh_deg + 1):
            if i_gh == 0 and j_gh == 0: continue
            param_mapping.append(f'GH-{i_gh}-{j_gh}')
    bundle = psf.params_of_bundles.get(bundle_id)
    fibers = [s['fiber'] for s in spots]; waves = [s['wave'] for s in spots]
    mon_np = np.array(monomials)
    rows = []
    for name in param_mapping:
        models = bundle.param_models.get(name) if bundle is not None else None
        if models is None:
            rows.append(np.zeros(mon_np.shape[1]))
            continue
        vals = np.array([float(models[f].value(w)) for f, w in zip(fibers, waves)])
        coeff, _, _, _ = np.linalg.lstsq(mon_np, vals, rcond=None)
        rows.append(coeff)
    return np.array(rows)

def get_bundle_monomials_jnp(psf, bundle_id, spots, wdeg=3):
    import jax.numpy as jnp
    from .math import legendre_pol_jnp
    bundle = psf.params_of_bundles[bundle_id]
    xdeg = 1; nz = get_sparse_nz(xdeg, wdeg)
    fmin, fmax = bundle.fiber_min, bundle.fiber_max
    wmin, wmax = psf.fiber_traces[fmin]['X_vs_W'].xmin, psf.fiber_traces[fmin]['X_vs_W'].xmax
    fiber_vals = jnp.array([s['fiber'] for s in spots]); wave_vals = jnp.array([s['wave'] for s in spots])
    rf = 2 * (fiber_vals - fmin) / (fmax - fmin) - 1; rw = 2 * (wave_vals - wmin) / (wmax - wmin) - 1
    mx = [legendre_pol_jnp(i, rf) for i in range(xdeg + 1)]; mw = [legendre_pol_jnp(j, rw) for j in range(wdeg + 1)]
    m = []
    for k in nz: i, j = k % (xdeg + 1), k // (xdeg + 1); m.append(mx[i] * mw[j])
    return jnp.stack(m, axis=1)

def get_bundle_block_diagonal_trace_monomials(psf, bundle_id, spots, trace_deg):
    """
    Per-fiber-independent trace design matrix (stage 1 of the full
    per-fiber trace redesign, see porting-notes.md) -- block-diagonal by
    fiber, each of the bundle's fibers getting its own (trace_deg+1)
    wavelength-Legendre columns with zero cross-fiber sharing, mirroring
    C++'s per-fiber independent Y_vs_W/X_vs_W refit (specex_psf_fitter.cc:
    1213-1238) exactly in the DOF sense. Expressed as a *correction* on
    top of xc_init/yc_init (already-close anchors from spot
    selection/--force-spots) rather than replacing the trace outright, so
    no changes are needed to PSF_Fitter.fit()'s anchor+correction
    architecture, jax jit kernel signatures, or the line-search/step
    logic -- this is a drop-in replacement for get_bundle_monomials_jnp's
    output, just built from a different (structurally sparse, densely
    stored) basis. One shared matrix serves both X and Y trace_coeffs
    (same wavelength basis, different coefficient values), same
    convention as the shared-basis path.
    """
    import jax.numpy as jnp
    from .math import legendre_pol_jnp
    bundle = psf.params_of_bundles[bundle_id]
    fmin, fmax = bundle.fiber_min, bundle.fiber_max
    n_fibers = fmax - fmin + 1
    wmin, wmax = psf.fiber_traces[fmin]['X_vs_W'].xmin, psf.fiber_traces[fmin]['X_vs_W'].xmax
    fiber_vals = jnp.array([s['fiber'] for s in spots]); wave_vals = jnp.array([s['wave'] for s in spots])
    local_fiber = (fiber_vals - fmin).astype(jnp.int32)
    rw = 2 * (wave_vals - wmin) / (wmax - wmin) - 1
    wave_mono = jnp.stack([legendre_pol_jnp(k, rw) for k in range(trace_deg + 1)], axis=1)  # (Ns, trace_deg+1)
    onehot = (jnp.arange(n_fibers)[jnp.newaxis, :] == local_fiber[:, jnp.newaxis]).astype(wave_mono.dtype)  # (Ns, n_fibers)
    block = onehot[:, :, jnp.newaxis] * wave_mono[:, jnp.newaxis, :]  # (Ns, n_fibers, trace_deg+1)
    return block.reshape(block.shape[0], n_fibers * (trace_deg + 1))


def compute_fiber_ndead(psf, fiber, weight):
    """Port of C++'s per-fiber dead-column diagnostic
    (specex_psf_fitter.cc:2333-2355, the source of its own logged
    "fiber N ndead=..." lines): counts zero-weight pixels in a +/-3-column
    window around the fiber's trace center, across the trace's full
    wavelength-defined row range. `weight` is indexed [x, y] (same
    convention as fit()'s xpix/ypix/idx_map usage). Vectorized over rows
    (a Python per-row loop with Legendre1DPol.invert()'s 1000-point grid
    interpolation inside would be ~4000 calls/fiber -- too slow to run
    routinely); only the +/-3 window is a short explicit loop.
    """
    trace = psf.fiber_traces[fiber]
    y_vs_w = trace['Y_vs_W']
    nx, ny = weight.shape
    begin_j = max(0, int(np.floor(y_vs_w.value(y_vs_w.xmin))))
    end_j = min(ny, int(np.floor(y_vs_w.value(y_vs_w.xmax))) + 1)
    if end_j <= begin_j:
        return 0
    j_arr = np.arange(begin_j, end_j)
    w_arr = y_vs_w.invert(j_arr.astype(float))
    i_center = np.round(np.asarray(psf.x_ccd(fiber, w_arr))).astype(int)
    ndead = 0
    for d in range(-3, 4):
        i_idx = i_center + d
        valid = (i_idx >= 0) & (i_idx < nx)
        ndead += int(np.sum(weight[i_idx[valid], j_arr[valid]] == 0))
    return ndead


def build_trace_prior_hessian(n_fibers, ndeg, prior_deg, weight, fiber_flag=None):
    """Port of C++'s trace-coefficient prior (specex_psf_fitter.cc:759-857,
    gated there by trace_prior_deg>0, off by default and not enabled by real
    DESI production -- see porting-notes.md's 2026-08-05 ndead investigation).
    For each wavelength-Legendre degree d >= prior_deg, C++ adds
    chi2 += weight*sum_i(c_i - mean_{j!=i}(c_j))**2 across the bundle's
    fibers -- a soft constraint pulling each fiber's HIGH-order per-fiber
    trace coefficients toward cross-fiber consensus, while leaving degrees
    below prior_deg (and, implicitly, everything when this is off) fully
    per-fiber independent. This is the mechanism C++ has, but doesn't use,
    for damping a single noisy/dead-column fiber's high-order coefficients
    under an otherwise fully-independent per-fiber trace basis (see
    z5@20220408 bundle 6 fiber 163 and z2@20220314 bundle 9 fiber 238,
    ndead=11829/2289 respectively, both driving their whole bundle's
    regression alone under --trace-per-fiber-deg with no such damping).

    fiber_flag (optional length-n_fibers 0/1 array): restricts the chi2 SUM
    above to i in flagged fibers only -- unflagged fibers never get their
    own residual/penalty term (no pull toward anything), though their
    coefficients still appear inside a flagged fiber's own "mean of the
    others" target. This matters because C++'s own literal prior, applied
    to every fiber unconditionally (fiber_flag=None, i.e. all-ones), was
    found to measurably HURT already-healthy bundles when tested here
    (b1@20260401 bundle 0: xrms 0.0087->0.0152, yrms 0.0068->0.0136 at
    C++'s weight=1e8) -- expected, since C++ never actually runs with this
    prior on in production, so that weight was never tuned against real
    per-fiber data. Gating activation by ndead (see fit()'s
    SPECEX_TRACE_PRIOR_NDEAD_THRESHOLD) makes this a no-op on the ~590/600
    bundles that don't have a bad fiber, by construction.

    Returns the (n_fibers*ndeg, n_fibers*ndeg) Hessian contribution H such
    that, for a coefficient vector c flattened as c[fiber*ndeg+d] (matching
    get_bundle_block_diagonal_trace_monomials' block.reshape layout), the
    prior's Gauss-Newton contribution is A += H, B += -H @ c (see the
    inline derivation at the fit() call site: this is J^T W J for a linear
    "residual" r = -L_F @ c pulling toward L_F @ c = 0 where L_F is L with
    only the flagged fibers' rows kept, i.e. flagged fibers' coefficients
    get pulled toward the mean of the OTHER fibers' coefficients at that
    degree, while unflagged fibers contribute no residual of their own
    (L_F^T L_F = L @ diag(fiber_flag) @ L for symmetric L, since a 0/1
    diagonal is idempotent)).
    """
    import jax.numpy as jnp
    if n_fibers < 2:
        return jnp.zeros((n_fibers * ndeg, n_fibers * ndeg))
    L = (n_fibers / (n_fibers - 1)) * jnp.eye(n_fibers) - (1.0 / (n_fibers - 1)) * jnp.ones((n_fibers, n_fibers))
    if fiber_flag is None:
        LFL = L @ L  # L symmetric, all fibers flagged: L^T L == L @ L
    else:
        flag_diag = jnp.diag(jnp.asarray(fiber_flag, dtype=L.dtype))
        LFL = L @ flag_diag @ L
    mask = jnp.array([1.0 if d >= prior_deg else 0.0 for d in range(ndeg)])
    return weight * jnp.kron(LFL, jnp.diag(mask))

# --- Global JIT Kernels (Standardized signatures) ---

def _predict_bundle_jax(flux, psf_coeffs, trace_coeffs, continuum_coeffs,
                        xc_init, yc_init, psf_monomials, trace_monomials, xpix, ypix,
                        sx_g, sy_g, idx_gg, degree,
                        tx_g, tw_g, wmin_c, wmax_c, image_data, weight_data):
    import jax.numpy as jnp
    from jax import vmap
    from .math import legendre_pol_jnp

    Ns = flux.shape[0]; Np = xpix.shape[0]
    gh_all = jnp.dot(psf_monomials, psf_coeffs.T)
    dx, dy = jnp.dot(trace_monomials, trace_coeffs[0]), jnp.dot(trace_monomials, trace_coeffs[1])
    xc_all, yc_all = xc_init + dx, yc_init + dy
    
    def spot_sig(i):
        Bx, By = GaussHermitePSF.get_gh_basis(xc_all[i], yc_all[i], sx_g[i], sy_g[i], gh_all[i], degree)
        psf_v = By[0] * Bx[0] 
        nx_p, ny_p = degree + 1, degree + 1; k = 2
        for j_p in range(ny_p):
            bj = By[j_p]; imin_p = 1 if j_p == 0 else 0
            for i_p in range(imin_p, nx_p):
                psf_v += gh_all[i, k] * bj * Bx[i_p]; k += 1
        return psf_v
    
    unit_sigs = vmap(spot_sig)(jnp.arange(Ns))
    tsig_p = jnp.zeros(Np + 1).at[idx_gg.flatten()].add((unit_sigs * flux[:, jnp.newaxis]).flatten())
    
    Ncont = continuum_coeffs.shape[0]
    rw = 2 * (tw_g - wmin_c) / (wmax_c - wmin_c) - 1
    m_c = jnp.stack([legendre_pol_jnp(k, rw) for k in range(Ncont)], axis=0)
    f_c = jnp.tensordot(continuum_coeffs, m_c, axes=([0], [0]))
    striped_cont = jnp.sum(vmap(lambda fi: f_c[fi] * jnp.exp(-0.5 * (xpix - tx_g[fi])**2) / jnp.sqrt(2 * jnp.pi))(jnp.arange(25)), axis=0)
    
    total_sig = tsig_p[:Np] + striped_cont
    return jnp.sum(weight_data * (image_data - total_sig)**2)

from jax import jit
_predict_bundle_jax_jit = jit(_predict_bundle_jax, static_argnums=(13,))

def _accumulate_bundle_jax(flux, psf_coeffs, trace_coeffs, continuum_coeffs,
                               xc_init, yc_init, psf_monomials, trace_monomials, xpix, ypix,
                               sx_g, sy_g, idx_gg, degree,
                               tx_g, tw_g, wmin_c, wmax_c, image_data, weight_data,
                               gain, psf_error, wscale):

    import jax
    import jax.numpy as jnp
    from jax import vmap, lax

    # psf_monomials and trace_monomials are deliberately separate design
    # matrices (possibly different wavelength degree/Npoly) -- see
    # PSF_Fitter.fit()'s trace_wdeg parameter and porting-notes.md "trace
    # correction basis" entries. Before this split they were one shared
    # array, which coupled the PSF-shape and trace-position fits through a
    # shared basis in a way C++ never does (its trace fit is a fully
    # independent per-fiber polynomial refit) -- raising the shared wdeg to
    # fix trace's missing wavelength curvature also handed the PSF-shape
    # terms more freedom at the same time, opening a trace-position/
    # PSF-asymmetry degeneracy that inflated xrms.
    Ns = flux.shape[0]; Npoly_psf = psf_monomials.shape[1]; Npoly_trace = trace_monomials.shape[1]
    Np = xpix.shape[0]; Ncont = continuum_coeffs.shape[0]
    Nparams = psf_coeffs.shape[0]; Nsh = Nparams * Npoly_psf + 2 * Npoly_trace; Ntot = Ns + Nsh + Ncont; stamp_area = sx_g.shape[1]

    gh_all = jnp.dot(psf_monomials, psf_coeffs.T); dx, dy = jnp.dot(trace_monomials, trace_coeffs[0]), jnp.dot(trace_monomials, trace_coeffs[1])
    xc_all, yc_all = xc_init + dx, yc_init + dy
    
    # No padding: batch_size = Ns exactly. This used to be a flat 2000
    # regardless of the real spot count (~600-1600 in cases seen so far).
    # Every input to this jitted function (flux, xc_init, sx_g, ...) already
    # has Ns baked into its shape, so a fresh XLA compile happens per
    # distinct Ns either way -- there's no compile-cache reuse the old flat
    # constant was protecting. Measured effect (single-GPU nvidia-smi
    # profiling, see porting-notes.md "task 22 follow-up"): this removes the
    # padding waste cleanly and never costs more than the old flat 2000, but
    # it did NOT reduce the observed process-level GPU memory peak for the
    # dominant ~1500-1600-spot bundle case even at this most-aggressive
    # (zero-padding) setting -- so this alone does not close the 8.65GB ->
    # <8.0GB gap needed for 5 workers/GPU. The real driver of that peak is
    # still unidentified (likely an XLA compile-time/kernel-selection
    # scratch allocation that only shrinks for much smaller problem sizes,
    # not the live tensor footprint) and needs real profiler-based
    # investigation (jax.profiler device memory profile), not more
    # code-reading guesses.
    batch_size = Ns; n_pad = batch_size - Ns
    gh_p = jnp.pad(gh_all, ((0, n_pad), (0, 0))); xc_p, yc_p, f_p = jnp.pad(xc_all, (0, n_pad)), jnp.pad(yc_all, (0, n_pad)), jnp.pad(flux, (0, n_pad))
    m_p = jnp.pad(psf_monomials, ((0, n_pad), (0, 0))); m_trace_p = jnp.pad(trace_monomials, ((0, n_pad), (0, 0)))
    sx_p = jnp.pad(sx_g, ((0, n_pad), (0, 0))); sy_p = jnp.pad(sy_g, ((0, n_pad), (0, 0)))
    idx_p = jnp.pad(idx_gg, ((0, n_pad), (0, 0)), constant_values=Np); mask_p = jnp.arange(batch_size) < Ns

    def get_all_grads(xi, yi, gi, si_x, si_y):
        sigx = jnp.maximum(gi[0], 0.1); sigy = jnp.maximum(gi[1], 0.1)
        isx = 1.0 / sigx; isy = 1.0 / sigy
        x1 = (jnp.floor(si_x + 0.5) - xi - 0.5) * isx; x2 = (jnp.floor(si_x + 0.5) - xi + 0.5) * isx
        y1 = (jnp.floor(si_y + 0.5) - yi - 0.5) * isy; y2 = (jnp.floor(si_y + 0.5) - yi + 0.5) * isy
        isq2 = 1.0 / jnp.sqrt(2.0); isq2pi = 1.0 / jnp.sqrt(2.0 * jnp.pi)
        gx1 = isq2pi * isx * jnp.exp(-0.5 * x1**2); gx2 = isq2pi * isx * jnp.exp(-0.5 * x2**2)
        gy1 = isq2pi * isy * jnp.exp(-0.5 * y1**2); gy2 = isq2pi * isy * jnp.exp(-0.5 * y2**2)
        ex = 0.5 * (jax.scipy.special.erf(x2 * isq2) - jax.scipy.special.erf(x1 * isq2))
        ey = 0.5 * (jax.scipy.special.erf(y2 * isq2) - jax.scipy.special.erf(y1 * isq2))
        
        from .psf import hermite_pol_jnp
        nx_p = degree + 1; ny_p = degree + 1
        H1_u = jnp.stack([hermite_pol_jnp(n, x1) for n in range(nx_p)], axis=0)
        H2_u = jnp.stack([hermite_pol_jnp(n, x2) for n in range(nx_p)], axis=0)
        H1_v = jnp.stack([hermite_pol_jnp(n, y1) for n in range(ny_p)], axis=0)
        H2_v = jnp.stack([hermite_pol_jnp(n, y2) for n in range(ny_p)], axis=0)
        
        def get_P(n, val, g1, g2, h1, h2, sig):
            return jnp.where(n == 0, val, sig * (g1 * h1[n-1] - g2 * h2[n-1]))
        
        Bx = vmap(lambda n: get_P(n, ex, gx1, gx2, H1_u, H2_u, sigx))(jnp.arange(nx_p))
        By = vmap(lambda n: get_P(n, ey, gy1, gy2, H1_v, H2_v, sigy))(jnp.arange(ny_p))
        dexdsx = (x1 * gx1 - x2 * gx2); deydsy = (y1 * gy1 - y2 * gy2); dexdx = (gx1 - gx2); deydy = (gy1 - gy2)
        
        def get_dPds(n, g1, g2, h1, h2, x1, x2, val, dvalds, sig):
            h1_m1 = h1[jnp.maximum(n-1, 0)]; h2_m1 = h2[jnp.maximum(n-1, 0)]
            h1_m2 = h1[jnp.maximum(n-2, 0)]; h2_m2 = h2[jnp.maximum(n-2, 0)]
            t1 = sig * g1 * h1_m1; t2 = sig * g2 * h2_m1
            h_prime1 = (n-1)*h1_m2; h_prime2 = (n-1)*h2_m2
            res = (-g1*x1*h_prime1 + g2*x2*h_prime2) + (t1*x1*x1/sig - t2*x2*x2/sig)
            return jnp.where(n == 0, dvalds, res)

        def get_dPdx(n, g1, g2, h1, h2, x1, x2, val, dvaldx, sig):
            h1_m1 = h1[jnp.maximum(n-1, 0)]; h2_m1 = h2[jnp.maximum(n-1, 0)]
            h1_m2 = h1[jnp.maximum(n-2, 0)]; h2_m2 = h2[jnp.maximum(n-2, 0)]
            t1 = sig * g1 * h1_m1; t2 = sig * g2 * h2_m1
            h_prime1 = (n-1)*h1_m2; h_prime2 = (n-1)*h2_m2
            res = -1.0/sig * (sig*g1*h_prime1 - sig*g2*h_prime2 - x1*t1 + x2*t2)
            return jnp.where(n == 0, dvaldx, res)

        dBxdsx = vmap(lambda n: get_dPds(n, gx1, gx2, H1_u, H2_u, x1, x2, ex, dexdsx, sigx))(jnp.arange(nx_p))
        dBydsy = vmap(lambda n: get_dPds(n, gy1, gy2, H1_v, H2_v, y1, y2, ey, deydsy, sigy))(jnp.arange(ny_p))
        dBxdx = vmap(lambda n: get_dPdx(n, gx1, gx2, H1_u, H2_u, x1, x2, ex, dexdx, sigx))(jnp.arange(nx_p))
        dBydy = vmap(lambda n: get_dPdx(n, gy1, gy2, H1_v, H2_v, y1, y2, ey, deydy, sigy))(jnp.arange(ny_p))

        psf_v = By[0] * Bx[0]; basis_gh = []; k = 2
        dsigx = dBxdsx[0] * By[0]; dsigy = Bx[0] * dBydsy[0]
        dxc = dBxdx[0] * By[0]; dyc = Bx[0] * dBydy[0]
        for j in range(ny_p):
            bj = By[j]; dbjdsy = dBydsy[j]; dbjdy = dBydy[j]
            imin = 1 if j == 0 else 0
            for i in range(imin, nx_p):
                bi = Bx[i]; dbidsx = dBxdsx[i]; dbidx = dBxdx[i]
                term = bj * bi; psf_v += gi[k] * term; basis_gh.append(term)
                dsigx += gi[k] * bj * dbidsx; dsigy += gi[k] * bi * dbjdsy
                dxc += gi[k] * bj * dbidx; dyc += gi[k] * bi * dbjdy
                k += 1
        return psf_v, jnp.stack(basis_gh, axis=1), dsigx, dsigy, dxc, dyc

    b_psf, b_gh_basis, b_jsigx, b_jsigy, b_jx, b_jy = vmap(get_all_grads)(xc_p, yc_p, gh_p, sx_p, sy_p)
    bm = mask_p[:, jnp.newaxis]; b_psf *= bm

    # --- Mixed precision (default on; set SPECEX_MIXED_PRECISION=0, or the
    # CLI's --double-precision, to force full float64): the (batch_size,
    # stamp_area, Nsh) Jacobian terms below are the single largest tensor in
    # this function (b_jac and its j_sx/j_sy/j_gh/j_xc/j_yc constituents) --
    # build them in float32 instead of the ambient float64 (jax_enable_x64,
    # psf.py) to roughly halve their footprint and the transient scratch XLA
    # needs to construct/contract them. Everything upstream (get_all_grads'
    # erf/exp/Hermite-recurrence math, which is numerically delicate and not
    # the memory driver) and everything downstream (the small Ntot x Ntot
    # accumulated normal-equations matrix A/B and the linear solve in
    # PSF_Fitter.fit) stays float64 -- only the wide-but-shallow per-spot
    # per-pixel per-parameter broadcast products are narrowed. Validated:
    # single-bundle chi2 relative error 2.4e-6, full-CCD wavelength RMS vs
    # line-list truth matches the float64 pipeline to 4 decimals, 71% GPU
    # memory cut (8657MiB -> 2513MiB/worker) -- see porting-notes.md.
    _mp = os.environ.get("SPECEX_MIXED_PRECISION", "1") != "0"
    _jdt = jnp.float32 if _mp else jnp.float64
    f_p_j, m_p_j, m_trace_p_j = f_p.astype(_jdt), m_p.astype(_jdt), m_trace_p.astype(_jdt)
    b_gh_basis_j = b_gh_basis.astype(_jdt); b_jsigx_j, b_jsigy_j = b_jsigx.astype(_jdt), b_jsigy.astype(_jdt)
    b_jx_j, b_jy_j = b_jx.astype(_jdt), b_jy.astype(_jdt)
    j_sx = (f_p_j[:, jnp.newaxis, jnp.newaxis] * b_jsigx_j[:, :, jnp.newaxis] * m_p_j[:, jnp.newaxis, :])
    j_sy = (f_p_j[:, jnp.newaxis, jnp.newaxis] * b_jsigy_j[:, :, jnp.newaxis] * m_p_j[:, jnp.newaxis, :])
    j_gh = (f_p_j[:, jnp.newaxis, jnp.newaxis, jnp.newaxis] * b_gh_basis_j[:, :, :, jnp.newaxis] * m_p_j[:, jnp.newaxis, jnp.newaxis, :]).reshape(batch_size, stamp_area, -1)
    j_xc = (f_p_j[:, jnp.newaxis, jnp.newaxis] * b_jx_j[:, :, jnp.newaxis] * m_trace_p_j[:, jnp.newaxis, :])
    j_yc = (f_p_j[:, jnp.newaxis, jnp.newaxis] * b_jy_j[:, :, jnp.newaxis] * m_trace_p_j[:, jnp.newaxis, :])
    b_jac = jnp.concatenate([j_sx, j_sy, j_gh, j_xc, j_yc], axis=2) * bm[:, jnp.newaxis].astype(_jdt)

    flat_idx = idx_p.flatten(); valid = flat_idx < Np
    total_sig = jnp.zeros(Np + 1).at[flat_idx].add((b_psf * f_p[:, jnp.newaxis]).flatten())
    from .math import legendre_pol_jnp
    rw = 2 * (tw_g - wmin_c) / (wmax_c - wmin_c) - 1
    m_cont = jnp.stack([legendre_pol_jnp(k, rw) for k in range(Ncont)], axis=0)
    f_cont = jnp.tensordot(continuum_coeffs, m_cont, axes=([0], [0]))
    striped_cont = jnp.sum(vmap(lambda fi: f_cont[fi] * jnp.exp(-0.5 * (xpix - tx_g[fi])**2) / jnp.sqrt(2 * jnp.pi))(jnp.arange(25)), axis=0)
    res = image_data - (total_sig[:Np] + striped_cont); chi2 = jnp.sum(weight_data * res**2)
    h_cont = vmap(lambda k: jnp.sum(vmap(lambda fi: m_cont[k, fi] * jnp.exp(-0.5 * (xpix - tx_g[fi])**2) / jnp.sqrt(2 * jnp.pi))(jnp.arange(25)), axis=0))(jnp.arange(Ncont)).T

    b_res = jnp.where(valid, res[jnp.where(valid, flat_idx, 0)], 0.0).reshape(batch_size, stamp_area)
    b_w = jnp.where(valid, weight_data[jnp.where(valid, flat_idx, 0)], 0.0).reshape(batch_size, stamp_area); wr = b_w * b_res

    # No B-vector Poisson/signal-dependent correction here (specex_psf_fitter.cc:670-673's
    # `bfact += (1/wscale)*0.5*(w*res)^2*(1/gain+2*psf_error^2*signal)`).
    # That C++ line is gated by `recompute_weight_in_fit`, which is
    # declared `false` in the PSF_Fitter constructor
    # (specex_psf_fitter.h:151) and is never once assigned `true` anywhere
    # else in the C++ codebase (confirmed by grepping the whole source
    # tree) -- structurally unreachable dead code in the real reference
    # implementation, not just "off by default." An earlier version of
    # this port applied the correction unconditionally, every mode, every
    # iteration, for this whole project -- a real, if narrow, C++/Python
    # mismatch. Tested disabling it across a 15-bundle sample
    # (2026-07-30, porting-notes.md): mean xrms 0.0428->0.0420, mean yrms
    # 0.0466->0.0458 -- a real, modest improvement concentrated in the 2
    # cases with the largest residual*signal product (the disabled term
    # scaled with (w*res)^2 * signal, so bundles with small residuals were
    # unaffected regardless of trace difficulty), no regression anywhere
    # in the sample.
    wr_corrected = wr
    
    # b_jac is (float32 if SPECEX_MIXED_PRECISION else float64) -- cast its
    # co-operands to match so these contractions don't get silently
    # upcast/downcast-mismatched or promoted back to float64 by JAX's
    # type-promotion rules, then cast each (small) result back to float64
    # before it's written into the float64 accumulated A/B.
    wr_corrected_j, b_w_j, b_psf_j = wr_corrected.astype(_jdt), b_w.astype(_jdt), b_psf.astype(_jdt)
    B = jnp.zeros(Ntot).at[:Ns].set(jnp.sum(wr * b_psf, axis=1)[:Ns]).at[Ns:Ns+Nsh].set(jnp.sum(jnp.sum(b_jac * wr_corrected_j[:, :, jnp.newaxis], axis=1), axis=0).astype(jnp.float64)).at[-Ncont:].set(jnp.dot(h_cont.T, weight_data * res))
    A = jnp.zeros((Ntot, Ntot)).at[jnp.arange(Ns), jnp.arange(Ns)].set(jnp.sum(b_w * b_psf**2, axis=1)[:Ns])
    A = A.at[Ns:Ns+Nsh, Ns:Ns+Nsh].set(jnp.einsum('bij,bi,bik->jk', b_jac, b_w_j, b_jac).astype(jnp.float64))
    A_fs = jnp.einsum('bij,bi,bi->bj', b_jac, b_w_j, b_psf_j).astype(jnp.float64); A = A.at[Ns:Ns+Nsh, :Ns].set(A_fs[:Ns].T).at[:Ns, Ns:Ns+Nsh].set(A_fs[:Ns])
    h_lookup = h_cont[jnp.where(valid, flat_idx, 0)]; b_hcont = jnp.where(valid[:, jnp.newaxis], h_lookup, 0.0).reshape(batch_size, stamp_area, Ncont)
    A_fc = jnp.einsum('bi,bi,bik->bk', b_psf, b_w, b_hcont); A = A.at[:Ns, -Ncont:].set(A_fc[:Ns]).at[-Ncont:, :Ns].set(A_fc[:Ns].T)
    A_sc = jnp.einsum('bij,bi,bik->jk', b_jac, b_w_j, b_hcont.astype(_jdt)).astype(jnp.float64); A = A.at[Ns:Ns+Nsh, -Ncont:].set(A_sc).at[-Ncont:, Ns:Ns+Nsh].set(A_sc.T)
    A = A.at[-Ncont:, -Ncont:].set(jnp.dot(h_cont.T * weight_data, h_cont))
    return chi2, A[:Ns+Nsh+Ncont, :Ns+Nsh+Ncont], B[:Ns+Nsh+Ncont]

_accumulate_bundle_jax_jit = jit(_accumulate_bundle_jax, static_argnums=(13, 20, 21, 22))

def apply_dead_column_mask(psf, fiber_min, fiber_max, weight):
    nx, ny = weight.shape
    w_new = weight.copy(); rows_j = np.arange(ny).astype(float)
    for fib in range(fiber_min, fiber_max + 1):
        try:
            w_vals = psf.fiber_traces[fib]['Y_vs_W'].invert(rows_j)
            x_vals = psf.x_ccd(fib, w_vals); i_centers = np.floor(x_vals + 0.5).astype(int)
            valid = (rows_j >= 0) & (rows_j < ny) & (i_centers >= 0) & (i_centers < nx)
            bad = weight[i_centers[valid], rows_j[valid].astype(int)] == 0
            for j in np.where(valid)[0][bad]:
                ic = i_centers[j]; i_s, i_e = max(0, ic - 4), min(nx, ic + 5)
                w_new[i_s:i_e, j] = 0.0
        except: continue
    return w_new

def filter_dead_column_spots(spots, weight):
    """
    Port of C++'s per-spot can_measure_flux/ignore exclusion
    (specex_psf_fitter.cc:244-258, InitTmpData -- gated there behind
    spots.size()>1 and confirmed, by reading FitIndividualSpotFluxes/
    FitOneSpot, to NEVER fire during the SNR-based selection stage, only
    during the final joint fit's FitSeveralSpots(selected_spots,...) calls,
    same pipeline stage as this function's caller): a spot is dropped from
    the fit entirely if more than 5 of the 25 pixels in a 5x5 window
    centered on its own position have zero weight ("can survive one dead
    column = 5pix, not more"). C++'s center is int(floor(x)+0.5), which for
    positive x truncates the .5 straight back off and is therefore exactly
    floor(x) -- NOT the round-to-nearest-pixel floor(x+0.5) convention used
    elsewhere in this file (e.g. apply_dead_column_mask) -- replicated
    faithfully here rather than substituting the other convention.

    A genuinely different mechanism from apply_dead_column_mask: that
    function broadens zero-weight regions (a +/-4px column band) affecting
    every spot whose stamp overlaps it; this one leaves the weight array
    completely untouched and instead makes a binary per-spot keep/drop
    decision using the RAW (unbroadened) weight -- so this must run BEFORE
    apply_dead_column_mask, on the original weight array.
    """
    if len(spots) <= 1:
        return spots
    nx, ny = weight.shape
    kept = []
    for s in spots:
        ic, jc = int(np.floor(s['xc_init'])), int(np.floor(s['yc_init']))
        i0, i1 = max(0, ic - 2), min(nx, ic + 3)
        j0, j1 = max(0, jc - 2), min(ny, jc + 3)
        nbad = int(np.sum(weight[i0:i1, j0:j1] == 0))
        if nbad <= 5:
            kept.append(s)
    return kept

def get_bundle_footprint(psf, spots, fiber_min, fiber_max, weight=None):
    t0 = time.time(); nx, ny = (4114, 4128)
    if weight is not None: nx, ny = weight.shape
    rows_j = np.arange(ny).astype(float)
    w1 = psf.fiber_traces[fiber_min]['Y_vs_W'].invert(rows_j); w2 = psf.fiber_traces[fiber_max]['Y_vs_W'].invert(rows_j)
    x1, x2 = psf.x_ccd(fiber_min, w1), psf.x_ccd(fiber_max, w2)
    xmin_env, xmax_env = np.floor(np.minimum(x1, x2) + 0.5).astype(int), np.floor(np.maximum(x1, x2) + 0.5).astype(int) + 1
    j_min_all = max(0, min(s['stamp_jmin'] for s in spots)); j_max_all = min(ny, max(s['stamp_jmax'] for s in spots))
    i_min_all = max(0, min(s['stamp_imin'] for s in spots)); i_max_all = min(nx, max(s['stamp_imax'] for s in spots))
    active_mask = np.zeros((i_max_all - i_min_all, j_max_all - j_min_all), dtype=bool)
    for s in spots:
        i0, i1 = max(i_min_all, s['stamp_imin']), min(i_max_all, s['stamp_imax'])
        j0, j1 = max(j_min_all, s['stamp_jmin']), min(j_max_all, s['stamp_jmax'])
        active_mask[i0 - i_min_all : i1 - i_min_all, j0 - j_min_all : j1 - j_min_all] = True
    ii, jj = np.where(active_mask); global_i, global_j = ii + i_min_all, jj + j_min_all
    valid_env = (global_i >= xmin_env[global_j]) & (global_i < xmax_env[global_j])
    if weight is not None:
        valid_weight = weight[global_i[valid_env], global_j[valid_env]] > 0
        fi, fj = global_i[valid_env][valid_weight], global_j[valid_env][valid_weight]
    else: fi, fj = global_i[valid_env], global_j[valid_env]
    pixels = np.stack([fi, fj], axis=1); pixels = pixels[np.lexsort((pixels[:, 1], pixels[:, 0]))]
    print(f"  Footprint generation took {time.time() - t0:.2f}s ({len(pixels)} pixels)", flush=True)
    return pixels[:, 0], pixels[:, 1], {(p[0], p[1]): i for i, p in enumerate(pixels)}

def select_spots_cpp(fibers, waves, snrs, xc, yc, image_shape,
                      sn_threshold, min_wave_dist, max_number_of_lines,
                      min_dwave=5.0, max_dwave=300.0):
    """
    Pure spot-selection logic mirroring C++ PSF_Fitter::select_spots
    (src/specex_psf_fitter.cc:1976-2249). Operates on parallel arrays of
    already-fit candidate stats and performs no fitting itself, so it can
    be called repeatedly with different thresholds as the fit iterates
    (matching FitEverything's multi-pass selection).

    Returns a 0/1 status array over the candidates.
    """
    fibers = np.asarray(fibers); waves = np.asarray(waves); snrs = np.asarray(snrs)
    xc = np.asarray(xc); yc = np.asarray(yc)
    nx, ny = image_shape
    Ns = len(fibers)
    status = np.zeros(Ns, dtype=int)

    # --- First pass: S/N, in-bounds, and min wavelength distance to nearest same-fiber neighbor ---
    in_bounds = (xc >= 0) & (xc < nx) & (yc >= 0) & (yc < ny)
    passes_snr = snrs >= sn_threshold
    for i in range(Ns):
        if not (passes_snr[i] and in_bounds[i]):
            continue
        if min_wave_dist > 0:
            same_fiber = fibers == fibers[i]
            same_fiber[i] = False
            dist = np.min(np.abs(waves[same_fiber] - waves[i])) if np.any(same_fiber) else 1000.0
            if dist < min_wave_dist:
                continue
        status[i] = 1

    if max_number_of_lines <= 0:
        return status

    # --- Second pass: coverage-limited line pruning (mirrors lines 2063-2190) ---
    # C++ bins wavelength with a truncating int() cast, not rounding.
    wave_ids = (waves * 10).astype(np.int64)
    selected_mask = status == 1
    unique_ids = np.unique(wave_ids[selected_mask])
    if unique_ids.size == 0:
        return status

    nspots_per_wave = {int(wid): int(np.sum(selected_mask & (wave_ids == wid))) for wid in unique_ids}
    snr_per_wave = {int(wid): float(np.mean(snrs[selected_mask & (wave_ids == wid)])) for wid in unique_ids}

    sorted_ids = sorted(nspots_per_wave.keys())
    selected_waves = {wid: 1 for wid in sorted_ids}
    max_fibers = max(nspots_per_wave.values())

    begin_id = 0
    has_found_max = False
    for wid in sorted_ids:
        if has_found_max:
            begin_id = wid
            break
        if nspots_per_wave[wid] == max_fibers:
            has_found_max = True
    end_id = 0
    for wid in reversed(sorted_ids):
        if nspots_per_wave[wid] == max_fibers:
            end_id = wid - 1
            break

    # Pass 2a: remove low-S/N lines one at a time, skipping any whose removal
    # would leave a gap wider than max_dwave between its selected neighbors.
    while True:
        num_lines = sum(1 for wid in sorted_ids
                         if selected_waves[wid] == 1 and nspots_per_wave[wid] >= max_fibers * 0.6)
        if num_lines <= max_number_of_lines:
            break

        wave_vs_snr = sorted((snr_per_wave[wid], wid) for wid in sorted_ids if selected_waves[wid] == 1)

        waveid_to_remove = None
        for _, wid in wave_vs_snr:
            if wid <= begin_id or wid >= end_id:
                continue
            previous_selected = begin_id
            next_selected = end_id
            for swid in sorted_ids:
                if selected_waves[swid] == 1:
                    if swid < wid and swid > previous_selected:
                        previous_selected = swid
                    if swid > wid and swid < next_selected:
                        next_selected = swid
            dwave = (next_selected - previous_selected) / 10.0
            if dwave > max_dwave:
                continue
            waveid_to_remove = wid
            break

        if waveid_to_remove is None:
            break
        selected_waves[waveid_to_remove] = 0

    # Pass 2b: bring back any unselected line within min_dwave of a selected one.
    while True:
        brought_back = False
        for wid_s in sorted_ids:
            if selected_waves[wid_s] == 0:
                continue
            for wid_j in sorted_ids:
                if selected_waves[wid_j] == 1:
                    continue
                if abs(wid_s - wid_j) / 10.0 < min_dwave:
                    selected_waves[wid_j] = 1
                    brought_back = True
        if not brought_back:
            break

    for i in range(Ns):
        wid = int(wave_ids[i])
        if selected_waves.get(wid, 0) == 0:
            status[i] = 0

    return status

def generate_bundle_candidates(psf, fiber_min, fiber_max, lamp_lines, image_shape, broken_fibers=None):
    """
    Build the raw spot-candidate list for a bundle from the lamp line list and
    the current trace model, mirroring the candidate list C++ builds before
    any fitting (fiber x wavelength grid, filtered to fall on the CCD).
    """
    lines_sc = [l for l in lamp_lines if 1 <= l.get('score', 1) <= 4]
    broken_list = []
    if broken_fibers:
        if isinstance(broken_fibers, str): broken_list = [int(f) for f in broken_fibers.split(",") if f.strip()]
        else: broken_list = list(broken_fibers)
    nx, ny = image_shape
    candidates = []
    for fiber in range(fiber_min, fiber_max + 1):
        if fiber in broken_list or fiber not in psf.fiber_traces: continue
        trace = psf.fiber_traces[fiber]
        x_vs_w, y_vs_w = trace['X_vs_W'], trace['Y_vs_W']
        for line in lines_sc:
            wave = line['wave']
            # Port of C++'s wavelength-domain gate (specex_lamp_lines_utils.cc:
            # 65-72, allocate_spots_of_bundle): reject any line-list entry
            # outside this fiber's OWN trace's fitted wavelength domain,
            # before ever evaluating a CCD position for it. Without this,
            # psf.x_ccd/y_ccd (a Legendre polynomial evaluation) happily
            # *extrapolates* a plausible-looking (xc,yc) for a wavelength far
            # outside where the trace was ever actually calibrated -- found
            # via a handful of real-but-out-of-band XeI lines
            # (specex_linelist_desi.txt's ~9660-9802A high-order-diffraction-
            # ghost entries, "added by hand... from KPNO data inspection")
            # landing inside a b-band CCD's pixel bounds by extrapolation
            # coincidence and corrupting --trace-per-fiber-deg's per-fiber
            # polynomial fit for every fiber in the bundle, not just the 1-3
            # fibers the spurious candidate actually appears on (see
            # porting-notes.md's 2026-08-05 b2@20260401/b6@20250125
            # investigation). C++ never has this problem because this check
            # runs unconditionally, independent of --trace-per-fiber-deg.
            if wave < x_vs_w.xmin or wave > x_vs_w.xmax: continue
            if wave < y_vs_w.xmin or wave > y_vs_w.xmax: continue
            xc = psf.x_ccd(fiber, wave); yc = psf.y_ccd(fiber, wave)
            if 0 <= xc < nx and -4 <= yc < ny + 4:
                candidates.append({'fiber': fiber, 'wave': wave, 'xc_init': xc, 'yc_init': yc})
    return candidates


def fit_candidate_fluxes(psf, candidates, image, weight):
    """
    Individual-spot flux fit for every candidate, holding position fixed
    (mirrors C++ FitIndividualSpotFluxes: fit_flux=true, fit_position=false).
    Uses the housekeeping-phase stamp cap hSizeX/Y=min(3, ...)
    (specex_psf_fitter.cc:2500-2501). Updates each candidate dict in place
    with flux/eflux/snr/chi2 and returns the parallel numpy arrays.
    """
    import jax.numpy as jnp
    c_xc = jnp.array([s['xc_init'] for s in candidates]); c_yc = jnp.array([s['yc_init'] for s in candidates])
    gh_all = jnp.array([psf.gh_params(s['fiber'], s['wave']) for s in candidates])
    housekeeping_hsize_x = min(3, psf.h_size_x); housekeeping_hsize_y = min(3, psf.h_size_y)
    fluxes, snrs, chi2s, efluxes = _get_spot_stats_jax(jnp.array(image), jnp.array(weight), c_xc, c_yc, gh_all,
                                                        psf.gh_psf.degree, housekeeping_hsize_x, housekeeping_hsize_y)
    fluxes, snrs, chi2s, efluxes = np.array(fluxes), np.array(snrs), np.array(chi2s), np.array(efluxes)
    for i, s in enumerate(candidates):
        s['flux'] = float(fluxes[i]); s['eflux'] = float(efluxes[i])
        s['snr'] = float(snrs[i]); s['chi2'] = float(chi2s[i])
    return fluxes, efluxes, snrs, chi2s


def _finalize_selected(psf, candidates, status):
    selected = []
    for s, keep in zip(candidates, status):
        if keep:
            s['stamp_imin'] = int(np.floor(s['xc_init'] + 0.5)) - psf.h_size_x
            s['stamp_imax'] = int(np.floor(s['xc_init'] + 0.5)) + psf.h_size_x + 1
            s['stamp_jmin'] = int(np.floor(s['yc_init'] + 0.5)) - psf.h_size_y
            s['stamp_jmax'] = int(np.floor(s['yc_init'] + 0.5)) + psf.h_size_y + 1
            selected.append(s)
    return selected


def select_bundle_spots_iterative(psf, fiber_min, fiber_max, lamp_lines, image, weight, bundle_id,
                                   broken_fibers=None, max_number_of_lines=200, wdeg=3, fit_continuum=True):
    """
    Mirrors the housekeeping/selection phase of C++ FitEverything
    (specex_psf_fitter.cc:2493-2783): repeatedly fits individual candidate
    fluxes and reselects from the *full* raw candidate list, interleaved
    with a trace-only warm-up fit that updates every candidate's xc/yc from
    the fitted trace model, before a final loose-threshold selection pass.

    min_snr_non_linear_terms=5, min_wave_dist_non_linear_terms=4A ("strict")
    and min_snr_linear_terms=3, min_wave_dist_linear_terms=0 ("loose") match
    the C++ constants of the same name.
    """
    t0 = time.time()
    nx, ny = image.shape
    candidates = generate_bundle_candidates(psf, fiber_min, fiber_max, lamp_lines, (nx, ny), broken_fibers)
    print(f"PYTHON SPOT SELECTION: Initial candidates: {len(candidates)}", flush=True)
    if not candidates:
        return []

    fibers_arr = np.array([s['fiber'] for s in candidates])
    waves_arr = np.array([s['wave'] for s in candidates])

    def strict_select():
        fluxes, efluxes, snrs, chi2s = fit_candidate_fluxes(psf, candidates, image, weight)
        xc = np.array([s['xc_init'] for s in candidates]); yc = np.array([s['yc_init'] for s in candidates])
        status = select_spots_cpp(fibers_arr, waves_arr, snrs, xc, yc, (nx, ny),
                                   sn_threshold=5.0, min_wave_dist=4.0, max_number_of_lines=max_number_of_lines)
        return status

    if hasattr(psf, 'output_psf_path') and psf.output_psf_path:
        raw_path = psf.output_psf_path.replace('.fits', '.pyrawspots.txt')

    def write_pass_checkpoint(status, pass_name):
        if not (hasattr(psf, 'output_psf_path') and psf.output_psf_path):
            return
        path = psf.output_psf_path.replace('.fits', f'.pyspots_{pass_name}.txt')
        with open(path, 'w') as f:
            for s, keep in zip(candidates, status):
                if keep:
                    f.write(f"{s['fiber']},{s['wave']:.7f},{s['xc_init']:.7f},{s['yc_init']:.7f}\n")
        print(f"  Written {int(status.sum())} spots to {path}", flush=True)

    # Pass 1: individual flux fit + strict select on raw (un-refit) candidates.
    # This is the direct analog of C++'s cpp_cp0_pass1.txt/cppspots_pass1.txt -
    # first FitIndividualSpotFluxes call, before any trace warm-up.
    status = strict_select()
    print(f"  Pass 1 (strict): {int(status.sum())}/{len(candidates)} spots", flush=True)
    write_pass_checkpoint(status, 'pass1')
    if hasattr(psf, 'output_psf_path') and psf.output_psf_path:
        with open(raw_path, 'w') as f:
            for s in candidates:
                f.write(f"{s['fiber']},{s['wave']:.7f},{s['xc_init']:.7f},{s['yc_init']:.7f},{s['flux']:.7f},{s['eflux']:.7f},{s['snr']:.7f}\n")
        print(f"  Written {len(candidates)} raw (pass-1) candidates to {raw_path}", flush=True)

    # Trace-only warm-up loop: refit fluxes, strict-reselect, fit FLUX+TRACE on
    # the selected set, then snap every candidate's xc/yc to the new trace
    # model. Break once the max centroid shift is small, matching C++'s
    # trace_loop (up to 5 iterations, break when max_delta < 0.5px).
    for trace_loop in range(5):
        status = strict_select()
        if trace_loop == 0:
            write_pass_checkpoint(status, 'pass2')
        selected = _finalize_selected(psf, candidates, status)
        if not selected:
            break
        fitter = PSF_Fitter(psf)
        _, pc, tc, cc, flux_fit, xc_final, yc_final, _ = fitter.fit(image, weight, selected, bundle_id, max_iter=5, wdeg=wdeg, fit_continuum=fit_continuum)

        max_delta = 0.0
        for s in candidates:
            nxc = psf.x_ccd(s['fiber'], s['wave'], tc_x=tc[0], wdeg=wdeg)
            nyc = psf.y_ccd(s['fiber'], s['wave'], tc_y=tc[1], wdeg=wdeg)
            max_delta = max(max_delta, ((s['xc_init'] - nxc)**2 + (s['yc_init'] - nyc)**2) ** 0.5)
            s['xc_init'] = float(nxc); s['yc_init'] = float(nyc)
        print(f"  Trace warm-up {trace_loop}: max centroid shift = {max_delta:.4f}px", flush=True)
        if max_delta < 0.5:
            break

    # One more strict pass on the trace-refined candidates.
    status = strict_select()
    print(f"  Pass 3 (strict, post-trace): {int(status.sum())}/{len(candidates)} spots", flush=True)
    write_pass_checkpoint(status, 'pass3')

    # Final pass: loose thresholds (min_snr_linear_terms=3, min_wave_dist=0) -
    # this is the list C++ hands to the final joint PSF+FLUX fit.
    fluxes, efluxes, snrs, chi2s = fit_candidate_fluxes(psf, candidates, image, weight)

    if hasattr(psf, 'output_psf_path') and psf.output_psf_path:
        final_raw_path = psf.output_psf_path.replace('.fits', '.pyrawspots_final.txt')
        with open(final_raw_path, 'w') as f:
            for i, s in enumerate(candidates):
                f.write(f"{s['fiber']},{s['wave']:.7f},{s['xc_init']:.7f},{s['yc_init']:.7f},{fluxes[i]:.7f},{efluxes[i]:.7f},{snrs[i]:.7f}\n")
        print(f"  Written {len(candidates)} trace-refined candidates to {final_raw_path}", flush=True)

    xc = np.array([s['xc_init'] for s in candidates]); yc = np.array([s['yc_init'] for s in candidates])
    status = select_spots_cpp(fibers_arr, waves_arr, snrs, xc, yc, (nx, ny),
                               sn_threshold=3.0, min_wave_dist=0.0, max_number_of_lines=max_number_of_lines)
    selected = _finalize_selected(psf, candidates, status)

    if hasattr(psf, 'output_psf_path') and psf.output_psf_path:
        spots_path = psf.output_psf_path.replace('.fits', '.pyspots.txt')
        with open(spots_path, 'w') as f:
            for s in selected:
                f.write(f"{s['fiber']},{s['wave']:.7f},{s['xc_init']:.7f},{s['yc_init']:.7f}\n")
        print(f"  Written {len(selected)} spots to {spots_path}", flush=True)

    print(f"  Iterative spot selection took {time.time() - t0:.2f}s ({len(selected)} spots)", flush=True)
    return selected


def get_bundle_spots(psf, fiber_min, fiber_max, lamp_lines, image=None, weight=None,
                       sn_threshold=3.0, min_dist_angstrom=0.0, wave_min=None, wave_max=None,
                       broken_fibers=None, max_number_of_lines=100, provided_candidates=None):
    """
    Single-pass candidate generation + selection (no trace warm-up / re-fit
    loop). Kept as a lighter-weight building block; prefer
    select_bundle_spots_iterative for full C++ FitEverything parity.
    """
    t0 = time.time()
    nx, ny = (4114, 4128)
    if image is not None: nx, ny = image.shape

    if provided_candidates is not None:
        candidates = provided_candidates
    else:
        candidates = generate_bundle_candidates(psf, fiber_min, fiber_max, lamp_lines, (nx, ny), broken_fibers)

    print(f"PYTHON SPOT SELECTION: Starting first pass. Initial candidates: {len(candidates)}", flush=True)

    if not candidates: return []

    Ns = len(candidates)
    if image is not None:
        fluxes, efluxes, snrs, chi2s = fit_candidate_fluxes(psf, candidates, image, weight)
        waves = np.array([s['wave'] for s in candidates])

        if hasattr(psf, 'output_psf_path') and psf.output_psf_path:
            raw_path = psf.output_psf_path.replace('.fits', '.pyrawspots.txt')
            with open(raw_path, 'w') as f:
                for i in range(Ns):
                    s = candidates[i]
                    f.write(f"{s['fiber']},{s['wave']:.7f},{s['xc_init']:.7f},{s['yc_init']:.7f},{fluxes[i]:.7f},{efluxes[i]:.7f},{snrs[i]:.7f}\n")
            print(f"  Written {len(candidates)} raw candidates to {raw_path}", flush=True)

        fibers_arr = np.array([s['fiber'] for s in candidates])
        xc = np.array([s['xc_init'] for s in candidates]); yc = np.array([s['yc_init'] for s in candidates])
        status = select_spots_cpp(fibers_arr, waves, snrs, xc, yc, (nx, ny),
                                   sn_threshold=sn_threshold, min_wave_dist=min_dist_angstrom,
                                   max_number_of_lines=max_number_of_lines)
        print(f"  Selection: {int(np.sum(status))}/{Ns} spots survive "
              f"(SNR>={sn_threshold}, min_dwave={min_dist_angstrom}A, max_lines={max_number_of_lines})", flush=True)

        selected = _finalize_selected(psf, candidates, status)

        if hasattr(psf, 'output_psf_path') and psf.output_psf_path:
            spots_path = psf.output_psf_path.replace('.fits', '.pyspots.txt')
            with open(spots_path, 'w') as f:
                for s in selected:
                    f.write(f"{s['fiber']},{s['wave']:.7f},{s['xc_init']:.7f},{s['yc_init']:.7f}\n")
            print(f"  Written {len(selected)} spots to {spots_path}", flush=True)

        print(f"  Spot selection took {time.time() - t0:.2f}s ({len(selected)} spots)", flush=True)
        return selected
    return candidates


def _fit_one_spot_jax(image, weight, xc, yc, gh, degree, hsize_x, hsize_y):
    import jax.numpy as jnp
    from jax import jit
    nx, ny = image.shape
    ix_rel, iy_rel = jnp.meshgrid(jnp.arange(2*hsize_x+1), jnp.arange(2*hsize_y+1), indexing='ij')
    dx = ix_rel.flatten() - hsize_x; dy = iy_rel.flatten() - hsize_y
    
    def fit_loop(flux, x_shift, y_shift):
        # Simple 3-parameter fit (flux, dx, dy) for selection purposes
        # This mimics C++ FitOneSpot's initial phase
        im = jnp.floor(xc + 0.5).astype(int); jm = jnp.floor(yc + 0.5).astype(int)
        gix = im + dx + x_shift; giy = jm + dy + y_shift
        valid = (gix >= 0) & (gix < nx) & (giy >= 0) & (giy < ny)
        p_val = GaussHermitePSF.single_pix_value_jnp(xc + x_shift, yc + y_shift, gix, giy, gh, degree)
        p_val = jnp.where(valid, p_val, 0.0)
        d_val = image[jnp.clip(gix, 0, nx-1), jnp.clip(giy, 0, ny-1)]
        w_val = weight[jnp.clip(gix, 0, nx-1), jnp.clip(giy, 0, ny-1)]
        w_val = jnp.where(valid, w_val, 0.0)
        
        A = jnp.sum(w_val * p_val**2)
        B = jnp.sum(w_val * d_val * p_val)
        flux_est = jnp.where(A > 0, B/A, 0.0)
        return flux_est, A

    return jit(fit_loop)(0.0, 0.0, 0.0)

def _fit_all_spots_batch(cand_xc, cand_yc, gh_params, image, weight, hsize_x, hsize_y, degree):
    """Per-candidate closed-form flux/S/N fit (position held fixed, mirrors
    C++ FitIndividualSpotFluxes), vmapped over the candidate batch. `image`/
    `weight` are real traced arguments (not closures) so this module-level
    jitted function is cached purely by (shape, dtype), same as
    _accumulate_bundle_jax_jit/_predict_bundle_jax_jit -- letting it be
    reused across different cameras/exposures that share the same CCD shape,
    not just repeated calls with the identical image object.
    """
    import jax.numpy as jnp
    from jax import vmap
    nx, ny = image.shape
    ix_rel, iy_rel = jnp.meshgrid(jnp.arange(2*hsize_x+1), jnp.arange(2*hsize_y+1), indexing='ij')
    dx = ix_rel.flatten() - hsize_x; dy = iy_rel.flatten() - hsize_y

    def fit_spot(xc, yc, gh):
        im = jnp.floor(xc + 0.5).astype(int); jm = jnp.floor(yc + 0.5).astype(int)
        gix = im + dx; giy = jm + dy
        valid = (gix >= 0) & (gix < nx) & (giy >= 0) & (giy < ny)
        p_val = GaussHermitePSF.single_pix_value_jnp(xc, yc, gix, giy, gh, degree)
        p_val = jnp.where(valid, p_val, 0.0)
        d_val = image[jnp.clip(gix, 0, nx-1), jnp.clip(giy, 0, ny-1)]
        w_val = weight[jnp.clip(gix, 0, nx-1), jnp.clip(giy, 0, ny-1)]
        w_val = jnp.where(valid, w_val, 0.0)

        A = jnp.sum(w_val * p_val**2)
        B = jnp.sum(w_val * d_val * p_val)
        flux = jnp.where(A > 0, B/A, 0.0)
        # C++ parity: eflux = sqrt(cov) where cov = 1/A is the inverse-Hessian
        # diagonal for the flux parameter (specex_psf_fitter.cc:1721-1724).
        eflux = jnp.where(A > 0, 1.0 / jnp.sqrt(A), 0.0)
        snr = jnp.where(A > 0, flux / eflux, -1.0)
        chi2 = jnp.sum(w_val * (d_val - flux * p_val)**2)
        return flux, snr, chi2, eflux

    return vmap(fit_spot)(cand_xc, cand_yc, gh_params)

from jax import jit as _jit
_fit_all_spots_batch_jit = _jit(_fit_all_spots_batch, static_argnums=(5, 6, 7))

def _get_spot_stats_jax(image, weight, cand_xc, cand_yc, gh_params, degree, hsize_x, hsize_y):
    import jax.numpy as jnp
    # Pad the candidate batch to a power-of-2 bucket so JAX reuses one
    # compiled shape across bundles/cameras instead of recompiling for every
    # distinct raw-candidate count (campaign-wide these range ~1100-1800 and
    # collapse into a single bucket -- see porting-notes.md "power-of-2
    # shape bucketing"). Safe by construction: fit_spot is vmapped, so each
    # padding row is fit fully independently and cannot influence any real
    # candidate's result -- padding rows are simply sliced off below.
    Ns_real = cand_xc.shape[0]
    Ns_pad = next_pow2_bucket(Ns_real)
    if Ns_pad > Ns_real:
        n_extra = Ns_pad - Ns_real
        cand_xc = jnp.pad(cand_xc, (0, n_extra))
        cand_yc = jnp.pad(cand_yc, (0, n_extra))
        gh_params = jnp.pad(gh_params, ((0, n_extra), (0, 0)))

    flux, snr, chi2, eflux = _fit_all_spots_batch_jit(cand_xc, cand_yc, gh_params, image, weight, hsize_x, hsize_y, degree)
    return flux[:Ns_real], snr[:Ns_real], chi2[:Ns_real], eflux[:Ns_real]


class PSF_Fitter:
    def __init__(self, psf):
        self.psf = psf
        # EXPERIMENT (SPECEX_CHI2_PRECISION_OVERRIDE): C++'s equivalent
        # per-stage convergence threshold is a looser 0.1 (specex_psf_fitter.cc:
        # 2710,2790) vs this 0.01 -- direction argues Python already does
        # *more* refinement, not less, so unlikely to explain worse
        # hard-bundle accuracy, but never empirically swept (python-vs-cpp-
        # diff.txt section 1f). Opt-in override for that sweep.
        self.chi2_precision = float(os.environ.get("SPECEX_CHI2_PRECISION_OVERRIDE", 0.01))
    def fit(self, image, weight, spots, bundle_id, fit_type='full', max_iter=20, wdeg=3, fit_continuum=True, trace_wdeg=None, trace_wdeg_x=None, trace_wdeg_y=None, trace_per_fiber_deg=None, line_search='grid', trace_prior_deg=None):
        import jax.numpy as jnp
        print(f"Starting HIGH-PERFORMANCE OPTIMIZED fit for bundle {bundle_id}...")
        # EXPERIMENT (SPECEX_TRACE_PRIOR_DEG / SPECEX_TRACE_PRIOR_WEIGHT):
        # opt-in port of C++'s trace prior (see build_trace_prior_hessian's
        # docstring) -- only meaningful alongside trace_per_fiber_deg.
        # Env-var default follows this branch's established convention
        # (SPECEX_MATCH_CPP_DEAD_COLUMN, SPECEX_TRACE_MAX_ITERS_OVERRIDE)
        # for experimental knobs not yet promoted to a first-class CLI flag.
        if trace_prior_deg is None and os.environ.get("SPECEX_TRACE_PRIOR_DEG"):
            trace_prior_deg = int(os.environ["SPECEX_TRACE_PRIOR_DEG"])
        # Default 1e5, not C++'s literal 1e8 -- the weight sweep run
        # alongside the ndead-gating work (porting-notes.md, 2026-08-05)
        # found 1e8 measurably over-smooths even the fibers it's meant to
        # fix relative to 1e5/1e6 (more collateral pull on a bad bundle's
        # healthy fibers, no extra benefit to the bad fiber itself), and
        # 1e5 fully resolves z2@20220314 bundle 9's fiber 238 while giving
        # z5@20220408 bundle 6's far-more-extreme fiber 163 (ndead=11829,
        # largely a real data-floor problem no reweighting fixes) the same
        # benefit as any higher weight tested.
        trace_prior_weight = float(os.environ.get("SPECEX_TRACE_PRIOR_WEIGHT", 1e5))
        # trace_wdeg (a shared X/Y default) falls back to wdeg -- see
        # porting-notes.md's r2@20250109 investigation. trace_wdeg_x/
        # trace_wdeg_y independently override it per axis, falling back to
        # trace_wdeg in turn -- added after finding X didn't need (and was
        # mildly destabilized by) the same extra wavelength curvature that
        # closed the Y gap: giving X unused extra freedom let its own
        # (mostly noise-driven) high-order coefficients drift.
        trace_wdeg = wdeg if trace_wdeg is None else trace_wdeg
        trace_wdeg_x = trace_wdeg if trace_wdeg_x is None else trace_wdeg_x
        trace_wdeg_y = trace_wdeg if trace_wdeg_y is None else trace_wdeg_y
        # min()/max() over all spots, not spots[0]/spots[-1] -- the incoming
        # list is only reliably fiber-sorted when it comes from this
        # process's own select_bundle_spots_iterative(); a --force-spots
        # file (e.g. C++'s pass4 spot dump) is ordered by selection pass,
        # not by fiber, so spots[0]/[-1] silently narrowed the inferred
        # bundle range (e.g. [10,24] instead of the true [0,24]). That
        # narrowed range then clips get_bundle_footprint's trace envelope,
        # dropping the excluded fibers' own pixels from the fit entirely --
        # root cause of the r2@20250109 footprint-size anomaly in
        # porting-notes.md (Python: 66506 px vs C++: 121360 px for the same
        # forced spot list; fixed here gives 113242, matching C++ to ~7%,
        # the residual being get_bundle_footprint's own missing x-margin).
        fmin, fmax = min(s['fiber'] for s in spots), max(s['fiber'] for s in spots)
        # EXPERIMENT (SPECEX_MATCH_CPP_DEAD_COLUMN): C++'s per-spot
        # can_measure_flux/ignore exclusion (see filter_dead_column_spots'
        # docstring) -- must run on the RAW weight, before
        # apply_dead_column_mask broadens it. Opt-in for now, same pattern
        # as SPECEX_MATCH_CPP_FLUX_CLAMP, pending a correctness check on the
        # known hard bundles (porting-notes.md, 2026-08-04).
        if os.environ.get("SPECEX_MATCH_CPP_DEAD_COLUMN"):
            n_before = len(spots)
            spots = filter_dead_column_spots(spots, weight)
            if len(spots) < n_before:
                print(f"  SPECEX_MATCH_CPP_DEAD_COLUMN: dropped {n_before - len(spots)}/{n_before} spots (dead-column can_measure_flux)", flush=True)
        weight = apply_dead_column_mask(self.psf, fmin, fmax, weight)
        # EXPERIMENT (SPECEX_TRACE_PRIOR_DEG): decide, once per bundle,
        # which fibers (if any) get the trace prior activated, using C++'s
        # own per-fiber ndead diagnostic (see compute_fiber_ndead) against
        # a threshold. Computed here (using the bundle's TRUE fixed fiber
        # span, not the spots-derived fmin/fmax which shrinks for broken
        # fibers) so it's available once for every 'trace'-mode iteration
        # below rather than recomputed per-iteration -- ndead only depends
        # on `weight` and trace geometry, both fixed for this bundle.
        # Threshold default (500) matches the only precedent for this
        # exact diagnostic anywhere in the C++ codebase (its own
        # number_of_fibers_with_dead_columns gate, specex_psf_fitter.cc:
        # 2354) -- comfortably above the ~20-120 seen on this bundle's
        # normal fibers and comfortably below both known bad cases
        # (fiber163=11829, fiber238=2289).
        trace_prior_fiber_flag = None
        if trace_per_fiber_deg is not None and trace_prior_deg is not None:
            _bundle_span0 = self.psf.params_of_bundles[bundle_id]
            _fmin0, _fmax0 = _bundle_span0.fiber_min, _bundle_span0.fiber_max
            _ndead_threshold = int(os.environ.get("SPECEX_TRACE_PRIOR_NDEAD_THRESHOLD", 500))
            _ndeads = [compute_fiber_ndead(self.psf, f, weight) for f in range(_fmin0, _fmax0 + 1)]
            trace_prior_fiber_flag = np.array([1.0 if nd > _ndead_threshold else 0.0 for nd in _ndeads])
            if trace_prior_fiber_flag.any():
                _flagged = [(_fmin0 + i, nd) for i, nd in enumerate(_ndeads) if nd > _ndead_threshold]
                print(f"  SPECEX_TRACE_PRIOR_DEG: activating trace prior for {len(_flagged)} fiber(s) "
                      f"(ndead>{_ndead_threshold}): {_flagged}", flush=True)
        xpix, ypix, pix_idx = get_bundle_footprint(self.psf, spots, fmin, fmax, weight)
        Np = len(xpix); area = (2*self.psf.h_size_x+1)*(2*self.psf.h_size_y+1); Ns = len(spots)
        # Pixel-footprint padding: pad the pixel dimension fed to the JIT
        # calls below to a power-of-2 bucket (real footprints range ~33k-
        # 125k pixels campaign-wide and collapse into 2 buckets -- see
        # porting-notes.md "power-of-2 shape bucketing"). idx_g's
        # out-of-footprint sentinel must equal Np_pad (not Np): the JIT
        # functions gate "valid" as flat_idx < Np, computed internally from
        # the padded array's own (Np_pad) length, so the sentinel has to
        # track it exactly or genuinely-invalid entries would be
        # misclassified as valid.
        Np_pad = next_pow2_bucket(Np)
        t0 = time.time(); sx, sy = np.zeros((Ns, area)), np.zeros((Ns, area)); idx_g = np.full((Ns, area), Np_pad, dtype=np.int32)
        nx, ny = image.shape; idx_map = np.full((nx, ny), -1, dtype=np.int32); idx_map[xpix, ypix] = np.arange(Np)
        for s_i, s in enumerate(spots):
            ix, iy = np.meshgrid(np.arange(s['stamp_imin'], s['stamp_imax']), np.arange(s['stamp_jmin'], s['stamp_jmax']), indexing='ij')
            sx[s_i], sy[s_i] = ix.flatten(), iy.flatten()
            # Spot stamps near the CCD edge (stamp_imin/imax computed as xc_init +/- h_size_x
            # with no clamping) can extend outside [0,nx)x[0,ny) -- clip only for the idx_map
            # lookup and mask those pixels out (same -1-sentinel convention as in-bounds pixels
            # that fall outside any spot's footprint), rather than indexing idx_map out of bounds.
            in_bounds = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
            ix_c, iy_c = np.clip(ix, 0, nx - 1), np.clip(iy, 0, ny - 1)
            st_idx = idx_map[ix_c.astype(int), iy_c.astype(int)].flatten()
            mask = in_bounds.flatten() & (st_idx >= 0)
            idx_g[s_i, mask] = st_idx[mask]
        print(f"  Stamp indexing took {time.time() - t0:.2f}s", flush=True)
        rows_u = np.unique(ypix); row_m = {j: i for i, j in enumerate(rows_u)}
        tx_j, tw_j = np.zeros((25, len(rows_u))), np.zeros((25, len(rows_u)))
        for f_i in range(25):
            fib = fmin + f_i; w_v = self.psf.fiber_traces[fib]['Y_vs_W'].invert(rows_u.astype(float)); tw_j[f_i], tx_j[f_i] = w_v, self.psf.x_ccd(fib, w_v)
        ix_r = np.array([row_m[j] for j in ypix])
        # Extend xpix/ypix/ix_r to Np_pad by repeating pixel 0's coordinates.
        # This reuses an already-real row, so rows_u (and tx_g/tw_g's shape)
        # is unaffected, and no spot's idx_g ever points into the appended
        # range (idx_map only ever maps the real xpix/ypix positions above)
        # -- so the padding entries are only ever touched by the whole-array
        # image_data/weight_data/chi2 terms in the JIT functions, never by
        # any per-spot Jacobian/Hessian accumulation. w_d is forced to 0 on
        # the padding entries below so those terms are always exactly zero
        # regardless of what image value ends up duplicated there.
        n_pix_extra = Np_pad - Np
        if n_pix_extra > 0:
            xpix_p = np.concatenate([xpix, np.full(n_pix_extra, xpix[0], dtype=xpix.dtype)])
            ypix_p = np.concatenate([ypix, np.full(n_pix_extra, ypix[0], dtype=ypix.dtype)])
            ix_r_p = np.concatenate([ix_r, np.full(n_pix_extra, ix_r[0], dtype=ix_r.dtype)])
        else:
            xpix_p, ypix_p, ix_r_p = xpix, ypix, ix_r
        tx_g, tw_g = jnp.array(tx_j[:, ix_r_p]), jnp.array(tw_j[:, ix_r_p])
        flux = jnp.array([s['flux'] for s in spots]); xc_init, yc_init = jnp.array([s['xc_init'] for s in spots]), jnp.array([s['yc_init'] for s in spots])
        psf_monomials = get_bundle_monomials_jnp(self.psf, bundle_id, spots, wdeg=wdeg)
        # trace_per_fiber_deg (stage 1 of the full per-fiber redesign, see
        # porting-notes.md) swaps the shared low-degree basis for a
        # block-diagonal-by-fiber one -- each fiber gets its own
        # (trace_per_fiber_deg+1) wavelength columns with zero cross-fiber
        # sharing, matching C++'s per-fiber independent trace refit's
        # degrees of freedom. When set, it takes over entirely from
        # trace_wdeg_x/trace_wdeg_y (no freeze-masking -- both axes get
        # the full per-fiber basis).
        if trace_per_fiber_deg is not None:
            trace_monomials = get_bundle_block_diagonal_trace_monomials(self.psf, bundle_id, spots, trace_per_fiber_deg)
            Npoly_trace_x = Npoly_trace_y = trace_monomials.shape[1]
        else:
            # Single shared trace_monomials sized at max(trace_wdeg_x,
            # trace_wdeg_y): get_sparse_nz(1, d)'s output is a strict
            # prefix of get_sparse_nz(1, d+1)'s (each higher wdeg only
            # ever *appends* terms), so a lower-degree axis's basis is
            # exactly the leading columns of this shared matrix. tc's
            # unused trailing columns for whichever axis has the smaller
            # degree are frozen at zero (see the "freeze" step below,
            # right after each Newton step is computed) instead of
            # building a second differently-sized monomials array --
            # avoids threading a third matrix through the JIT kernels for
            # what's structurally just a masked subset of the one already
            # there.
            trace_wdeg_shared = max(trace_wdeg_x, trace_wdeg_y)
            trace_monomials = get_bundle_monomials_jnp(self.psf, bundle_id, spots, wdeg=trace_wdeg_shared)
            Npoly_trace_x = len(get_sparse_nz(1, trace_wdeg_x)); Npoly_trace_y = len(get_sparse_nz(1, trace_wdeg_y))
        gh_deg = self.psf.gh_psf.degree; n_gh = (gh_deg + 1) * (gh_deg + 1) - 1
        pc = jnp.array(build_warm_start_pc(self.psf, bundle_id, spots, gh_deg, psf_monomials))
        # pc0/asym_gh_rows support the anti-drift correction below: a
        # near-null Hessian direction mixes the trace correction with the
        # GH-i-0/GH-0-j ("pure x"/"pure y", i.e. antisymmetric-in-one-axis)
        # shape terms (see porting-notes.md's b-band anomaly writeup --
        # correlation of -0.93 between the leading trace_x coefficient and
        # GH-1-0). C++ never encounters this at all (it never solves trace
        # and PSF shape jointly), so there's nothing to port structurally;
        # instead we damp those specific rows back toward their warm-start
        # value every 'full'-mode iteration below, which bounds how far
        # this degenerate combination can drift over many iterations
        # without constraining any well-determined shape parameter.
        pc0 = pc
        _r = 2; _gh_row_of = {}
        for _j in range(gh_deg + 1):
            for _i in range(gh_deg + 1):
                if _i == 0 and _j == 0: continue
                _gh_row_of[(_i, _j)] = _r; _r += 1
        asym_gh_rows = jnp.array([_gh_row_of[(_i, 0)] for _i in range(1, gh_deg + 1)] +
                                  [_gh_row_of[(0, _j)] for _j in range(1, gh_deg + 1)], dtype=jnp.int32)
        tc = jnp.zeros((2, trace_monomials.shape[1])); Ncont = 4; cc = jnp.zeros(Ncont)
        img_d, w_d = jnp.array(image[xpix_p, ypix_p]), jnp.array(weight[xpix_p, ypix_p])
        if n_pix_extra > 0:
            w_d = w_d.at[Np:].set(0.0)
        xpix_j, ypix_j = jnp.array(xpix_p), jnp.array(ypix_p)
        wmin_c, wmax_c = float(self.psf.fiber_traces[fmin]['X_vs_W'].xmin), float(self.psf.fiber_traces[fmin]['X_vs_W'].xmax); old_chi2 = 1e30; prev_mode = None
        # REVERTED (branch experiment/cpp-alternating-solve): tried making
        # trace mode's exit convergence-based instead of a fixed 3
        # iterations (up to a 20-iteration cap), hypothesizing trace was
        # under-converged and that was causing the ~3% wrms (truth
        # comparison) regression seen after excluding trace from 'full'
        # mode. Result on a 9-case check: xrms/yrms/wrms were unchanged to
        # 3-4 decimal places from the fixed-3 version, case by case --
        # trace was already fully converged in 3 iterations for every case
        # tested, so the extra iterations (typically 8-11, one case ~20)
        # bought nothing. Cost ~1.5x more wall time for zero benefit.
        # Reverted cleanly; the wrms gap is not an iteration-budget problem
        # -- see porting-notes.md for the next hypothesis to test instead
        # (the unreplicated stricter-SNR trace-specific selection pass
        # C++ uses, not iteration count).

        best_chi2 = 1e30
        best_tc = tc.copy()
        best_pc = pc.copy()
        best_cc = cc.copy()
        best_flux = flux.copy()
        
        sx_g, sy_g, idx_gg = jnp.array(sx), jnp.array(sy), jnp.array(idx_g)
        # Stage/mode is now tracked as mutable state across iterations
        # (rather than derived purely from `i`) so 'trace' mode's exit can
        # be convergence-based specifically when trace_per_fiber_deg is set.
        # The REVERTED note above already showed convergence-based exit is a
        # pure-cost no-benefit change for the shared basis (9-14 params,
        # always converges within the fixed 3 iterations) -- this scopes the
        # extension to only the ~350-param per-fiber basis (2026-08-04
        # regression finding, porting-notes.md), which does NOT reliably
        # converge from zero-init in 3 iterations. Non-per-fiber runs keep
        # the exact fixed-3-iteration schedule (min_it == max_it == 3 below).
        mode = 'flux'
        stage_iter = 0
        # EXPERIMENT (SPECEX_TRACE_MAX_ITERS_OVERRIDE): quick knob for
        # testing whether a tradeoff case (e.g. r9@20220120:17, section 1a/
        # 4e) is still iteration-budget-limited at the default cap of 20,
        # opt-in, default unchanged.
        trace_min_iters, trace_max_iters = 3, int(os.environ.get("SPECEX_TRACE_MAX_ITERS_OVERRIDE", 20))
        for i in range(max_iter):
            chi2, A, B = _accumulate_bundle_jax_jit(flux, pc, tc, cc, xc_init, yc_init, psf_monomials, trace_monomials, xpix_j, ypix_j, sx_g, sy_g, idx_gg, gh_deg, tx_g, tw_g, wmin_c, wmax_c, img_d, w_d, self.psf.gain, self.psf.psf_error, 1.0)

            if i == 0 and os.environ.get("SPECEX_DEBUG_MEM"):
                import jax
                Nparams_dbg = pc.shape[0]; Npoly_psf_dbg = psf_monomials.shape[1]; Npoly_trace_dbg = trace_monomials.shape[1]; stamp_area_dbg = sx_g.shape[1]
                Nsh_dbg = Nparams_dbg * Npoly_psf_dbg + 2 * Npoly_trace_dbg
                print(f"  DEBUG_MEM: Ns={flux.shape[0]} Nparams={Nparams_dbg} Npoly_psf={Npoly_psf_dbg} Npoly_trace={Npoly_trace_dbg} "
                      f"stamp_area={stamp_area_dbg} Nsh={Nsh_dbg} Np(footprint)={xpix.shape[0]}", flush=True)
                arrs = sorted(jax.live_arrays(), key=lambda a: -a.nbytes)
                total = sum(a.nbytes for a in arrs)
                print(f"  DEBUG_MEM: {len(arrs)} live arrays, {total/1e9:.3f} GB total", flush=True)
                for a in arrs[:25]:
                    print(f"  DEBUG_MEM:   shape={a.shape} dtype={a.dtype} nbytes={a.nbytes/1e6:.2f}MB", flush=True)

            if chi2 < best_chi2:
                best_chi2 = chi2
                best_tc = tc.copy()
                best_pc = pc.copy()
                best_cc = cc.copy()
                best_flux = flux.copy()
                
            # 'sigma' mode (i=5..7, before 'full'): matches C++'s
            # scheduled_fit_of_sigmas stage exactly -- GHSIGX/GHSIGY (pc
            # rows 0-1) are fit alone (with flux), then PERMANENTLY
            # EXCLUDED from 'full' mode's idx below. This is a genuine
            # freeze, not just a warm-start: C++'s FitParPolXW for its own
            # main PSF-fit stage (specex_psf_fitter.cc:2870-2880) explicitly
            # filters out GHSIGX/GHSIGY/GHNSIG/tail terms by name, so they
            # are *never* free at the same time as the higher-order GH
            # terms in any single least-squares solve. A prior attempt
            # (see porting-notes.md, same date) only warm-started
            # GHSIGX/GHSIGY from a separate strict-selected pre-fit but
            # then let 'full' mode re-solve them jointly with everything
            # else anyway -- a coefficient-level comparison against cached
            # C++ output (also this date) showed why that couldn't work:
            # GHSIGX/GHSIGY and GH-2-0/GH-0-2 are a classic near-degenerate
            # pair (width vs. 2nd-order shape term), and Python's joint
            # solve was resolving that degeneracy differently from C++ in
            # every single one of 9 test cases (GH-2-0 larger in Python in
            # all 9; GHSIGX smaller in 8/9) regardless of warm start --
            # only *structurally excluding* sigma from the joint solve, not
            # just seeding it well, can match C++ here.
            print(f"Iter {i}: chi2 = {float(chi2):.4f} [Mode: {mode}]", flush=True)
            Npoly_psf = psf_monomials.shape[1]; Npoly_trace = trace_monomials.shape[1]; Ns_l = len(flux); n_psf_tot = (n_gh + 2) * Npoly_psf
            # EXPERIMENT (SPECEX_TRACE_PRIOR_DEG): add C++'s trace-prior
            # Gauss-Newton contribution directly to A/B (see
            # build_trace_prior_hessian's docstring for the derivation).
            # Only meaningful during 'trace' mode -- trace params are
            # excluded from idx (frozen) in every other mode on this
            # branch, so A/B changes in those rows/cols are inert there;
            # gating on mode=='trace' just avoids the wasted compute.
            # residual convention: r = -L@c (target: L@c == 0, i.e. each
            # fiber's coefficient at degree>=trace_prior_deg equals the
            # mean of the other fibers' coefficients there), so
            # A_prior = J^T W J = weight*L^T L and B_prior = J^T W r =
            # -A_prior @ c -- the sign is deliberately opposite the more
            # common "B pulls toward zero" pattern: this prior pulls
            # toward cross-fiber *consensus*, not toward zero.
            # trace_prior_fiber_flag is None whenever trace_prior_deg isn't
            # set at all, and all-zeros whenever no fiber in this bundle
            # cleared the ndead threshold -- skip entirely in either case
            # (not just for wasted compute: this is what makes the prior a
            # guaranteed no-op on every bundle without a flagged fiber).
            if (trace_per_fiber_deg is not None and trace_prior_deg is not None and mode == 'trace'
                    and trace_prior_fiber_flag is not None and trace_prior_fiber_flag.any()):
                # Bundle's true fixed fiber span (always e.g. 25), NOT
                # fmax-fmin+1 over `spots` -- that range shrinks whenever a
                # fiber has zero selected spots (broken fibers, e.g. z5
                # bundle 6's 164-174), which silently mismatched
                # get_bundle_block_diagonal_trace_monomials' own n_fibers
                # (built from the bundle's fixed fiber_min/fiber_max) and
                # broke the shapes of trace_monomials vs H_prior.
                _bundle_span = self.psf.params_of_bundles[bundle_id]
                n_fibers_bundle = _bundle_span.fiber_max - _bundle_span.fiber_min + 1
                ndeg = trace_per_fiber_deg + 1
                H_prior = build_trace_prior_hessian(n_fibers_bundle, ndeg, trace_prior_deg, trace_prior_weight,
                                                     fiber_flag=trace_prior_fiber_flag)
                trace_start = Ns_l + n_psf_tot
                x_sl = slice(trace_start, trace_start + Npoly_trace)
                y_sl = slice(trace_start + Npoly_trace, trace_start + 2 * Npoly_trace)
                A = A.at[x_sl, x_sl].add(H_prior).at[y_sl, y_sl].add(H_prior)
                B = B.at[x_sl].add(-H_prior @ tc[0]).at[y_sl].add(-H_prior @ tc[1])
            if mode == 'full' and os.environ.get("SPECEX_DEBUG_DUMP_A"):
                np.savez(os.environ["SPECEX_DEBUG_DUMP_A"], A=np.array(A), B=np.array(B),
                         Ns_l=Ns_l, n_gh=n_gh, Npoly_psf=Npoly_psf, Npoly_trace=Npoly_trace,
                         Ncont=Ncont, iter=i, chi2=float(chi2))
            # EXPERIMENT (SPECEX_FREEZE_GH10): b2@20241105:17's coefficient-
            # level disagreement is concentrated in GHSIGX and the low-order
            # pure-X shape terms (GH-1-0, GH-2-0, GH-3-0), not spread evenly
            # -- exactly the signature the existing GHSIGX/GHSIGY freeze
            # already fixed for a Y-axis-adjacent degeneracy (2026-07-29).
            # GH-1-0 (pc row 2, the lowest-order antisymmetric-in-X term,
            # near-degenerate with a small X-trace-position shift to first
            # order) has no analogous early, isolated stage -- it's only
            # ever fit jointly with every other shape term in 'full' mode.
            # This extends the existing sigma/full split to also cover row
            # 2: 'sigma' mode gains GH-1-0 alongside GHSIGX/GHSIGY, 'full'
            # mode then excludes it (frozen), same treatment.
            _gh10_extra = 1 if os.environ.get("SPECEX_FREEZE_GH10") else 0
            if mode == 'flux': idx = jnp.concatenate([jnp.arange(Ns_l), jnp.arange(A.shape[0]-Ncont, A.shape[0])])
            elif mode == 'trace': idx = jnp.concatenate([jnp.arange(Ns_l), jnp.arange(Ns_l + n_psf_tot, Ns_l + n_psf_tot + 2*Npoly_trace), jnp.arange(A.shape[0]-Ncont, A.shape[0])])
            elif mode == 'sigma': idx = jnp.concatenate([jnp.arange(Ns_l + (2+_gh10_extra)*Npoly_psf), jnp.arange(A.shape[0]-Ncont, A.shape[0])])
            # EXPERIMENT (branch experiment/cpp-alternating-solve): 'full'
            # mode excludes trace entirely instead of solving trace+shape
            # jointly -- matches C++'s real structure exactly (it never
            # solves trace and PSF shape together at all; trace is finalized
            # in its own stage and never revisited -- see porting-notes.md's
            # b-band degeneracy investigation). Trace is fit only during
            # 'trace' mode (i=2..4) and stays frozen for the remainder of
            # the fit, same as C++ freezing it after its own TRACE stage.
            # Also excludes GHSIGX/GHSIGY (pc rows 0-1, frozen after 'sigma'
            # mode above) for the same reason, plus GH-1-0 (row 2) too when
            # SPECEX_FREEZE_GH10 is set.
            else: idx = jnp.concatenate([jnp.arange(Ns_l), jnp.arange(Ns_l + (2+_gh10_extra)*Npoly_psf, Ns_l + n_psf_tot), jnp.arange(A.shape[0]-Ncont, A.shape[0])])
            A_sub, B_sub = A[jnp.ix_(idx, idx)], B[idx]; diag = jnp.diag(A_sub); S = jnp.sqrt(diag); S = jnp.where(S < 1e-12, 1.0, S)
            A_reg = (A_sub / jnp.outer(S, S)) + 1e-8 * jnp.eye(A_sub.shape[0])
            try: ds = jnp.linalg.solve(A_reg, B_sub / S); d_p = jnp.zeros(A.shape[0]).at[idx].set(ds / S)
            except: d_p = jnp.zeros(A.shape[0])
            # Freeze tc's trailing columns for whichever trace axis has the
            # smaller degree (see trace_monomials' construction above) --
            # tc starts at exactly zero and never gets a nonzero step in
            # these entries, so it stays exactly zero for the rest of the
            # fit, equivalent to that axis never having had those extra
            # basis columns at all.
            if Npoly_trace_x < Npoly_trace or Npoly_trace_y < Npoly_trace:
                trace_start = Ns_l + n_psf_tot
                d_trace = d_p[trace_start:trace_start + 2 * Npoly_trace].reshape(2, Npoly_trace)
                d_trace = d_trace.at[0, Npoly_trace_x:].set(0.0).at[1, Npoly_trace_y:].set(0.0)
                d_p = d_p.at[trace_start:trace_start + 2 * Npoly_trace].set(d_trace.flatten())
            # Ncont stays structurally 4 regardless of fit_continuum (changing
            # it to 0 would make every "-Ncont:"-style slice above/below turn
            # into "-0:", which numpy/jax silently reinterpret as the whole
            # array, not zero elements). Instead, when continuum fitting is
            # disabled, simply never let the step move cc away from its zero
            # init -- equivalent to not fitting a continuum at all.
            if not fit_continuum: d_p = d_p.at[-Ncont:].set(0.0)
            # Cap the trace-position step implied by this iteration's raw
            # Newton step at 0.5px, mirroring C++'s own per-iteration trace
            # step limiter (specex_psf_fitter.cc:1516-1533, "don't want a
            # step larger than N pix for all spots"). EXPERIMENT (branch
            # experiment/cpp-alternating-solve): restricted to 'trace' mode
            # only -- 'full' mode's d_p now has zero trace entries by
            # construction (trace excluded from idx above), so this was
            # already a guaranteed no-op there; narrowed for clarity, not a
            # behavior change from that.
            if mode == 'trace':
                trace_start = Ns_l + n_psf_tot
                d_trace_step = d_p[trace_start:trace_start + 2 * Npoly_trace].reshape(2, Npoly_trace)
                dx_spots = trace_monomials @ d_trace_step[0]
                dy_spots = trace_monomials @ d_trace_step[1]
                max_dist = jnp.max(jnp.sqrt(dx_spots ** 2 + dy_spots ** 2))
                max_trace_step = 0.5
                trace_scale = jnp.where(max_dist > max_trace_step, max_trace_step / max_dist, 1.0)
                d_p = d_p * trace_scale
            # Line search. NOTE: ls_chi2 must be a separate variable from
            # best_chi2 (the global best-state tracker above). They were
            # previously the same variable, so after each step best_chi2 held
            # the new params' predicted chi2 and the next iteration's
            # top-of-loop check (accumulate-chi2 of the SAME params) passed or
            # failed on pure float reduction-order noise between the two
            # kernels - when it failed every iteration, best_tc/best_pc froze
            # at their INITIAL values and the whole converged fit was
            # silently discarded (seen as dx/dy_final == 0 with a perfectly
            # healthy chi2 trajectory).
            # EXPERIMENT (SPECEX_MATCH_CPP_FLUX_CLAMP): C++'s real
            # FitEverything (specex_psf_fitter.cc) only sets
            # force_positive_flux=true starting at the sigma-fitting stage
            # (line ~2791, right before "PSF+FLUX only gaussian terms") --
            # during the earlier flux/trace stages (force_positive_flux
            # defaults false, never set otherwise in the real
            # direct_simultaneous_fit=true production path) C++ allows flux
            # to go transiently negative. This branch's fit() has always
            # clamped flux >= 0 on every iteration unconditionally; gated
            # here to test whether relaxing that during 'flux'/'trace'
            # modes changes hard-bundle behavior.
            clamp_flux = not (os.environ.get("SPECEX_MATCH_CPP_FLUX_CLAMP") and mode in ('flux', 'trace'))
            best_alpha, ls_chi2 = 0.0, float(chi2)
            if jnp.any(d_p != 0):
                def _ls_chi2_at(alpha):
                    f_try = flux + alpha * d_p[:Ns_l]
                    if clamp_flux: f_try = jnp.maximum(f_try, 0.0)
                    p_try = pc + alpha * d_p[Ns_l : Ns_l + n_psf_tot].reshape(n_gh + 2, Npoly_psf); t_try = tc + alpha * d_p[Ns_l + n_psf_tot : Ns_l + n_psf_tot + 2*Npoly_trace].reshape(2, Npoly_trace); c_try = cc + alpha * d_p[-Ncont:]
                    return float(_predict_bundle_jax_jit(f_try, p_try, t_try, c_try, xc_init, yc_init, psf_monomials, trace_monomials, xpix_j, ypix_j, sx_g, sy_g, idx_gg, gh_deg, tx_g, tw_g, wmin_c, wmax_c, img_d, w_d))
                if line_search == 'cpp':
                    # EXPERIMENT (--line-search cpp): a *faithful* replica of
                    # C++'s actual line-search logic
                    # (specex_psf_fitter.cc:1497-1642, specex_brent.cc --
                    # literally Numerical Recipes' brent()), not just "a
                    # continuous search instead of a grid" (that weaker
                    # version, --line-search brent below, tested negative and
                    # used the wrong bracket/method entirely).
                    # C++'s real behavior, mapped onto this branch's modes:
                    #  - 'flux' mode is a pure linear least-squares problem in
                    #    (flux, continuum) alone -- matches C++'s `linear`
                    #    condition (fit_flux && !fit_position && !fit_trace
                    #    && !fit_psf) exactly. C++ takes the raw Newton step
                    #    with NO search at all in this case.
                    #  - 'trace' mode: C++ always uses brent when fit_trace is
                    #    true, regardless of whether the raw step already
                    #    helps ("we use brent anyway for the fit of traces").
                    #  - 'sigma'/'full' modes: not linear (flux is jointly fit
                    #    with genuinely nonlinear/bilinear terms) and not
                    #    fitting trace, so C++ first checks whether the raw
                    #    Newton step (alpha=1) already decreases chi2; if so
                    #    it's taken directly with no search at all; only if it
                    #    increases chi2 does brent get called.
                    # Brent bracket/tolerance match C++ exactly: bracketing
                    # triple (-0.05, 1, 1.001) -- note this allows a small
                    # *negative* step but essentially never exceeds 1.0 -- and
                    # brent_precision=0.01 (C++ uses this as both the chi2
                    # convergence threshold AND brent's own tol argument).
                    # Uses _cpp_brent (a direct NR port, see its own
                    # docstring) rather than scipy.optimize.minimize_scalar --
                    # scipy raises on a "loose" bracket where f(xb) isn't
                    # strictly below f(xa)/f(xc), which happens routinely
                    # here (this codepath only runs when the raw step alpha=1
                    # already made chi2 *worse*, so bx=1 need not be the best
                    # of the three points at all) -- C++'s raw NR
                    # implementation has no such precondition check and just
                    # proceeds via golden-section fallback, which is what
                    # _cpp_brent replicates.
                    if mode == 'flux':
                        best_alpha, ls_chi2 = 1.0, _ls_chi2_at(1.0)
                    elif mode == 'trace':
                        best_alpha, ls_chi2 = _cpp_brent(_ls_chi2_at, -0.05, 1.0, 1.001, 0.01)
                    else:
                        chi2_1 = _ls_chi2_at(1.0)
                        if chi2_1 <= ls_chi2:
                            best_alpha, ls_chi2 = 1.0, chi2_1
                        else:
                            a_x, f_x = _cpp_brent(_ls_chi2_at, -0.05, 1.0, 1.001, 0.01)
                            if f_x < ls_chi2: best_alpha, ls_chi2 = a_x, f_x
                elif line_search == 'brent':
                    # EXPERIMENT (--line-search brent, earlier/weaker version,
                    # kept for reference): a continuous but NOT C++-faithful
                    # search -- wrong bracket/method (bounded [0,1.5] instead
                    # of C++'s real (-0.05,1,1.001) triple), wrong tolerance
                    # semantics. Tested negative (no change vs the 3-point
                    # grid on 4 cases spanning hard/normal). Superseded by
                    # --line-search cpp above for an honest comparison.
                    from scipy.optimize import minimize_scalar
                    sol = minimize_scalar(_ls_chi2_at, bounds=(0.0, 1.5), method='bounded', options={'xatol': 1e-4})
                    if sol.fun < ls_chi2: best_alpha, ls_chi2 = float(sol.x), float(sol.fun)
                else:
                    for alpha in [0.2, 0.5, 1.0]:
                        c2 = _ls_chi2_at(alpha)
                        if c2 < ls_chi2: best_alpha, ls_chi2 = alpha, c2
            if os.environ.get("SPECEX_DEBUG_ALPHA"):
                print(f"  ALPHA_DEBUG iter={i} mode={mode} best_alpha={float(best_alpha)} ls_chi2={float(ls_chi2)}", flush=True)
            # A line-search failure (no tried alpha improves chi2) only means
            # "no progress in THIS mode's active parameter subset" -- flux/
            # trace/sigma/full each free a different subset (see idx above),
            # so a failure in e.g. 'sigma' mode says nothing about whether
            # 'full' mode's still-untouched higher-order GH shape terms have
            # room to improve. Restricting the early-exit to 'full' mode
            # (the terminal stage, where nothing new unlocks afterward) keeps
            # the "genuinely converged, stop" behavior while fixing a real
            # bug where ~half of all bundles were exiting mid-'sigma' mode
            # and never fitting the higher-order shape terms at all.
            if best_alpha == 0 and mode == 'full' and i > 5: break
            if best_alpha == 0: best_alpha = 0.1
            flux = flux + best_alpha * d_p[:Ns_l]
            if clamp_flux: flux = jnp.maximum(flux, 0.0)
            pc = pc + best_alpha * d_p[Ns_l : Ns_l + n_psf_tot].reshape(n_gh + 2, Npoly_psf); tc = tc + best_alpha * d_p[Ns_l + n_psf_tot : Ns_l + n_psf_tot + 2*Npoly_trace].reshape(2, Npoly_trace); cc = cc + best_alpha * d_p[-Ncont:]
            # Anti-drift damping for the trace/GH degenerate direction --
            # DISABLED on this branch (experiment/cpp-alternating-solve).
            # The degeneracy this guards against requires trace and the
            # antisymmetric GH-shape rows to be free *simultaneously*;
            # since 'full' mode no longer includes trace in idx at all,
            # that direction is structurally unreachable here, matching
            # C++ exactly, and damping would just be inert dead weight.
            # Left commented (not deleted) so a diff against the main
            # branch stays easy to read. See main branch / porting-notes.md
            # for the live version of this fix.
            # if mode == 'full':
            #     pc = pc.at[asym_gh_rows].set(pc0[asym_gh_rows] + 0.9 * (pc[asym_gh_rows] - pc0[asym_gh_rows]))

            # --- Iterative Snapping: Update xc_init/yc_init to the current model prediction ---
            # We use a staged approach to prevent oscillation.
            # Flux mode: No snapping. Trace mode: Full snapping. Full mode: Increased frequency.
            if mode == 'trace':
                # Removed iterative snapping to align with C++ logic.
                # xc_init and yc_init must remain constant anchors.
                pass
            
            # REMOVED: All iterative snapping in 'full' mode. 
            # We keep xc_init fixed from the end of trace mode.
            # tc will now naturally accumulate the total shift from this anchor.
            
            # prev_mode == 'full' (not just mode == 'full') -- old_chi2 is
            # whatever the PREVIOUS iteration measured *before its own
            # step*, so at the very first full-mode iteration old_chi2 is
            # actually trace-mode's pre-step chi2, and this check would be
            # judging trace-mode's last step's improvement, not the
            # just-applied full-mode step's (which unlocks PSF-shape
            # parameters for the first time and deserves its own iteration
            # to be judged on). If trace mode had already nearly stalled by
            # the time full mode starts -- common, since 'trace' mode's own
            # 3 iterations are often enough to mostly converge trace given
            # how few parameters it has -- this let the loop declare
            # "converged" after a single, barely-evaluated full-mode step.
            # Confirmed: r2@20250109 bundle 16 stopped after 6 iterations at
            # a *worse* chi2 than trace_wdeg=1's 13-iteration result: see
            # porting-notes.md.
            if mode == 'full' and prev_mode == 'full' and jnp.abs(old_chi2 - chi2) < self.chi2_precision: break

            # Stage advancement (flux -> trace -> sigma -> full). stage_iter
            # counts iterations completed so far in the CURRENT mode. 'trace'
            # mode's exit is convergence-based (min 3 / max 20 iterations)
            # only when trace_per_fiber_deg is set -- see the comment above
            # the loop for why this is scoped that way (shared-basis case
            # already confirmed to need no more than a fixed 3, tested and
            # reverted). All other stage transitions keep their original
            # fixed iteration counts (2 for flux, 3 for sigma) unchanged.
            # mode_used (not the possibly-just-advanced `mode`) is what
            # prev_mode must record below -- it's this iteration's mode, the
            # one every check above (idx selection, line search, the
            # full-mode break just above) actually ran under.
            mode_used = mode
            stage_iter += 1
            if mode == 'flux' and stage_iter >= 2:
                mode, stage_iter = 'trace', 0
            elif mode == 'trace':
                per_fiber = trace_per_fiber_deg is not None
                max_it = trace_max_iters if per_fiber else trace_min_iters
                trace_stalled = (per_fiber and prev_mode == 'trace' and
                                  stage_iter >= trace_min_iters and
                                  jnp.abs(old_chi2 - chi2) < self.chi2_precision)
                if stage_iter >= max_it or trace_stalled:
                    mode, stage_iter = 'sigma', 0
            elif mode == 'sigma' and stage_iter >= 3:
                mode, stage_iter = 'full', 0

            old_chi2 = chi2; prev_mode = mode_used

        # The state after the last applied step is never seen by the
        # top-of-loop best-state check - evaluate it explicitly so a
        # converged final state cannot lose to a stale intermediate one.
        final_chi2 = float(_predict_bundle_jax_jit(flux, pc, tc, cc, xc_init, yc_init, psf_monomials, trace_monomials, xpix_j, ypix_j, sx_g, sy_g, idx_gg, gh_deg, tx_g, tw_g, wmin_c, wmax_c, img_d, w_d))
        if final_chi2 < best_chi2:
            best_chi2 = final_chi2; best_tc = tc.copy(); best_pc = pc.copy(); best_cc = cc.copy(); best_flux = flux.copy()

        # --- C++ Parity: Snap centroids to the final optimized model ---
        # Use the best coefficients found during the optimization process
        import jax.numpy as jnp
        dx_final = jnp.dot(trace_monomials, best_tc[0])
        dy_final = jnp.dot(trace_monomials, best_tc[1])
        
        # DEBUG: Check if the shifts are actually non-zero
        print(f"  DEBUG: dx_final mean={np.mean(np.abs(dx_final)):.6f}, dy_final mean={np.mean(np.abs(dy_final)):.6f}", flush=True)
        
        # The final position is the anchor (xc_init) plus the optimized shift
        # xc_init remained constant throughout the fit, mirroring C++ logic
        xc_final = np.array(xc_init + dx_final)
        yc_final = np.array(yc_init + dy_final)
        
        # Return `spots` too (not just the fitted arrays): when
        # SPECEX_MATCH_CPP_DEAD_COLUMN drops spots at the top of this
        # function, that reassignment is local to fit() and never
        # propagates back to the caller's own `spots` list -- xc_final/
        # yc_final are sized to THIS (possibly-filtered) list, so any
        # caller-side per-spot computation (x_orig/y_orig, etc.) must use
        # this returned list, not its own original one, or the array
        # lengths silently mismatch. Found 2026-08-05 via a real crash
        # (ValueError: operands could not be broadcast together with
        # shapes (1130,) (1131,)) the first time this filter actually
        # dropped a spot in a real production-scale run.
        return float(best_chi2), np.array(best_pc), np.array(best_tc), np.array(best_cc), np.array(best_flux), xc_final, yc_final, spots
