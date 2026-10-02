import numpy as np
import jax.numpy as jnp
from jax import jit, vmap, config, jacfwd
from functools import partial
from jax.scipy.special import erf as jax_erf
from .math import SparseLegendre2DPol, Legendre1DPol, hermite_pol_jnp, hermite_pol_np

config.update("jax_enable_x64", True)

class GaussHermitePSF:
    """The Gauss-Hermite PSF model: a 2D Gaussian core times a Gauss-Hermite polynomial expansion (plus an optional power-law tail), matching C++'s GaussHermitePSF (specex_gauss_hermite_psf.cc). Pixel-value evaluation is pixel-integrated (erf-based), not point-sampled, matching C++'s own convention.

    Attributes:
        degree (int): Gauss-Hermite expansion degree (both axes).

    Status: ACTIVE (production default path) -- the only PSF model this port
    implements.
    """
    def __init__(self, degree=6):
        """Args:
            degree (int): Gauss-Hermite expansion degree, default 6.

        Status: ACTIVE (production default path).
        """
        self.degree = degree

    @staticmethod
    def single_pix_value_jnp(xc, yc, x, y, params, degree):
        """Evaluate the pixel-integrated Gauss-Hermite PSF value at one pixel (x, y) for a spot centered at (xc, yc), JAX-backed (jnp), differentiable/jittable.

        The 0th-order (i=0, j=0) term is the implicit unit-amplitude core (ex*ey),
        never stored as a fit parameter -- see PSF.canonical_param_names. If
        `params` carries 5 extra trailing entries beyond the GH terms
        (TAILAMP/TAILCORE/TAILXSCA/TAILYSCA/TAILINDE), an additional power-law
        tail term is added.

        Args:
            xc, yc (float or jnp.ndarray): spot center, CCD pixel coordinates.
            x, y (float or jnp.ndarray): pixel coordinate(s) to evaluate at.
            params (jnp.ndarray): shape (2 + n_gh_terms[+5],) -- [sigma_x, sigma_y,
                GH coefficients in canonical_param_names order, optional tail
                params].
            degree (int): Gauss-Hermite expansion degree.

        Returns:
            Same shape as x/y: the PSF value at each pixel.

        Status: ACTIVE (production default path) -- the core PSF evaluation used
        throughout the fit (via pix_value_jnp) and by _fit_one_spot_jax (itself
        dead, see fitter.py).
        """
        from jax.scipy.special import erf as jax_erf
        sx = jnp.maximum(params[0], 0.1); sy = jnp.maximum(params[1], 0.1)
        isx = 1.0 / sx; isy = 1.0 / sy
        x1 = (jnp.floor(x + 0.5) - xc - 0.5) * isx; x2 = (jnp.floor(x + 0.5) - xc + 0.5) * isx
        y1 = (jnp.floor(y + 0.5) - yc - 0.5) * isy; y2 = (jnp.floor(y + 0.5) - yc + 0.5) * isy
        isq2 = 1.0 / jnp.sqrt(2.0); isq2pi = 1.0 / jnp.sqrt(2.0 * jnp.pi)
        gx1 = isq2pi * isx * jnp.exp(-0.5 * x1**2); gx2 = isq2pi * isx * jnp.exp(-0.5 * x2**2)
        gy1 = isq2pi * isy * jnp.exp(-0.5 * y1**2); gy2 = isq2pi * isy * jnp.exp(-0.5 * y2**2)
        ex = 0.5 * (jax_erf(x2 * isq2) - jax_erf(x1 * isq2)); ey = 0.5 * (jax_erf(y2 * isq2) - jax_erf(y1 * isq2))
        psfval = ex * ey; param_index = 2
        for j in range(degree + 1):
            pj = ey if j == 0 else sy * (gy1 * hermite_pol_jnp(j-1, y1) - gy2 * hermite_pol_jnp(j-1, y2))
            imin = 1 if j == 0 else 0
            for i in range(imin, degree + 1):
                pi = ex if i == 0 else sx * (gx1 * hermite_pol_jnp(i-1, x1) - gx2 * hermite_pol_jnp(i-1, x2))
                psfval += params[param_index] * pj * pi; param_index += 1
        if params.shape[0] > param_index:
            t_amp, t_core, t_xsca, t_ysca, t_inde = params[param_index:param_index+5]
            r2 = ((x-xc)*t_xsca)**2 + ((y-yc)*t_ysca)**2; denom = t_core**2 + r2 + 1e-10
            psfval += t_amp * r2 / denom * (denom)**(-t_inde/2.0)
        return psfval

    @staticmethod
    def get_gh_basis(xi, yi, sx, sy, params, degree):
        """Compute the per-order pixel-integrated 1D Gauss-Hermite basis values (Bx, By) at (xi, yi) for a spot centered at (sx, sy), so the full 2D PSF value / Jacobian can be assembled as an outer-product sum by the caller instead of re-evaluating single_pix_value_jnp per parameter.

        Args:
            xi, yi (jnp.ndarray): pixel coordinates to evaluate at.
            sx, sy (float or jnp.ndarray): spot center, CCD pixel coordinates.
            params (jnp.ndarray): shape (>=2,) -- at least [sigma_x, sigma_y].
            degree (int): Gauss-Hermite expansion degree.

        Returns:
            tuple[jnp.ndarray, jnp.ndarray]: (Bx, By), each shape
            (degree+1, *xi.shape) -- the per-order 1D basis values along X and Y.

        Status: ACTIVE (production default path) -- used by fitter.py's bundle
        Jacobian assembly (_accumulate_bundle_jax/_predict_bundle_jax).
        """
        from jax.scipy.special import erf as jax_erf
        sigx = jnp.maximum(params[0], 0.1); sigy = jnp.maximum(params[1], 0.1)
        isx = 1.0 / sigx; isy = 1.0 / sigy
        x1 = (jnp.floor(sx + 0.5) - xi - 0.5) * isx; x2 = (jnp.floor(sx + 0.5) - xi + 0.5) * isx
        y1 = (jnp.floor(sy + 0.5) - yi - 0.5) * isy; y2 = (jnp.floor(sy + 0.5) - yi + 0.5) * isy
        isq2 = 1.0 / jnp.sqrt(2.0); isq2pi = 1.0 / jnp.sqrt(2.0 * jnp.pi)
        gx1 = isq2pi * isx * jnp.exp(-0.5 * x1**2); gx2 = isq2pi * isx * jnp.exp(-0.5 * x2**2)
        gy1 = isq2pi * isy * jnp.exp(-0.5 * y1**2); gy2 = isq2pi * isy * jnp.exp(-0.5 * y2**2)
        ex = 0.5 * (jax_erf(x2 * isq2) - jax_erf(x1 * isq2)); ey = 0.5 * (jax_erf(y2 * isq2) - jax_erf(y1 * isq2))
        
        H_u = jnp.stack([hermite_pol_jnp(n, x1) for n in range(degree + 1)], axis=0)
        H2_u = jnp.stack([hermite_pol_jnp(n, x2) for n in range(degree + 1)], axis=0)
        H_v = jnp.stack([hermite_pol_jnp(n, y1) for n in range(degree + 1)], axis=0)
        H2_v = jnp.stack([hermite_pol_jnp(n, y2) for n in range(degree + 1)], axis=0)
        
        def get_basis(n, val, g1, g2, h1, h2, sig):
            """Return the n-th order 1D basis term: val (n==0, the implicit unit term) or sig*(g1*h1[n-1] - g2*h2[n-1]) (n>0), vmapped over n by the caller.

            Status: ACTIVE (production default path) -- internal helper closure of
            get_gh_basis.
            """
            return jnp.where(n == 0, val, sig * (g1 * h1[n-1] - g2 * h2[n-1]))
        
        Bx = vmap(lambda n: get_basis(n, ex, gx1, gx2, H_u, H2_u, sigx))(jnp.arange(degree + 1))
        By = vmap(lambda n: get_basis(n, ey, gy1, gy2, H_v, H2_v, sigy))(jnp.arange(degree + 1))
        return Bx, By


    @staticmethod
    def single_pix_value_np(xc, yc, xpix, ypix, params, degree):
        """NumPy reference implementation of single_pix_value_jnp (same math, scipy.special.erf instead of jax.scipy.special.erf).

        Args:
            xc, yc (float or np.ndarray): spot center, CCD pixel coordinates.
            xpix, ypix (float or np.ndarray): pixel coordinate(s) to evaluate at.
            params (np.ndarray): shape (2 + n_gh_terms[+5],), same layout as
                single_pix_value_jnp.
            degree (int): Gauss-Hermite expansion degree.

        Returns:
            Same shape as xpix/ypix: the PSF value at each pixel.

        Status: DEAD -- no callers anywhere in the active package; a NumPy
        reference implementation left alongside single_pix_value_jnp, which is
        what every real caller (pix_value_jnp, get_gh_basis's basis-decomposed
        path) actually uses.
        """
        from scipy.special import erf as scipy_erf
        sx = np.maximum(params[0], 0.1)
        sy = np.maximum(params[1], 0.1)
        isx = 1.0 / sx; isy = 1.0 / sy
        x_rel = xpix - xc; y_rel = ypix - yc
        x1 = (np.floor(xpix + 0.5) - xc - 0.5) * isx
        x2 = (np.floor(xpix + 0.5) - xc + 0.5) * isx
        y1 = (np.floor(ypix + 0.5) - yc - 0.5) * isy
        y2 = (np.floor(ypix + 0.5) - yc + 0.5) * isy
        isq2 = 1.0 / np.sqrt(2.0); isq2pi = 1.0 / np.sqrt(2.0 * np.pi)
        gx1 = isq2pi * isx * np.exp(-0.5 * x1**2); gx2 = isq2pi * isx * np.exp(-0.5 * x2**2)
        gy1 = isq2pi * isy * np.exp(-0.5 * y1**2); gy2 = isq2pi * isy * np.exp(-0.5 * y2**2)
        ex = 0.5 * (scipy_erf(x2 * isq2) - scipy_erf(x1 * isq2))
        ey = 0.5 * (scipy_erf(y2 * isq2) - scipy_erf(y1 * isq2))
        psfval = ex * ey
        param_index = 2
        for j in range(degree + 1):
            pj = ey if j == 0 else sy * (gy1 * hermite_pol_np(j-1, y1) - gy2 * hermite_pol_np(j-1, y2))
            imin = 1 if j == 0 else 0
            for i in range(imin, degree + 1):
                pi = ex if i == 0 else sx * (gx1 * hermite_pol_np(i-1, x1) - gx2 * hermite_pol_np(i-1, x2))
                psfval += params[param_index] * pj * pi
                param_index += 1
        if params.shape[0] > param_index:
            t_amp, t_core, t_xsca, t_ysca, t_inde = params[param_index:param_index+5]
            r2_core = t_core**2
            r2 = (x_rel * t_xsca)**2 + (y_rel * t_ysca)**2
            denom = r2_core + r2 + 1e-10
            tail = t_amp * r2 / denom * (denom)**(-t_inde/2.0)
            psfval += tail
        return psfval

    @staticmethod
    @partial(jit, static_argnums=(5,))
    def pix_value_jnp(xc, yc, xpix, ypix, params, degree):
        # xc, yc, params are scalars (for AD) or arrays (for bulk)
        # Handle scalar inputs correctly for jax.jacobian
        """Vectorized, jitted wrapper around single_pix_value_jnp: evaluate the PSF over a whole pixel stamp (xpix, ypix) for one or many spots at once, handling both a single scalar (xc, yc, params) triple (used under jax.jacobian/jacfwd for autodiff) and batched arrays (one spot per leading-axis element) transparently.

        Args:
            xc, yc (float or jnp.ndarray): spot center(s); scalar for a single
                spot (AD path) or shape (n_spots,) for a batch.
            xpix, ypix (jnp.ndarray): shared pixel-stamp coordinates, shape
                (n_pix,), evaluated for every spot.
            params (jnp.ndarray): shape (n_params,) for a single spot, or
                (n_spots, n_params) for a batch.
            degree (int): Gauss-Hermite expansion degree (static arg for jit).

        Returns:
            jnp.ndarray: shape (n_pix,) for a single spot, or (n_spots, n_pix) for
            a batch.

        Status: CI-TEST-ONLY -- not called by the production fit path
        (`fit_ccd_native`/`fit_bundle_task`/`PSF_Fitter.fit`, which all call
        `single_pix_value_jnp` directly via `_predict_bundle_jax`/
        `_fit_all_spots_batch`), but real: exercised by the pytest suites
        `testing/test_math_psf.py` and `testing/test_vectorization.py`,
        which CI actually runs (see `docs/python-port/code-reading-guide.md`).
        """
        is_scalar = jnp.ndim(xc) == 0
        if is_scalar:
            if jnp.ndim(xpix) == 0:
                # Single spot, single pixel -- nothing to vmap over (vmap
                # requires rank >= 1 along the mapped axis; a bare scalar
                # xpix/ypix has rank 0). Call the scalar kernel directly.
                return GaussHermitePSF.single_pix_value_jnp(xc, yc, xpix, ypix, params, degree)
            return vmap(GaussHermitePSF.single_pix_value_jnp, in_axes=(None, None, 0, 0, None, None))(xc, yc, xpix, ypix, params, degree)
        else:
            def spot_pix(xc_s, yc_s, p_s):
                """Evaluate single_pix_value_jnp for one spot's (xc_s, yc_s, p_s) over the shared pixel-stamp grid, vmapped by the caller over the batch axis.

                Status: CI-TEST-ONLY -- internal helper closure of
                pix_value_jnp's batched-input branch, see its docstring.
                """
                if jnp.ndim(xpix) == 0:
                    # Same rank-0 case as the is_scalar branch above, just
                    # for a batch of spots evaluated at a single pixel.
                    return GaussHermitePSF.single_pix_value_jnp(xc_s, yc_s, xpix, ypix, p_s, degree)
                return vmap(GaussHermitePSF.single_pix_value_jnp, in_axes=(None, None, 0, 0, None, None))(xc_s, yc_s, xpix, ypix, p_s, degree)
            return vmap(spot_pix)(xc, yc, params)

    def pix_value(self, xc, yc, xpix, ypix, params, use_jax=True):
        """Plain-NumPy-facing wrapper around pix_value_jnp: convert inputs to jnp arrays, evaluate, convert the result back to np.ndarray.

        Args:
            xc, yc (float or array-like): spot center(s).
            xpix, ypix (array-like): pixel-stamp coordinates.
            params (array-like): GH parameters, see pix_value_jnp.
            use_jax (bool): accepted but unused -- pix_value_jnp is always called
                regardless of this flag's value.

        Returns:
            float or np.ndarray: a plain Python float for fully scalar
            (xc, yc, xpix, ypix) input (matching `isinstance(..., float)`,
            as `testing/test_vectorization.py` checks), else an
            `np.ndarray` -- see `pix_value_jnp`'s docstring for shapes.

        Status: CI-TEST-ONLY, see pix_value_jnp's docstring.
        """
        result = np.array(GaussHermitePSF.pix_value_jnp(jnp.array(xc), jnp.array(yc), jnp.array(xpix), jnp.array(ypix), jnp.array(params), self.degree))
        if result.ndim == 0:
            return float(result)
        return result

class PSF_Params:
    """Per-bundle PSF-fit parameter/state container: which fibers the bundle spans, the fitted 2D polynomial models for each shape parameter, and continuum-fit state.

    Attributes:
        bundle_id (int): bundle index (0-19).
        fiber_min, fiber_max (int): inclusive fiber range of the bundle.
        param_names (list[str]): names of parameters actually being fit for
            this bundle (a subset/ordering of canonical_param_names).
        all_par_pol_xw, fit_par_pol_xw (list): per-parameter
            SparseLegendre2DPol models (all vs. actively-fit subset).
        continuum_pol: continuum-background polynomial model, or None if
            continuum fitting is off for this bundle.
        continuum_sigma_x (float): continuum model's X width parameter.

    Status: ACTIVE (production default path).
    """
    def __init__(self, bundle_id, fiber_min, fiber_max):
        """Args:
            bundle_id (int): bundle index (0-19).
            fiber_min, fiber_max (int): inclusive fiber range of the bundle.

        Status: ACTIVE (production default path).
        """
        self.bundle_id = bundle_id; self.fiber_min = fiber_min; self.fiber_max = fiber_max
        self.param_names = []; self.all_par_pol_xw = []; self.fit_par_pol_xw = []; self.continuum_pol = None; self.continuum_sigma_x = 1.0

class PSF:
    """Per-camera PSF container: the Gauss-Hermite shape model, per-bundle fitted parameters (params_of_bundles), and per-fiber trace polynomials (fiber_traces) -- the in-memory representation loaded from/written to a PSF FITS file by io.py's load_python_psf/write_python_psf.

    Attributes:
        name (str): PSF model name, always "GaussHermitePSF".
        h_size_x, h_size_y (int): PSF evaluation/fit stamp half-sizes in X/Y.
        gain (float): detector gain (currently unused downstream, kept for
            parity with the C++ PSF object).
        psf_error (float): PSF model error floor (currently unused downstream).
        fiber_min (int): lowest fiber index covered by this PSF object.
        params_of_bundles (dict[int, PSF_Params]): per-bundle fitted state,
            keyed by bundle_id.
        fiber_traces (dict[int, dict]): per-fiber trace polynomials, keyed by
            absolute fiber index; each value has 'X_vs_W'/'Y_vs_W'
            Legendre1DPol models.
        gh_psf (GaussHermitePSF): the shape-evaluation model.

    Status: ACTIVE (production default path).
    """
    def __init__(self, degree=6):
        """Args:
            degree (int): Gauss-Hermite expansion degree for gh_psf, default 6.

        Status: ACTIVE (production default path).
        """
        self.name = "GaussHermitePSF"
        self.h_size_x = 8 
        self.h_size_y = 8
        self.gain = 1.0
        self.psf_error = 0.0
        self.fiber_min = 0; self.params_of_bundles = {}; self.fiber_traces = {}
        self.gh_psf = GaussHermitePSF(degree=degree)
    
    def canonical_param_names(self):
        """Names of the parameters GaussHermitePSF.single_pix_value_jnp actually consumes, in the exact order it expects them (matches C++ GaussHermitePSF::DefaultParamNames, specex_gauss_hermite_psf.cc:394-418).

        The (i=0, j=0) GH term is intentionally excluded: it's the implicit unit-
        amplitude 0th-order term, hardcoded as ex*ey in the PSF evaluation and
        never stored as a fit parameter. FITS PSF tables carry an explicit
        'GH-0-0' row (fixed at 1.0) plus non-shape bookkeeping rows
        (BUNDLE/STATUS/CONT) that must not be fed into the parameter array.

        Returns:
            list[str]: ['GHSIGX', 'GHSIGY', 'GH-<i>-<j>', ..., 'TAILAMP',
            'TAILCORE', 'TAILXSCA', 'TAILYSCA', 'TAILINDE'].

        Status: ACTIVE (production default path).
        """
        degree = self.gh_psf.degree
        names = ['GHSIGX', 'GHSIGY']
        for j in range(degree + 1):
            for i in range(degree + 1):
                if i == 0 and j == 0:
                    continue
                names.append(f'GH-{i}-{j}')
        names += ['TAILAMP', 'TAILCORE', 'TAILXSCA', 'TAILYSCA', 'TAILINDE']
        return names

    def gh_params(self, fiber, wave):
        """Evaluate this fiber's Gauss-Hermite shape parameters at a given wavelength, by evaluating each canonical parameter's fitted 2D polynomial model.

        Args:
            fiber (int): absolute fiber index.
            wave (float): wavelength (Angstrom).

        Returns:
            np.ndarray: shape (55,) parameter vector in canonical_param_names
            order. If `fiber`'s bundle has no fitted model yet (e.g. before the
            first fit, using only the initial shift), returns a default vector
            (sigma_x=sigma_y=1.1, all else 0).

        Status: ACTIVE (production default path).
        """
        bundle_id = self.get_bundle_of_fiber(fiber)
        if bundle_id not in self.params_of_bundles:
            # Default params if no model exists (e.g. initial shift)
            p = np.zeros(55)
            p[0] = 1.1 # sigma_x
            p[1] = 1.1 # sigma_y
            return p

        params = self.params_of_bundles[bundle_id]
        # param_models lists are indexed by ABSOLUTE fiber (load_python_psf
        # builds them with `for fib in range(500)`) - do not use a
        # bundle-relative index here; that returned another bundle's shape
        # (e.g. fiber 130 got fiber 5's coefficients).
        p = []
        for name in self.canonical_param_names():
            if name in params.param_models:
                p.append(params.param_models[name][fiber].value(wave))
            else:
                # Default for GH terms not in model
                if name == 'GHSIGX' or name == 'GHSIGY': p.append(1.1)
                else: p.append(0.0)
        return np.array(p)

    def x_ccd(self, fiber, wave, tc_x=None, wdeg=3):
        """Evaluate this fiber's fitted X trace position at a given wavelength, optionally adding a joint-fit trace correction (tc_x) on top of the base per-fiber trace polynomial.

        Args:
            fiber (int): absolute fiber index.
            wave (float): wavelength (Angstrom).
            tc_x (array-like or None): trace-correction coefficients in the
                shared bundle-relative-fiber x wavelength sparse basis (see
                get_sparse_nz(1, wdeg)); None means no correction applied.
            wdeg (int): wavelength degree of the tc_x basis (only used if tc_x is
                given).

        Returns:
            float: X CCD position, or 0.0 if `fiber` has no trace loaded.

        Status: ACTIVE (production default path).
        """
        if fiber in self.fiber_traces:
            val = self.fiber_traces[fiber]['X_vs_W'].value(wave)
            if tc_x is not None:
                bundle_id = self.get_bundle_of_fiber(fiber)
                bundle = self.params_of_bundles[bundle_id]
                fmin, fmax = bundle.fiber_min, bundle.fiber_max
                wmin, wmax = self.fiber_traces[fmin]['X_vs_W'].xmin, self.fiber_traces[fmin]['X_vs_W'].xmax
                rf = 2 * (fiber - fmin) / (fmax - fmin) - 1
                rw = 2 * (wave - wmin) / (wmax - wmin) - 1
                # NumPy, not JAX: called from a plain-Python per-candidate hot
                # loop (trace warm-up in select_bundle_spots_iterative) where
                # JAX's eager-mode per-op GPU dispatch overhead (~1-2ms/op)
                # dominates runtime for what's otherwise a handful of scalar
                # flops -- same class of bottleneck as gh_params/monomials().
                from .math import legendre_pol
                xdeg = 1
                mx = [legendre_pol(i, rf) for i in range(xdeg + 1)]
                mw = [legendre_pol(j, rw) for j in range(wdeg + 1)]
                nz = []
                for j in range(wdeg + 1):
                    for i in range(xdeg + 1):
                        if i == 0: nz.append(i + j*(xdeg + 1))
                        elif i == 1 and j < 2: nz.append(i + j*(xdeg + 1))
                        elif i > 1 and j == 0: nz.append(i + j*(xdeg + 1))
                m = []
                for k in nz:
                    i, j = k % (xdeg + 1), k // (xdeg + 1)
                    m.append(mx[i] * mw[j])
                val += np.dot(m, tc_x)
            return val
        return 0.0
    def y_ccd(self, fiber, wave, tc_y=None, wdeg=3):
        """Evaluate this fiber's fitted Y trace position at a given wavelength, optionally adding a joint-fit trace correction (tc_y). Same structure as x_ccd, for the Y axis.

        Args:
            fiber (int): absolute fiber index.
            wave (float): wavelength (Angstrom).
            tc_y (array-like or None): trace-correction coefficients, see x_ccd.
            wdeg (int): wavelength degree of the tc_y basis.

        Returns:
            float: Y CCD position, or 0.0 if `fiber` has no trace loaded.

        Status: ACTIVE (production default path).
        """
        if fiber in self.fiber_traces:
            val = self.fiber_traces[fiber]['Y_vs_W'].value(wave)
            if tc_y is not None:
                bundle_id = self.get_bundle_of_fiber(fiber)
                bundle = self.params_of_bundles[bundle_id]
                fmin, fmax = bundle.fiber_min, bundle.fiber_max
                wmin, wmax = self.fiber_traces[fmin]['X_vs_W'].xmin, self.fiber_traces[fmin]['X_vs_W'].xmax
                rf = 2 * (fiber - fmin) / (fmax - fmin) - 1
                rw = 2 * (wave - wmin) / (wmax - wmin) - 1
                from .math import legendre_pol
                xdeg = 1
                mx = [legendre_pol(i, rf) for i in range(xdeg + 1)]
                mw = [legendre_pol(j, rw) for j in range(wdeg + 1)]
                nz = []
                for j in range(wdeg + 1):
                    for i in range(xdeg + 1):
                        if i == 0: nz.append(i + j*(xdeg + 1))
                        elif i == 1 and j < 2: nz.append(i + j*(xdeg + 1))
                        elif i > 1 and j == 0: nz.append(i + j*(xdeg + 1))
                m = []
                for k in nz:
                    i, j = k % (xdeg + 1), k // (xdeg + 1)
                    m.append(mx[i] * mw[j])
                val += np.dot(m, tc_y)
            return val
        return 0.0
    def get_bundle_of_fiber(self, fiber):
        """Look up which bundle a fiber belongs to, by linear scan of params_of_bundles' fiber_min/fiber_max ranges.

        Args:
            fiber (int): absolute fiber index.

        Returns:
            int: bundle_id, or -1 if no loaded bundle covers this fiber.

        Status: ACTIVE (production default path).
        """
        for bundle_id, params in self.params_of_bundles.items():
            if params.fiber_min <= fiber <= params.fiber_max: return bundle_id
        return -1
