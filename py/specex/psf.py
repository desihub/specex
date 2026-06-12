import numpy as np
import jax.numpy as jnp
from jax import jit, vmap, config, jacfwd
from functools import partial
from jax.scipy.special import erf as jax_erf
from .math import SparseLegendre2DPol, Legendre1DPol, hermite_pol_jnp, hermite_pol_np

config.update("jax_enable_x64", True)

class GaussHermitePSF:
    def __init__(self, degree=6):
        self.degree = degree

    @staticmethod
    def single_pix_value_jnp(xc, yc, x, y, params, degree):
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
            return jnp.where(n == 0, val, sig * (g1 * h1[n-1] - g2 * h2[n-1]))
        
        Bx = vmap(lambda n: get_basis(n, ex, gx1, gx2, H_u, H2_u, sigx))(jnp.arange(degree + 1))
        By = vmap(lambda n: get_basis(n, ey, gy1, gy2, H_v, H2_v, sigy))(jnp.arange(degree + 1))
        return Bx, By


    @staticmethod
    def single_pix_value_np(xc, yc, xpix, ypix, params, degree):
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
        is_scalar = jnp.ndim(xc) == 0
        if is_scalar:
            return vmap(GaussHermitePSF.single_pix_value_jnp, in_axes=(None, None, 0, 0, None, None))(xc, yc, xpix, ypix, params, degree)
        else:
            def spot_pix(xc_s, yc_s, p_s):
                return vmap(GaussHermitePSF.single_pix_value_jnp, in_axes=(None, None, 0, 0, None, None))(xc_s, yc_s, xpix, ypix, p_s, degree)
            return vmap(spot_pix)(xc, yc, params)

    def pix_value(self, xc, yc, xpix, ypix, params, use_jax=True):
        return np.array(GaussHermitePSF.pix_value_jnp(jnp.array(xc), jnp.array(yc), jnp.array(xpix), jnp.array(ypix), jnp.array(params), self.degree))


    def pix_value(self, xc, yc, xpix, ypix, params, use_jax=True):
        return np.array(GaussHermitePSF.pix_value_jnp(jnp.array(xc), jnp.array(yc), jnp.array(xpix), jnp.array(ypix), jnp.array(params), self.degree))

class PSF_Params:
    def __init__(self, bundle_id, fiber_min, fiber_max):
        self.bundle_id = bundle_id; self.fiber_min = fiber_min; self.fiber_max = fiber_max
        self.param_names = []; self.all_par_pol_xw = []; self.fit_par_pol_xw = []; self.continuum_pol = None; self.continuum_sigma_x = 1.0

class PSF:
    def __init__(self, degree=6):
        self.name = "GaussHermitePSF"
        self.h_size_x = 8 
        self.h_size_y = 8
        self.gain = 1.0
        self.fiber_min = 0; self.params_of_bundles = {}; self.fiber_traces = {}
        self.gh_psf = GaussHermitePSF(degree=degree)
    
    def gh_params(self, fiber, wave):
        bundle_id = self.get_bundle_of_fiber(fiber)
        if bundle_id not in self.params_of_bundles:
            # Default params if no model exists (e.g. initial shift)
            p = np.zeros(55)
            p[0] = 1.1 # sigma_x
            p[1] = 1.1 # sigma_y
            return p
            
        params = self.params_of_bundles[bundle_id]
        rel_fiber_idx = fiber - params.fiber_min
        p = []
        for name in params.param_names:
            if name in params.param_models:
                p.append(params.param_models[name][rel_fiber_idx].value(wave))
            else:
                # Default for GH terms not in model
                if name == 'GHSIGX' or name == 'GHSIGY': p.append(1.1)
                elif name == 'GH-0-0': p.append(1.0)
                else: p.append(0.0)
        return np.array(p)

    def x_ccd(self, fiber, wave):
        if fiber in self.fiber_traces: return self.fiber_traces[fiber]['X_vs_W'].value(wave)
        return 0.0
    def y_ccd(self, fiber, wave):
        if fiber in self.fiber_traces: return self.fiber_traces[fiber]['Y_vs_W'].value(wave)
        return 0.0
    def get_bundle_of_fiber(self, fiber):
        for bundle_id, params in self.params_of_bundles.items():
            if params.fiber_min <= fiber <= params.fiber_max: return bundle_id
        return -1
