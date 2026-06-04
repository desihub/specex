import numpy as np
import jax.numpy as jnp
from jax import jit, vmap, config, jacfwd
from functools import partial
from jax.scipy.special import erf as jax_erf
from .math import SparseLegendre2DPol, Legendre1DPol

config.update("jax_enable_x64", True)

def hermite_pol_jnp(degree, x):
    if degree == 0: return jnp.ones_like(x)
    if degree == 1: return x
    h_prev2 = jnp.ones_like(x)
    h_prev = x
    h_curr = x
    for i in range(2, degree + 1):
        h_curr = x * h_prev - (i - 1) * h_prev2
        h_prev2 = h_prev; h_prev = h_curr
    return h_curr

class GaussHermitePSF:
    def __init__(self, degree=6):
        self.degree = degree

    @staticmethod
    def single_pix_value_jnp(xc, yc, xpix, ypix, params, degree):
        sx = jnp.maximum(params[0], 0.1)
        sy = jnp.maximum(params[1], 0.1)
        isx = 1.0 / sx; isy = 1.0 / sy
        x_rel = xpix - xc; y_rel = ypix - yc
        x1 = (jnp.floor(xpix + 0.5) - xc - 0.5) * isx
        x2 = (jnp.floor(xpix + 0.5) - xc + 0.5) * isx
        y1 = (jnp.floor(ypix + 0.5) - yc - 0.5) * isy
        y2 = (jnp.floor(ypix + 0.5) - yc + 0.5) * isy
        isq2 = 1.0 / jnp.sqrt(2.0); isq2pi = 1.0 / jnp.sqrt(2.0 * jnp.pi)
        gx1 = isq2pi * isx * jnp.exp(-0.5 * x1**2); gx2 = isq2pi * isx * jnp.exp(-0.5 * x2**2)
        gy1 = isq2pi * isy * jnp.exp(-0.5 * y1**2); gy2 = isq2pi * isy * jnp.exp(-0.5 * y2**2)
        ex = 0.5 * (jax_erf(x2 * isq2) - jax_erf(x1 * isq2))
        ey = 0.5 * (jax_erf(y2 * isq2) - jax_erf(y1 * isq2))
        psfval = ex * ey
        param_index = 2
        for j in range(degree + 1):
            pj = ey if j == 0 else sy * (gy1 * hermite_pol_jnp(j-1, y1) - gy2 * hermite_pol_jnp(j-1, y2))
            imin = 1 if j == 0 else 0
            for i in range(imin, degree + 1):
                pi = ex if i == 0 else sx * (gx1 * hermite_pol_jnp(i-1, x1) - gx2 * hermite_pol_jnp(i-1, x2))
                psfval += params[param_index] * pj * pi
                param_index += 1
        if params.shape[0] > param_index:
            t_amp = params[param_index]
            t_core = params[param_index+1]
            t_xsca = params[param_index+2]
            t_ysca = params[param_index+3]
            t_inde = params[param_index+4]
            
            # C++ uses multiplication for scaling: r2 = square(dx*r_tail_x_scale)+square(dy*r_tail_y_scale)
            # and formula: r2/(r2_tail_core_size+r2)*pow(r2_tail_core_size+r2,-r_tail_power_law_index/2.)
            r2_core = t_core**2
            r2 = (x_rel * t_xsca)**2 + (y_rel * t_ysca)**2
            # Add epsilon to prevent division by zero at r=0
            denom = r2_core + r2 + 1e-10
            tail = t_amp * r2 / denom * (denom)**(-t_inde/2.0)
            psfval += tail
        return psfval

    @staticmethod
    @partial(jit, static_argnums=(5,))
    def pix_value_jnp(xc, yc, xpix, ypix, params, degree):
        vmap_spots = vmap(GaussHermitePSF.single_pix_value_jnp, in_axes=(0, 0, None, None, 0, None))
        vmap_all = vmap(vmap_spots, in_axes=(None, None, 0, 0, None, None))
        return vmap_all(xc, yc, xpix, ypix, params, degree)

    def pix_value(self, xc, yc, xpix, ypix, params, use_jax=True):
        return np.array(GaussHermitePSF.pix_value_jnp(jnp.array(xc), jnp.array(yc), jnp.array(xpix), jnp.array(ypix), jnp.array(params), self.degree))

class PSF_Params:
    def __init__(self, bundle_id, fiber_min, fiber_max):
        self.bundle_id = bundle_id; self.fiber_min = fiber_min; self.fiber_max = fiber_max
        self.param_names = []; self.all_par_pol_xw = []; self.fit_par_pol_xw = []; self.continuum_pol = None; self.continuum_sigma_x = 1.0

class PSF:
    def __init__(self, degree=6):
        self.name = "GaussHermitePSF"; self.h_size_x = 12; self.h_size_y = 12
        self.gain = 1.0; self.readout_noise = 1.0; self.psf_error = 0.0
        self.fiber_min = 0; self.params_of_bundles = {}; self.fiber_traces = {}
        self.gh_psf = GaussHermitePSF(degree=degree)

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
    def all_local_params_fw(self, fiber, wave, bundle_id=-1):
        if bundle_id == -1: bundle_id = self.get_bundle_of_fiber(fiber)
        params = self.params_of_bundles[bundle_id]
        # rel_fiber_idx is relative to the bundle's first fiber
        rel_fiber_idx = fiber - params.fiber_min
        local_params = np.zeros(len(params.param_names))
        for i, name in enumerate(params.param_names):
            local_params[i] = params.param_models[name][rel_fiber_idx].value(wave)
        return local_params
