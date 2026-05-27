import numpy as np
import jax.numpy as jnp
from jax import jit, vmap, config
from functools import partial
from scipy.special import erf as scipy_erf
from jax.scipy.special import erf as jax_erf
from .math import SparseLegendre2DPol, Legendre1DPol

config.update("jax_enable_x64", True)

def hermite_pol(degree, x):
    """
    NumPy version of hermite_pol.
    """
    if degree == 0:
        return np.ones_like(x)
    if degree == 1:
        return x
    
    h_prev2 = np.ones_like(x)
    h_prev = x
    h_curr = x
    for i in range(2, degree + 1):
        h_curr = x * h_prev - (i - 1) * h_prev2
        h_prev2 = h_prev
        h_prev = h_curr
    return h_curr

def hermite_pol_jnp(degree, x):
    """
    JAX version of hermite_pol.
    """
    if degree == 0:
        return jnp.ones_like(x)
    if degree == 1:
        return x
    
    h_prev2 = jnp.ones_like(x)
    h_prev = x
    h_curr = x
    for i in range(2, degree + 1):
        h_curr = x * h_prev - (i - 1) * h_prev2
        h_prev2 = h_prev
        h_prev = h_curr
    return h_curr

class GaussHermitePSF:
    def __init__(self, degree=6):
        self.degree = degree

    def pix_value(self, xc, yc, xpix, ypix, params):
        """
        Standard NumPy/SciPy version for CPU.
        """
        from scipy.special import erf
        
        xc = np.atleast_1d(xc)
        yc = np.atleast_1d(yc)
        xpix = np.atleast_1d(xpix)
        ypix = np.atleast_1d(ypix)
        params = np.atleast_2d(params)

        sx = np.maximum(params[:, 0], 0.1)
        sy = np.maximum(params[:, 1], 0.1)
        isx = 1.0 / sx
        isy = 1.0 / sy
        
        x_center = np.floor(xpix[:, np.newaxis] + 0.5)
        y_center = np.floor(ypix[:, np.newaxis] + 0.5)
        
        x1 = (x_center - xc - 0.5) * isx
        x2 = (x_center - xc + 0.5) * isx
        y1 = (y_center - yc - 0.5) * isy
        y2 = (y_center - yc + 0.5) * isy
        
        isq2 = 1.0 / np.sqrt(2.0)
        isq2pi = 1.0 / np.sqrt(2.0 * np.pi)
        
        gx1 = isq2pi * isx * np.exp(-0.5 * x1**2)
        gx2 = isq2pi * isx * np.exp(-0.5 * x2**2)
        gy1 = isq2pi * isy * np.exp(-0.5 * y1**2)
        gy2 = isq2pi * isy * np.exp(-0.5 * y2**2)
        
        ex = 0.5 * (erf(x2 * isq2) - erf(x1 * isq2))
        ey = 0.5 * (erf(y2 * isq2) - erf(y1 * isq2))
        
        nx = self.degree + 1
        ny = self.degree + 1
        
        psfval = ex * ey
        
        param_index = 2
        for j in range(ny):
            pj = ey if j == 0 else sy * (gy1 * hermite_pol(j-1, y1) - gy2 * hermite_pol(j-1, y2))
            imin = 1 if j == 0 else 0
            for i in range(imin, nx):
                pi = ex if i == 0 else sx * (gx1 * hermite_pol(i-1, x1) - gx2 * hermite_pol(i-1, x2))
                psfval += params[:, param_index] * pj * pi
                param_index += 1
                
        if psfval.size == 1:
            return psfval.item()
        return psfval

    @staticmethod
    @partial(jit, static_argnums=(5,))
    def pix_value_jnp(xc, yc, xpix, ypix, params, degree):
        """
        JAX-accelerated version of pix_value.
        """
        sx = jnp.maximum(params[:, 0], 0.1)
        sy = jnp.maximum(params[:, 1], 0.1)
        isx = 1.0 / sx
        isy = 1.0 / sy
        
        # Grid of (Np, Ns)
        x_center = jnp.floor(xpix[:, jnp.newaxis] + 0.5)
        y_center = jnp.floor(ypix[:, jnp.newaxis] + 0.5)
        
        x1 = (x_center - xc - 0.5) * isx
        x2 = (x_center - xc + 0.5) * isx
        y1 = (y_center - yc - 0.5) * isy
        y2 = (y_center - yc + 0.5) * isy
        
        isq2 = 1.0 / jnp.sqrt(2.0)
        isq2pi = 1.0 / jnp.sqrt(2.0 * jnp.pi)
        
        gx1 = isq2pi * isx * jnp.exp(-0.5 * x1**2)
        gx2 = isq2pi * isx * jnp.exp(-0.5 * x2**2)
        gy1 = isq2pi * isy * jnp.exp(-0.5 * y1**2)
        gy2 = isq2pi * isy * jnp.exp(-0.5 * y2**2)
        
        ex = 0.5 * (jax_erf(x2 * isq2) - jax_erf(x1 * isq2))
        ey = 0.5 * (jax_erf(y2 * isq2) - jax_erf(y1 * isq2))
        
        psfval = ex * ey
        
        param_index = 2
        for j in range(degree + 1):
            pj = ey if j == 0 else sy * (gy1 * hermite_pol_jnp(j-1, y1) - gy2 * hermite_pol_jnp(j-1, y2))
            imin = 1 if j == 0 else 0
            for i in range(imin, degree + 1):
                pi = ex if i == 0 else sx * (gx1 * hermite_pol_jnp(i-1, x1) - gx2 * hermite_pol_jnp(i-1, x2))
                psfval += params[:, param_index] * pj * pi
                param_index += 1
                
        return psfval

class PSF_Params:
    def __init__(self, bundle_id, fiber_min, fiber_max):
        self.bundle_id = bundle_id
        self.fiber_min = fiber_min
        self.fiber_max = fiber_max
        self.all_par_pol_xw = [] 
        self.fit_par_pol_xw = [] 
        self.continuum_pol = None 
        self.continuum_sigma_x = 1.0

class PSF:
    def __init__(self):
        self.name = "GaussHermitePSF"
        self.h_size_x = 12
        self.h_size_y = 12
        self.gain = 1.0
        self.readout_noise = 1.0
        self.psf_error = 0.0
        self.params_of_bundles = {} 
        self.fiber_traces = {} 
        self.gh_psf = GaussHermitePSF(degree=6)

    def x_ccd(self, fiber, wave):
        if fiber in self.fiber_traces:
            return self.fiber_traces[fiber]['X_vs_W'].value(wave)
        return 0.0

    def y_ccd(self, fiber, wave):
        if fiber in self.fiber_traces:
            return self.fiber_traces[fiber]['Y_vs_W'].value(wave)
        return 0.0

    def get_bundle_of_fiber(self, fiber):
        for bundle_id, params in self.params_of_bundles.items():
            if params.fiber_min <= fiber <= params.fiber_max:
                return bundle_id
        return -1

    def all_local_params_fw(self, fiber, wave, bundle_id=-1):
        if bundle_id == -1:
            bundle_id = self.get_bundle_of_fiber(fiber)
        
        if bundle_id not in self.params_of_bundles:
            raise ValueError(f"Bundle {bundle_id} not found")
        
        x = self.x_ccd(fiber, wave)
        params = self.params_of_bundles[bundle_id]
        
        local_params = np.zeros(len(params.all_par_pol_xw))
        for i, pol in enumerate(params.all_par_pol_xw):
            local_params[i] = pol.value(x, wave)
            
        return local_params

    def pix_value(self, xc, yc, xpix, ypix, params, use_jax=False):
        if use_jax:
            return GaussHermitePSF.pix_value_jnp(jnp.array(xc), jnp.array(yc), 
                                                jnp.array(xpix), jnp.array(ypix), 
                                                jnp.array(params), self.gh_psf.degree)
        return self.gh_psf.pix_value(xc, yc, xpix, ypix, params)
