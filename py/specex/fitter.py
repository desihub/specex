import numpy as np
import jax
import jax.numpy as jnp
from jax import jit, grad, vmap
from .psf import GaussHermitePSF

@jit
def compute_chi2_ab_jnp(xc, yc, flux, gh_params, xpix, ypix, image_data, weight_data, degree):
    """
    JAX implementation of ComputeChi2AB.
    """
    # 1. Evaluate PSF: shape (Np, Ns)
    psf_vals = GaussHermitePSF.pix_value_jnp(xc, yc, xpix, ypix, gh_params, degree)
    
    # 2. Compute signal: sum_s flux[s] * PSF(p, s)
    # shape (Np,)
    signal = jnp.dot(psf_vals, flux)
    
    # 3. Residuals and Chi2
    res = image_data - signal
    chi2 = jnp.sum(weight_data * res**2)
    
    # 4. Filling A and B (Flux-only for now)
    # B = H^T * W * res
    B = jnp.dot(psf_vals.T, weight_data * res)
    
    # A = H^T * W * H
    A = jnp.dot(psf_vals.T * weight_data, psf_vals)
    
    return chi2, A, B

def compute_chi2_ab(psf, image, weight, spots, use_jax=True):
    """
    Python orchestration for ComputeChi2AB.
    """
    n_spots = len(spots)
    xc = np.array([s['xc'] for s in spots])
    yc = np.array([s['yc'] for s in spots])
    flux = np.array([s['flux'] for s in spots])
    gh_params = np.array([s['gh_params'] for s in spots])
    
    i_min = int(np.min([s['stamp_imin'] for s in spots]))
    i_max = int(np.max([s['stamp_imax'] for s in spots]))
    j_min = int(np.min([s['stamp_jmin'] for s in spots]))
    j_max = int(np.max([s['stamp_jmax'] for s in spots]))
    
    ii, jj = np.meshgrid(np.arange(i_min, i_max), np.arange(j_min, j_max), indexing='ij')
    xpix = ii.flatten()
    ypix = jj.flatten()
    
    image_data = image[xpix, ypix]
    weight_data = weight[xpix, ypix]
    
    if use_jax:
        chi2, A, B = compute_chi2_ab_jnp(jnp.array(xc), jnp.array(yc), 
                                        jnp.array(flux), jnp.array(gh_params),
                                        jnp.array(xpix), jnp.array(ypix),
                                        jnp.array(image_data), jnp.array(weight_data),
                                        psf.gh_psf.degree)
        return float(chi2), np.array(A), np.array(B), len(xpix)
    else:
        # NumPy fallback (previous implementation)
        psf_vals = psf.pix_value(xc, yc, xpix, ypix, gh_params)
        signal = np.dot(psf_vals, flux)
        res = image_data - signal
        chi2 = np.sum(weight_data * res**2)
        B = np.dot(psf_vals.T, weight_data * res)
        A = np.dot(psf_vals.T * weight_data, psf_vals)
        return chi2, A, B, len(xpix)
