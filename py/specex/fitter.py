import numpy as np
import jax
import jax.numpy as jnp
from jax import jit, grad, vmap, jacfwd
from .psf import GaussHermitePSF

def model_function(xc, yc, flux, gh_params, xpix, ypix, degree):
    """
    Differentiable model function: signal = sum_s flux[s] * PSF(p, s)
    Returns shape (Np,)
    """
    # pix_value_jnp returns (Np, Ns)
    psf_vals = GaussHermitePSF.pix_value_jnp(xc, yc, xpix, ypix, gh_params, degree)
    return jnp.dot(psf_vals, flux)

@jit
def compute_chi2_ab_flux_only(xc, yc, flux, gh_params, xpix, ypix, image_data, weight_data, degree):
    """
    Optimized flux-only accumulation.
    """
    psf_vals = GaussHermitePSF.pix_value_jnp(xc, yc, xpix, ypix, gh_params, degree)
    signal = jnp.dot(psf_vals, flux)
    res = image_data - signal
    chi2 = jnp.sum(weight_data * res**2)
    B = jnp.dot(psf_vals.T, weight_data * res)
    A = jnp.dot(psf_vals.T * weight_data, psf_vals)
    return chi2, A, B

@jit
def compute_chi2_ab_full(xc, yc, flux, gh_params, xpix, ypix, image_data, weight_data, degree):
    """
    Full accumulation using Automatic Differentiation for xc, yc, and gh_params.
    """
    # 1. Compute signal and residuals
    psf_vals = GaussHermitePSF.pix_value_jnp(xc, yc, xpix, ypix, gh_params, degree)
    signal = jnp.dot(psf_vals, flux)
    res = image_data - signal
    chi2 = jnp.sum(weight_data * res**2)
    
    # 2. Jacobians via AD
    # We need d(signal)/d(flux), d(signal)/d(xc), d(signal)/d(yc), d(signal)/d(params)
    # H_flux is just psf_vals: (Np, Ns)
    
    # For xc, yc, gh_params, we use jacfwd on a per-spot basis or vmap
    # Let's use a helper that computes Jacobian of ONE spot's contribution to all pixels
    def spot_contribution(xc_s, yc_s, flux_s, gh_params_s):
        # Result is (Np,)
        # vmap over pixels for a single spot
        f = vmap(GaussHermitePSF.single_pix_value_jnp, in_axes=(None, None, 0, 0, None, None))
        return flux_s * f(xc_s, yc_s, xpix, ypix, gh_params_s, degree)

    # Jacobian wrt (xc_s, yc_s, gh_params_s)
    # Result is (Np, 1), (Np, 1), (Np, Nparams)
    jac_spot = vmap(jacfwd(spot_contribution, argnums=(0, 1, 3)), in_axes=(0, 0, 0, 0))
    
    d_xc, d_yc, d_gh = jac_spot(xc, yc, flux, gh_params)
    # d_xc: (Ns, Np), d_yc: (Ns, Np), d_gh: (Ns, Np, Nparams)
    
    # 3. Concatenate Jacobians into a full H matrix or accumulate directly
    # For now, let's just show how we'd get B for xc
    # B_xc = sum_p w_p * res_p * d(signal_p)/d(xc_s)
    # d_xc.T is (Np, Ns)
    B_xc = jnp.dot(d_xc, weight_data * res) # (Ns,)
    B_yc = jnp.dot(d_yc, weight_data * res) # (Ns,)
    
    # B_gh = sum_p w_p * res_p * d(signal_p)/d(gh_s,k)
    # d_gh is (Ns, Np, Nparams). We want (Ns, Nparams)
    B_gh = jnp.einsum('snp,p->sn', d_gh, weight_data * res)
    
    # Full B vector
    B_flux = jnp.dot(psf_vals.T, weight_data * res)
    B_full = jnp.concatenate([B_flux, B_xc, B_yc, B_gh.flatten()])
    
    # A matrix accumulation would follow similarly using jnp.dot or einsum
    # (Simplified for this step)
    
    return chi2, B_full

def compute_chi2_ab(psf, image, weight, spots, fit_type='flux', use_jax=True):
    """
    fit_type: 'flux', 'pos', 'psf', or 'full'
    """
    n_spots = len(spots)
    xc = jnp.array([s['xc'] for s in spots])
    yc = jnp.array([s['yc'] for s in spots])
    flux = jnp.array([s['flux'] for s in spots])
    gh_params = jnp.array([s['gh_params'] for s in spots])
    
    i_min = int(np.min([s['stamp_imin'] for s in spots]))
    i_max = int(np.max([s['stamp_imax'] for s in spots]))
    j_min = int(np.min([s['stamp_jmin'] for s in spots]))
    j_max = int(np.max([s['stamp_jmax'] for s in spots]))
    
    ii, jj = np.meshgrid(np.arange(i_min, i_max), np.arange(j_min, j_max), indexing='ij')
    xpix = jnp.array(ii.flatten())
    ypix = jnp.array(jj.flatten())
    
    image_data = jnp.array(image[xpix, ypix])
    weight_data = jnp.array(weight[xpix, ypix])
    
    if fit_type == 'flux':
        chi2, A, B = compute_chi2_ab_flux_only(xc, yc, flux, gh_params, xpix, ypix, 
                                              image_data, weight_data, psf.gh_psf.degree)
        return float(chi2), np.array(A), np.array(B), len(xpix)
    else:
        # Placeholder for full AD path
        chi2, B = compute_chi2_ab_full(xc, yc, flux, gh_params, xpix, ypix, 
                                      image_data, weight_data, psf.gh_psf.degree)
        return float(chi2), None, np.array(B), len(xpix)
