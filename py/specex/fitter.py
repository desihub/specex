import os
import numpy as np
import jax
import jax.numpy as jnp
from jax import jit, vmap, jacfwd, lax
from functools import partial
from .psf import GaussHermitePSF
from .math import legendre_pol_jnp

# --- Helpers ---

def get_sparse_nz(xdeg, ydeg):
    nz = []
    for j in range(ydeg + 1):
        for i in range(xdeg + 1):
            if i == 0: nz.append(i + j*(xdeg + 1))
            elif i == 1 and j < 2: nz.append(i + j*(xdeg + 1))
            elif i > 1 and j == 0: nz.append(i + j*(xdeg + 1))
    return nz

def read_lamp_lines(filename):
    lines = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'): continue
            parts = line.split()
            if len(parts) < 2: continue
            try:
                wave = float(parts[1]); name = parts[0]
                lines.append({'wave': wave, 'name': name})
            except (ValueError, IndexError): continue
    return lines

def get_bundle_spots(psf, fiber_min, fiber_max, lamp_lines, min_dist_angstrom=0.0, wave_min=None, wave_max=None):
    spots = []
    filtered_lines = [l for l in lamp_lines if (wave_min is None or l['wave'] >= wave_min) and (wave_max is None or l['wave'] <= wave_max)]
    for fiber in range(fiber_min, fiber_max + 1):
        for line in filtered_lines:
            wave = line['wave']
            xc = psf.x_ccd(fiber, wave); yc = psf.y_ccd(fiber, wave)
            if 0 <= xc < 4114 and 0 <= yc < 4128:
                hsize_x, hsize_y = psf.h_size_x, psf.h_size_y
                spots.append({
                    'fiber': fiber, 'wave': wave, 'xc_init': xc, 'yc_init': yc, 'flux': 1000.0,
                    'stamp_imin': int(np.floor(xc + 0.5)) - hsize_x, 'stamp_imax': int(np.floor(xc + 0.5)) + hsize_x + 1,
                    'stamp_jmin': int(np.floor(yc + 0.5)) - hsize_y, 'stamp_jmax': int(np.floor(yc + 0.5)) + hsize_y + 1
                })
    return spots

def get_bundle_footprint(spots):
    pixels = set()
    for s in spots:
        for i in range(s['stamp_imin'], s['stamp_imax']):
            for j in range(s['stamp_jmin'], s['stamp_jmax']):
                if 0 <= i < 4114 and 0 <= j < 4128: pixels.add((i, j))
    pixels = sorted(list(pixels))
    if not pixels: return np.array([]), np.array([])
    pixels = np.array(pixels); return pixels[:, 0], pixels[:, 1]

def get_bundle_monomials_jnp(psf, bundle_id, spots):
    bundle = psf.params_of_bundles[bundle_id]
    xdeg, wdeg = 1, 3 
    nz = get_sparse_nz(xdeg, wdeg)
    fmin, fmax = bundle.fiber_min, bundle.fiber_max
    wmin, wmax = psf.fiber_traces[fmin]['X_vs_W'].xmin, psf.fiber_traces[fmin]['X_vs_W'].xmax
    fiber_vals = jnp.array([s['fiber'] for s in spots])
    wave_vals = jnp.array([s['wave'] for s in spots])
    rf = 2 * (fiber_vals - fmin) / (fmax - fmin) - 1
    rw = 2 * (wave_vals - wmin) / (wmax - wmin) - 1
    mx = [legendre_pol_jnp(i, rf) for i in range(xdeg + 1)]
    mw = [legendre_pol_jnp(j, rw) for j in range(wdeg + 1)]
    m = []
    for k in nz:
        i, j = k % (xdeg + 1), k // (xdeg + 1)
        m.append(mx[i] * mw[j])
    return jnp.stack(m, axis=1)

# --- JAX Local Spot Jacobian ---

@partial(jit, static_argnums=(5,))
def compute_local_jac_batch(xc, yc, flux, gh_params, xpix, ypix, degree):
    def spot_fn(x, y, f, g, xp, yp):
        return f * GaussHermitePSF.single_pix_value_jnp(x, y, xp, yp, g, degree)
    v_pix = vmap(spot_fn, in_axes=(None, None, None, None, 0, 0))
    jac_fn = vmap(vmap(jacfwd(spot_fn, argnums=(0, 1, 2, 3)), in_axes=(None, None, None, None, 0, 0)), in_axes=(0, 0, 0, 0, None, None))
    df, dx, dy, dg = jac_fn(xc, yc, flux, gh_params, xpix, ypix)
    return df, dx, dy, dg

# --- PSF Fitter Class ---

class PSF_Fitter:
    def __init__(self, psf):
        self.psf = psf
        self.chi2_precision = 1e-4

    def fit(self, image, weight, spots, bundle_id, fit_type='flux', max_iter=15):
        print(f"Starting Python/JAX {fit_type} fit for bundle {bundle_id}...")
        
        flux = np.array([s['flux'] for s in spots])
        xc_init = jnp.array([s['xc_init'] for s in spots])
        yc_init = jnp.array([s['yc_init'] for s in spots])
        for s in spots:
            if 'gh_params' not in s: s['gh_params'] = self.psf.all_local_params_fw(s['fiber'], s['wave'], bundle_id)
        gh_params = np.array([s['gh_params'] for s in spots])
        
        xpix, ypix = get_bundle_footprint(spots)
        image_data, weight_data = np.array(image[xpix, ypix]), np.array(weight[xpix, ypix])
        degree = self.psf.gh_psf.degree; old_chi2 = 1e30

        monomials = np.array(get_bundle_monomials_jnp(self.psf, bundle_id, spots))
        Npoly = monomials.shape[1]
        
        # Internal params
        psf_coeffs = np.zeros((55, Npoly))
        psf_coeffs[0, 0] = 1.1; psf_coeffs[1, 0] = 1.1
        trace_coeffs = np.zeros((2, Npoly))
        continuum = 0.0

        for i in range(max_iter):
            # 1. Compute PSF Matrix and Signal
            gh_p_all = np.dot(monomials, psf_coeffs.T)
            dx_all = np.dot(monomials, trace_coeffs[0]); dy_all = np.dot(monomials, trace_coeffs[1])
            xc, yc = np.array(xc_init) + dx_all, np.array(yc_init) + dy_all
            
            psf_vals = np.array(GaussHermitePSF.pix_value_jnp(jnp.array(xc), jnp.array(yc), jnp.array(xpix), jnp.array(ypix), jnp.array(gh_p_all), degree))
            signal = np.dot(psf_vals, flux) + continuum
            res = image_data - signal
            chi2 = float(np.sum(weight_data * res**2))
            print(f"Iter {i}: chi2 = {chi2:.4f} npix = {len(xpix)}")
            if abs(old_chi2 - chi2) < self.chi2_precision: break
            old_chi2 = chi2

            # 2. Accumulate A and B
            Ntot = len(flux) + 55*Npoly + 2*Npoly + 1
            A = np.zeros((Ntot, Ntot))
            B = np.zeros(Ntot)
            
            # Linear Parts
            B[:len(flux)] = np.dot(psf_vals.T, weight_data * res)
            A[:len(flux), :len(flux)] = np.dot(psf_vals.T * weight_data, psf_vals)
            B[-1] = np.sum(weight_data * res)
            A[-1, -1] = np.sum(weight_data)
            
            # Cross terms Flux-Continuum
            A[:len(flux), -1] = np.sum(psf_vals.T * weight_data, axis=1)
            A[-1, :len(flux)] = A[:len(flux), -1]

            if fit_type == 'full':
                # Non-linear Parts: Project Spot Jacobians to Legendre space
                spot_chunk_size = 50
                for s_start in range(0, len(flux), spot_chunk_size):
                    s_end = min(s_start + spot_chunk_size, len(flux))
                    ldf, ldx, ldy, ldg = compute_local_jac_batch(
                        jnp.array(xc[s_start:s_end]), jnp.array(yc[s_start:s_end]), 
                        jnp.array(flux[s_start:s_end]), jnp.array(gh_p_all[s_start:s_end]), 
                        jnp.array(xpix), jnp.array(ypix), degree
                    )
                    m_chunk = monomials[s_start:s_end]
                    # Project and accumulate A and B for PSF and Trace coeffs...
                    # (Implementation of the full A projection matrix indexing...)
                    pass

            # Final Solve for this Turn (Flux+Continuum only for verification)
            idx_solve = np.concatenate([np.arange(len(flux)), [-1]])
            A_sub = A[np.ix_(idx_solve, idx_solve)]
            B_sub = B[idx_solve]
            delta_P = np.linalg.solve(A_sub + 1e-6*np.eye(len(idx_solve)), B_sub)
            
            flux += delta_P[:-1]
            continuum += delta_P[-1]
            for idx, s in enumerate(spots): s['flux'] = float(flux[idx])
            
        print(f"Final continuum: {continuum:.4f}")
        return old_chi2
