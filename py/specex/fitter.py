import os
import numpy as np
import jax
import jax.numpy as jnp
from jax import jit, vmap, jacfwd, lax
from functools import partial
from .psf import GaussHermitePSF

# --- Helpers ---

def read_lamp_lines(filename):
    lines = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'): continue
            parts = line.split()
            try:
                wave = float(parts[1])
                name = parts[0]
                lines.append({'wave': wave, 'name': name})
            except (ValueError, IndexError): continue
    return lines

def get_bundle_spots(psf, fiber_min, fiber_max, lamp_lines, min_dist_angstrom=4.0, wave_min=None, wave_max=None):
    spots = []
    filtered_lines = [l for l in lamp_lines if (wave_min is None or l['wave'] >= wave_min) and (wave_max is None or l['wave'] <= wave_max)]
    sorted_lines = sorted(filtered_lines, key=lambda x: x['wave'])
    valid_waves = []
    for i, line in enumerate(sorted_lines):
        dist_prev = line['wave'] - sorted_lines[i-1]['wave'] if i > 0 else 1e9
        dist_next = sorted_lines[i+1]['wave'] - line['wave'] if i < len(sorted_lines)-1 else 1e9
        if min(dist_prev, dist_next) >= min_dist_angstrom: valid_waves.append(line['wave'])
    for fiber in range(fiber_min, fiber_max + 1):
        for wave in valid_waves:
            xc = psf.x_ccd(fiber, wave); yc = psf.y_ccd(fiber, wave)
            if 0 <= xc < 4114 and 0 <= yc < 4128:
                hsize_x, hsize_y = psf.h_size_x, psf.h_size_y
                spots.append({
                    'fiber': fiber, 'wave': wave, 'xc': xc, 'yc': yc, 'flux': 1000.0,
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
    pixels = np.array(pixels)
    return pixels[:, 0], pixels[:, 1]

# --- JAX Differentiable Model ---

def spot_model_fn(xc, yc, flux, gh_params, xpix, ypix, degree):
    """
    Contribution of ONE spot to pixels.
    """
    p = vmap(GaussHermitePSF.single_pix_value_jnp, in_axes=(None, None, 0, 0, None, None))(xc, yc, xpix, ypix, gh_params, degree)
    return flux * p

# --- JAX Matrix Filling ---

@partial(jit, static_argnums=(8,))
def compute_chi2_ab_flux_only_jnp(xc, yc, flux, gh_params, xpix, ypix, image_data, weight_data, degree):
    psf_vals = GaussHermitePSF.pix_value_jnp(xc, yc, xpix, ypix, gh_params, degree)
    signal = jnp.dot(psf_vals, flux)
    res = image_data - signal
    chi2 = jnp.sum(weight_data * res**2)
    B = jnp.dot(psf_vals.T, weight_data * res)
    A = jnp.dot(psf_vals.T * weight_data, psf_vals)
    return chi2, A, B

@partial(jit, static_argnums=(9,))
def compute_chi2_ab_full_chunk_jnp(xc, yc, flux, gh_params, continuum, xpix, ypix, image_data, weight_data, degree):
    """
    Chunked Jacobian accumulation to avoid OOM.
    Fits: Flux, XC, YC per spot.
    """
    Ns = xc.shape[0]
    Np = xpix.shape[0]
    chunk_size = 5000
    n_chunks = (Np + chunk_size - 1) // chunk_size
    
    Ntot = 3 * Ns + 1 # flux, xc, yc, continuum
    
    def chunk_body(i, carry):
        A, B, total_chi2 = carry
        start = i * chunk_size
        end = jnp.minimum(start + chunk_size, Np)
        
        chunk_x = xpix[start:end]
        chunk_y = ypix[start:end]
        chunk_img = image_data[start:end]
        chunk_w = weight_data[start:end]
        
        # 1. Model for this chunk
        def chunk_signal_fn(f, x, y, c):
            p_vals = GaussHermitePSF.pix_value_jnp(x, y, chunk_x, chunk_y, gh_params, degree)
            return jnp.dot(p_vals, f) + c
        
        signal = chunk_signal_fn(flux, xc, yc, continuum)
        res = chunk_img - signal
        chi2 = jnp.sum(chunk_w * res**2)
        
        # 2. Jacobian of chunk signal wrt parameters
        # For simplicity, we use jacfwd on the chunk. 
        # (end-start, Ntot)
        jac_fn = jacfwd(chunk_signal_fn, argnums=(0, 1, 2, 3))
        df, dx, dy, dc = jac_fn(flux, xc, yc, continuum)
        # dc is (Np_chunk, 1)
        H = jnp.concatenate([df, dx, dy, dc], axis=1) # (Np_chunk, Ntot)
        
        # 3. Accumulate
        A = A + jnp.dot(H.T * chunk_w, H)
        B = B + jnp.dot(H.T, chunk_w * res)
        return A, B, total_chi2 + chi2

    A = jnp.zeros((Ntot, Ntot))
    B = jnp.zeros(Ntot)
    
    # We use a loop for now (easier to debug)
    # lax.fori_loop would be better for performance
    for i in range(n_chunks):
        A, B, chi2 = chunk_body(i, (A, B, 0.0))
        
    return chi2, A, B

# --- Brent's Method ---
CGOLD = 0.3819660
ZEPS = 1e-10

def brent(func, ax, bx, cx, tol=1e-2, itmax=100):
    a = min(ax, cx); b = max(ax, cx)
    x = w = v = bx; fw = fv = fx = func(x)
    d = e = 0.0
    for i in range(itmax):
        xm = 0.5 * (a + b); tol1 = tol * abs(x) + ZEPS; tol2 = 2.0 * tol1
        if abs(x - xm) <= (tol2 - 0.5 * (b - a)): return x, fx
        if abs(e) > tol1:
            r = (x - w) * (fx - fv); q = (x - v) * (fx - fw)
            p = (x - v) * q - (x - w) * r; q = 2.0 * (q - r)
            if q > 0.0: p = -p
            q = abs(q); etemp = e; e = d
            if abs(p) >= abs(0.5 * q * etemp) or p <= q * (a - x) or p >= q * (b - x):
                e = a - x if x >= xm else b - x; d = CGOLD * e
            else:
                d = p / q; u = x + d
                if (u - a) < tol2 or (b - u) < tol2: d = np.sign(xm - x) * tol1
        else:
            e = a - x if x >= xm else b - x; d = CGOLD * e
        u = x + d if abs(d) >= tol1 else x + np.sign(d) * tol1
        fu = func(u)
        if fu <= fx:
            if u >= x: a = x
            else: b = x
            v, w, x = w, x, u; fv, fw, fx = fw, fx, fu
        else:
            if u < x: a = u
            else: b = u
            if fu <= fw or w == x: v, w = w, u; fv, fw = fw, fu
            elif fu <= fv or v == x or v == w: v = u; fv = fu
    return x, fx

# --- PSF Fitter Class ---

class PSF_Fitter:
    def __init__(self, psf):
        self.psf = psf
        self.chi2_precision = 1e-4

    def fit(self, image, weight, spots, bundle_id, fit_type='flux', max_iter=15):
        print(f"Starting Python/JAX {fit_type} fit for bundle {bundle_id}...")
        flux = np.array([s['flux'] for s in spots])
        xc = np.array([s['xc'] for s in spots])
        yc = np.array([s['yc'] for s in spots])
        continuum = 0.0
        for s in spots:
            if 'gh_params' not in s: s['gh_params'] = self.psf.all_local_params_fw(s['fiber'], s['wave'], bundle_id)
        gh_params = np.array([s['gh_params'] for s in spots])
        xpix, ypix = get_bundle_footprint(spots)
        xpix = jnp.array(xpix); ypix = jnp.array(ypix)
        image_data = jnp.array(image[xpix, ypix]); weight_data = jnp.array(weight[xpix, ypix])
        degree = self.psf.gh_psf.degree; old_chi2 = 1e30
        for i in range(max_iter):
            if fit_type == 'flux':
                chi2, A, B = compute_chi2_ab_flux_only_jnp(jnp.array(xc), jnp.array(yc), jnp.array(flux), jnp.array(gh_params), xpix, ypix, image_data, weight_data, degree)
            else:
                chi2, A, B = compute_chi2_ab_full_chunk_jnp(jnp.array(xc), jnp.array(yc), jnp.array(flux), jnp.array(gh_params), continuum, xpix, ypix, image_data, weight_data, degree)
            current_chi2 = float(chi2); print(f"Iter {i}: chi2 = {current_chi2:.4f} npix = {len(xpix)}")
            if abs(old_chi2 - current_chi2) < self.chi2_precision: break
            old_chi2 = current_chi2
            A_reg = A + 1e-6 * jnp.eye(A.shape[0])
            try: delta_P = jnp.array(jax.lax.linalg.cholesky_solve(B, jax.lax.linalg.cholesky(A_reg, upper=False), left_side=True))
            except: delta_P = np.array(jnp.linalg.solve(A_reg, B))
            delta_P = np.array(delta_P)
            def step_func(step):
                Ns = len(flux)
                if fit_type == 'flux':
                    new_flux = flux + step * delta_P[:Ns]
                    new_xc, new_yc, new_cont = xc, yc, 0.0
                else:
                    new_flux = flux + step * delta_P[:Ns]
                    new_xc = xc + step * delta_P[Ns:2*Ns]
                    new_yc = yc + step * delta_P[2*Ns:3*Ns]
                    new_cont = continuum + step * delta_P[-1]
                c, _, _ = compute_chi2_ab_flux_only_jnp(jnp.array(new_xc), jnp.array(new_yc), jnp.array(new_flux), jnp.array(gh_params), xpix, ypix, image_data, weight_data, degree)
                return float(c)
            best_step = 1.0 if fit_type == 'flux' else brent(step_func, -0.05, 0.5, 1.05)[0]
            print(f"  Best step: {best_step:.4f}")
            Ns = len(flux); flux = flux + best_step * delta_P[:Ns]
            if fit_type == 'full':
                xc = xc + best_step * delta_P[Ns:2*Ns]
                yc = yc + best_step * delta_P[2*Ns:3*Ns]
                continuum = continuum + best_step * delta_P[-1]
            for idx, s in enumerate(spots):
                s['flux'] = float(flux[idx]); s['xc'] = float(xc[idx]); s['yc'] = float(yc[idx])
        return old_chi2
