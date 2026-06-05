import os
import numpy as np
import jax
import jax.numpy as jnp
from jax import jit, vmap, jacfwd, lax
from functools import partial
from .psf import GaussHermitePSF

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

def get_bundle_spots(psf, fiber_min, fiber_max, lamp_lines, image=None, weight=None, 
                     min_dist_angstrom=0.0, sn_threshold=0.0, wave_min=None, wave_max=None):
    lines = [l for l in lamp_lines if (wave_min is None or l['wave'] >= wave_min) and (wave_max is None or l['wave'] <= wave_max)]
    lines = sorted(lines, key=lambda x: x['wave'])
    selected_waves = []
    for i, line in enumerate(lines):
        dist_prev = line['wave'] - lines[i-1]['wave'] if i > 0 else 1e9
        dist_next = lines[i+1]['wave'] - line['wave'] if i < len(lines)-1 else 1e9
        if min(dist_prev, dist_next) >= min_dist_angstrom: selected_waves.append(line['wave'])
    spots = []
    for fiber in range(fiber_min, fiber_max + 1):
        for wave in selected_waves:
            xc = psf.x_ccd(fiber, wave); yc = psf.y_ccd(fiber, wave)
            if 0 <= xc < 4114 and 0 <= yc < 4128:
                hsize_x, hsize_y = psf.h_size_x, psf.h_size_y
                spot = {
                    'fiber': fiber, 'wave': wave, 'xc_init': xc, 'yc_init': yc, 'flux': 1000.0,
                    'stamp_imin': int(np.floor(xc + 0.5)) - hsize_x, 'stamp_imax': int(np.floor(xc + 0.5)) + hsize_x + 1,
                    'stamp_jmin': int(np.floor(yc + 0.5)) - hsize_y, 'stamp_jmax': int(np.floor(yc + 0.5)) + hsize_y + 1
                }
                spots.append(spot)
    return spots

def get_bundle_footprint(spots):
    pixels = set()
    for s in spots:
        for i in range(s['stamp_imin'], s['stamp_imax']):
            for j in range(s['stamp_jmin'], s['stamp_jmax']):
                if 0 <= i < 4114 and 0 <= j < 4128: pixels.add((i, j))
    pixels = sorted(list(pixels))
    if not pixels: return np.array([]), np.array([]), {}
    pixels = np.array(pixels); pix_to_idx = {(p[0], p[1]): i for i, p in enumerate(pixels)}
    return pixels[:, 0], pixels[:, 1], pix_to_idx

def get_bundle_monomials_jnp(psf, bundle_id, spots):
    from .math import legendre_pol_jnp
    bundle = psf.params_of_bundles[bundle_id]
    xdeg, wdeg = 1, 3; nz = get_sparse_nz(xdeg, wdeg)
    fmin, fmax = bundle.fiber_min, bundle.fiber_max
    wmin, wmax = psf.fiber_traces[fmin]['X_vs_W'].xmin, psf.fiber_traces[fmin]['X_vs_W'].xmax
    fiber_vals = jnp.array([s['fiber'] for s in spots]); wave_vals = jnp.array([s['wave'] for s in spots])
    rf = 2 * (fiber_vals - fmin) / (fmax - fmin) - 1; rw = 2 * (wave_vals - wmin) / (wmax - wmin) - 1
    mx = [legendre_pol_jnp(i, rf) for i in range(xdeg + 1)]; mw = [legendre_pol_jnp(j, rw) for j in range(wdeg + 1)]
    m = []
    for k in nz: i, j = k % (xdeg + 1), k // (xdeg + 1); m.append(mx[i] * mw[j])
    return jnp.stack(m, axis=1)

# --- Masked GPU Accumulation ---

@partial(jit, static_argnums=(15,))
def accumulate_bundle_gpu_jnp(flux, psf_coeffs, trace_coeffs, continuum, 
                             xc_init, yc_init, monomials, 
                             image_data, weight_data, xpix, ypix, 
                             imin_all, imax_all, jmin_all, jmax_all, degree):
    Ns = flux.shape[0]; Npoly = monomials.shape[1]
    Nparams = psf_coeffs.shape[0]; Nshared = (Nparams + 2) * Npoly
    Np = xpix.shape[0]

    gh_all = jnp.dot(monomials, psf_coeffs.T)
    dx = jnp.dot(monomials, trace_coeffs[0]); dy = jnp.dot(monomials, trace_coeffs[1])
    xc_all, yc_all = xc_init + dx, yc_init + dy

    def scan_body(carry, i):
        p_mat_sum, j_shared_sum = carry
        
        # Spot physics with Masking to match C++ stamp logic
        def spot_psf(x, y, g):
            val = GaussHermitePSF.single_pix_value_jnp(x, y, xpix, ypix, g, degree)
            # Mask out pixels outside this spot's specific stamp
            mask = (xpix >= imin_all[i]) & (xpix < imax_all[i]) & (ypix >= jmin_all[i]) & (ypix < jmax_all[i])
            return val * mask
        
        psf_s = spot_psf(xc_all[i], yc_all[i], gh_all[i])
        jx, jy, jg = jacfwd(spot_psf, argnums=(0, 1, 2))(xc_all[i], yc_all[i], gh_all[i])
        
        m = monomials[i]; f = flux[i]
        h_s_x = f * jx[:, jnp.newaxis] * m
        h_s_y = f * jy[:, jnp.newaxis] * m
        h_s_gh = (f * jg[:, :, jnp.newaxis] * m).reshape(Np, -1)
        js_s = jnp.concatenate([h_s_gh, h_s_x, h_s_y], axis=1)
        
        return (p_mat_sum.at[:, i].set(psf_s), j_shared_sum + js_s), None

    init_carry = (jnp.zeros((Np, Ns)), jnp.zeros((Np, Nshared)))
    (p_matrix, j_shared), _ = lax.scan(scan_body, init_carry, jnp.arange(Ns))

    signal = jnp.dot(p_matrix, flux) + continuum; res = image_data - signal
    chi2 = jnp.sum(weight_data * res**2)
    B_f = jnp.dot(p_matrix.T, weight_data * res)
    A_ff = jnp.dot(p_matrix.T * weight_data, p_matrix)
    B_c = jnp.sum(weight_data * res); A_cc = jnp.sum(weight_data)
    A_fc = jnp.dot(p_matrix.T, weight_data)
    B_s = jnp.dot(j_shared.T, weight_data * res); A_ss = jnp.dot(j_shared.T * weight_data, j_shared)
    A_sf = jnp.dot(j_shared.T * weight_data, p_matrix); A_sc = jnp.dot(j_shared.T, weight_data)
    Ntot = Ns + Nshared + 1; A = jnp.zeros((Ntot, Ntot)); B = jnp.zeros(Ntot)
    B = B.at[:Ns].set(B_f).at[Ns:-1].set(B_s).at[-1].set(jnp.sum(weight_data * res))
    A = A.at[:Ns, :Ns].set(A_ff).at[Ns:-1, Ns:-1].set(A_ss).at[-1, -1].set(jnp.sum(weight_data))
    A = A.at[Ns:-1, :Ns].set(A_sf).at[:Ns, Ns:-1].set(A_sf.T)
    cross_fc = jnp.dot(p_matrix.T, weight_data)
    A = A.at[:Ns, -1].set(cross_fc).at[-1, :Ns].set(cross_fc)
    cross_sc = jnp.dot(j_shared.T, weight_data)
    A = A.at[Ns:-1, -1].set(cross_sc).at[-1, Ns:-1].set(cross_sc)
    return chi2, A, B

# --- PSF Fitter Class ---

class PSF_Fitter:
    def __init__(self, psf):
        self.psf = psf
        self.chi2_precision = 1e-4

    def fit(self, image, weight, spots, bundle_id, fit_type='full', max_iter=50):
        print(f"Starting MASKED Spot-Scan fit for bundle {bundle_id}...")
        
        flux = jnp.array([s['flux'] for s in spots])
        xc_init = jnp.array([s['xc_init'] for s in spots]); yc_init = jnp.array([s['yc_init'] for s in spots])
        monomials = get_bundle_monomials_jnp(self.psf, bundle_id, spots)
        Nparams = 55; Npoly = monomials.shape[1]
        
        # Stamp boundaries for masking
        imin = jnp.array([s['stamp_imin'] for s in spots]); imax = jnp.array([s['stamp_imax'] for s in spots])
        jmin = jnp.array([s['stamp_jmin'] for s in spots]); jmax = jnp.array([s['stamp_jmax'] for s in spots])
        
        psf_coeffs = jnp.zeros((Nparams, Npoly)).at[0, 0].set(1.1).at[1, 0].set(1.1)
        trace_coeffs = jnp.zeros((2, Npoly)); continuum = 0.0
        xpix_np, ypix_np, _ = get_bundle_footprint(spots)
        xpix, ypix = jnp.array(xpix_np), jnp.array(ypix_np)
        image_data, weight_data = jnp.array(image[xpix_np, ypix_np]), jnp.array(weight[xpix_np, ypix_np])
        degree = self.psf.gh_psf.degree; old_chi2 = 1e30

        for i in range(max_iter):
            chi2, A, B = accumulate_bundle_gpu_jnp(flux, psf_coeffs, trace_coeffs, continuum, xc_init, yc_init, monomials, image_data, weight_data, xpix, ypix, imin, imax, jmin, jmax, degree)
            mode = 'flux'
            if i > 1: mode = 'trace'
            if i > 6: mode = 'full'
            print(f"Iter {i}: chi2 = {float(chi2):.4f} [Mode: {mode}]")
            Ns = len(flux)
            if mode == 'flux': idx = jnp.concatenate([jnp.arange(Ns), jnp.array([A.shape[0]-1])])
            elif mode == 'trace': idx = jnp.concatenate([jnp.arange(Ns), jnp.arange(Ns + 55*Npoly, Ns + 57*Npoly), jnp.array([A.shape[0]-1])])
            else: idx = jnp.arange(A.shape[0])
            A_sub = A[jnp.ix_(idx, idx)]; B_sub = B[idx]
            diag = jnp.diag(A_sub); S = jnp.sqrt(diag); S = jnp.where(S < 1e-12, 1.0, S)
            A_reg = (A_sub / jnp.outer(S, S)) + 1e-6 * jnp.eye(A_sub.shape[0])
            delta_scaled = jnp.linalg.solve(A_reg, B_sub / S); delta_sub = delta_scaled / S
            delta_P = jnp.zeros(A.shape[0]).at[idx].set(delta_sub)
            flux = jnp.maximum(flux + delta_P[:Ns], 0.0)
            psf_coeffs = psf_coeffs + delta_P[Ns : Ns + 55*Npoly].reshape(55, Npoly)
            trace_coeffs = trace_coeffs + delta_P[Ns + 55*Npoly : -1].reshape(2, Npoly)
            continuum = continuum + delta_P[-1]
            if mode == 'full' and jnp.abs(old_chi2 - chi2) < self.chi2_precision: break
            old_chi2 = chi2
        return float(old_chi2)
