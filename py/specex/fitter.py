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

# --- JAX Local Spot Derivatives ---

@partial(jit, static_argnums=(4,))
def compute_spot_stamp_jac_full(xc, yc, flux, gh_params, degree, sx, sy):
    def spot_fn(f, x, y, g): return f * GaussHermitePSF.single_pix_value_jnp(x, y, sx, sy, g, degree)
    jf, jx, jy, jg = jacfwd(spot_fn, argnums=(0, 1, 2, 3))(flux, xc, yc, gh_params)
    return jf, jx, jy, jg

# --- Evaluation Helper ---

def compute_bundle_chi2(flux, gh_params_all, dx_all, dy_all, continuum, xc_init, yc_init, image_data, weight_data, xpix, ypix, degree):
    xc, yc = np.array(xc_init) + dx_all, np.array(yc_init) + dy_all
    psf_vals = np.array(GaussHermitePSF.pix_value_jnp(jnp.array(xc), jnp.array(yc), jnp.array(xpix), jnp.array(ypix), jnp.array(gh_params_all), degree))
    signal = np.dot(psf_vals, flux) + continuum; res = image_data - signal
    res = np.nan_to_num(res, nan=0.0, posinf=0.0, neginf=0.0)
    return np.sum(weight_data * res**2), res, psf_vals

# --- Matrix Accumulation ---

def accumulate_bundle_ab(flux, gh_params_all, dx_all, dy_all, continuum, xc_init, yc_init, monomials, image_data, weight_data, xpix, ypix, degree, spots, pix_to_idx, res, psf_vals, mode):
    Ns = flux.shape[0]; Npoly = monomials.shape[1]; Nparams = gh_params_all.shape[1]; Nshared = (Nparams + 2) * Npoly; Ntot = Ns + Nshared + 1
    A = np.zeros((Ntot, Ntot)); B = np.zeros(Ntot)
    A[:Ns, :Ns] = np.dot(psf_vals.T * weight_data, psf_vals); B[:Ns] = np.dot(psf_vals.T, weight_data * res)
    B[-1] = np.sum(weight_data * res); A[-1, -1] = np.sum(weight_data)
    A[:Ns, -1] = np.sum(psf_vals.T * weight_data, axis=1); A[-1, :Ns] = A[:Ns, -1]
    if mode == 'flux': return A, B
    xc_all, yc_all = np.array(xc_init) + dx_all, np.array(yc_init) + dy_all
    H_shared = np.zeros((len(xpix), Nshared))
    for s in range(Ns):
        spot = spots[s]
        sx, sy, indices = [], [], []
        for i in range(spot['stamp_imin'], spot['stamp_imax']):
            for j in range(spot['stamp_jmin'], spot['stamp_jmax']):
                if (i, j) in pix_to_idx: sx.append(i); sy.append(j); indices.append(pix_to_idx[(i, j)])
        if not sx: continue
        sx, sy, indices = np.array(sx), np.array(sy), np.array(indices)
        jf, jx, jy, jg = compute_spot_stamp_jac_full(xc_all[s], yc_all[s], flux[s], gh_params_all[s], degree, jnp.array(sx), jnp.array(sy))
        jf, jx, jy, jg = np.array(jf), np.array(jx), np.array(jy), np.array(jg)
        m = monomials[s]
        h_s_x = jx[:, np.newaxis] * m; h_s_y = jy[:, np.newaxis] * m
        h_s_gh = (jg[:, :, np.newaxis] * m).reshape(len(sx), -1)
        h_s_shared = np.concatenate([h_s_gh, h_s_x, h_s_y], axis=1)
        H_shared[indices, :] += h_s_shared
    g_start, g_end = Ns, Ntot - 1
    B[g_start:g_end] = np.dot(H_shared.T, weight_data * res); A[g_start:g_end, g_start:g_end] = np.dot(H_shared.T * weight_data, H_shared)
    A[:Ns, g_start:g_end] = np.dot(psf_vals.T * weight_data, H_shared); A[g_start:g_end, :Ns] = A[:Ns, g_start:g_end].T
    A[g_start:g_end, -1] = np.dot(H_shared.T, weight_data); A[-1, g_start:g_end] = A[g_start:g_end, -1]
    return A, B

# --- Brent's Method ---
CGOLD = 0.3819660; ZEPS = 1e-10
def brent(func, ax, bx, cx, tol=1e-2, itmax=100):
    a = min(ax, cx); b = max(ax, cx); x = w = v = bx; fw = fv = fx = func(x); d = e = 0.0
    for i in range(itmax):
        xm = 0.5 * (a + b); tol1 = tol * abs(x) + ZEPS; tol2 = 2.0 * tol1
        if abs(x - xm) <= (tol2 - 0.5 * (b - a)): return x, fx
        if abs(e) > tol1:
            r = (x-w)*(fx-fv); q = (x-v)*(fx-fw); p = (x-v)*q-(x-w)*r; q = 2.0*(q-r)
            if q > 0.0: p = -p
            q = abs(q); etemp = e; e = d
            if abs(p) >= abs(0.5*q*etemp) or p <= q*(a-x) or p >= q*(b-x): e = a-x if x>=xm else b-x; d = CGOLD*e
            else: d = p/q; u = x+d;
            if (u-a)<tol2 or (b-u)<tol2: d = np.sign(xm-x)*tol1
        else: e = a-x if x>=xm else b-x; d = CGOLD*e
        u = x+d if abs(d)>=tol1 else x+np.sign(d)*tol1; fu = func(u)
        if fu <= fx:
            if u >= x: a = x
            else: b = x
            v,w,x = w,x,u; fv,fw,fx = fw,fx,fu
        else:
            if u < x: a = u
            else: b = u
            if fu <= fw or w == x: v,w = w,u; fv,fw = fw,fu
            elif fu <= fv or v == x or v == w: v = u; fv = fu
    return x, fx

# --- PSF Fitter Class ---

class PSF_Fitter:
    def __init__(self, psf):
        self.psf = psf
        self.chi2_precision = 1e-4

    def fit(self, image, weight, spots, bundle_id, fit_type='full', max_iter=20):
        print(f"Starting Python/JAX staged fit for bundle {bundle_id}...")
        flux = np.array([s['flux'] for s in spots])
        xc_init = np.array([s['xc_init'] for s in spots]); yc_init = np.array([s['yc_init'] for s in spots])
        monomials = np.array(get_bundle_monomials_jnp(self.psf, bundle_id, spots))
        Npoly = monomials.shape[1]
        gh_params_all = np.zeros((len(spots), 55)); gh_params_all[:, 0] = 1.1; gh_params_all[:, 1] = 1.1
        dx_all = np.zeros(len(spots)); dy_all = np.zeros(len(spots))
        psf_coeffs = np.zeros((55, Npoly)); psf_coeffs[0, 0] = 1.1; psf_coeffs[1, 0] = 1.1; trace_coeffs = np.zeros((2, Npoly))
        continuum = 0.0
        xpix, ypix, pix_to_idx = get_bundle_footprint(spots)
        image_data, weight_data = image[xpix, ypix], weight[xpix, ypix]
        degree = self.psf.gh_psf.degree; old_chi2 = 1e30

        for i in range(max_iter):
            if i > 1:
                gh_params_all = np.dot(monomials, psf_coeffs.T)
                dx_all = np.dot(monomials, trace_coeffs[0]); dy_all = np.dot(monomials, trace_coeffs[1])
            chi2, res, psf_vals = compute_bundle_chi2(flux, gh_params_all, dx_all, dy_all, continuum, xc_init, yc_init, image_data, weight_data, xpix, ypix, degree)
            
            mode = 'flux'
            if i > 1: mode = 'trace'
            if i > 6: mode = 'full'
            print(f"Iter {i}: chi2 = {chi2:.4f} [Mode: {mode}]")
            
            A, B = accumulate_bundle_ab(flux, gh_params_all, dx_all, dy_all, continuum, xc_init, yc_init, monomials, image_data, weight_data, xpix, ypix, degree, spots, pix_to_idx, res, psf_vals, mode)
            
            Ns = len(flux)
            idx_flux = np.arange(Ns); idx_psf = np.arange(Ns, Ns + 55*Npoly); idx_trace = np.arange(Ns + 55*Npoly, Ns + 57*Npoly); idx_cont = np.array([A.shape[0]-1])
            if mode == 'flux': idx_solve = np.concatenate([idx_flux, idx_cont])
            elif mode == 'trace': idx_solve = np.concatenate([idx_flux, idx_trace, idx_cont])
            else: idx_solve = np.arange(A.shape[0])
            
            A_sub = A[np.ix_(idx_solve, idx_solve)]; B_sub = B[idx_solve]
            diag = np.diag(A_sub); S = np.sqrt(diag); S[S < 1e-12] = 1.0
            A_scaled = A_sub / np.outer(S, S); B_scaled = B_sub / S
            A_reg = A_scaled + 1e-6 * np.eye(A_scaled.shape[0])
            try: delta_scaled = np.linalg.solve(A_reg, B_scaled)
            except: delta_scaled = np.linalg.lstsq(A_reg, B_scaled, rcond=1e-8)[0]
            delta_sub = delta_scaled / S
            delta_P = np.zeros(A.shape[0]); delta_P[idx_solve] = delta_sub
            
            def line_search_fn(step):
                t_flux = np.maximum(flux + step * delta_P[:Ns], 0.0); t_cont = continuum + step * delta_P[-1]
                if i <= 1: t_gh, t_dx, t_dy = gh_params_all, dx_all, dy_all; t_pc, t_tc = psf_coeffs, trace_coeffs
                else:
                    t_pc = psf_coeffs + step * delta_P[Ns : Ns + 55*Npoly].reshape(55, Npoly)
                    t_tc = trace_coeffs + step * delta_P[Ns + 55*Npoly : -1].reshape(2, Npoly)
                    t_gh = np.dot(monomials, t_pc.T); t_dx = np.dot(monomials, t_tc[0]); t_dy = np.dot(monomials, t_tc[1])
                tc2, _, _ = compute_bundle_chi2(t_flux, t_gh, t_dx, t_dy, t_cont, xc_init, yc_init, image_data, weight_data, xpix, ypix, degree)
                return tc2

            best_step, best_chi2 = brent(line_search_fn, -0.05, 1.0, 1.1, tol=1e-2)
            print(f"  Step={best_step:.4f}, Chi2_new={best_chi2:.4f}")
            
            # Apply best step
            flux = np.maximum(flux + best_step * delta_P[:Ns], 0.0)
            continuum += best_step * delta_P[-1]
            if i > 1:
                psf_coeffs += best_step * delta_P[Ns : Ns + 55*Npoly].reshape(55, Npoly)
                trace_coeffs += best_step * delta_P[Ns + 55*Npoly : -1].reshape(2, Npoly)
            
            if i == 1:
                print("  Initializing 2D Legendre model from local state...")
                for p in range(55): psf_coeffs[p, :] = np.linalg.lstsq(monomials, gh_params_all[:, p], rcond=None)[0]
                trace_coeffs[0, :] = np.linalg.lstsq(monomials, dx_all, rcond=None)[0]
                trace_coeffs[1, :] = np.linalg.lstsq(monomials, dy_all, rcond=None)[0]

            if mode == 'full' and abs(old_chi2 - best_chi2) < self.chi2_precision: 
                print("  Stagnation reached in full mode.")
                break
            old_chi2 = best_chi2
        return old_chi2
