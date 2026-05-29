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
    """
    Matches specex::SparseLegendre2DPol::Fill(true)
    """
    nz = []
    for j in range(ydeg + 1):
        for i in range(xdeg + 1):
            if i == 0: nz.append(i + j*(xdeg + 1))
            elif i == 1 and j < 2: nz.append(i + j*(xdeg + 1))
            elif i > 1 and j == 0: nz.append(i + j*(xdeg + 1))
    return nz

def get_bundle_monomials_jnp(psf, bundle_id, spots):
    """
    Pre-computes the 2D Legendre monomials for all spots in a bundle.
    """
    bundle = psf.params_of_bundles[bundle_id]
    # For now assume all params in bundle use same degrees
    # DESI baseline: xdeg=1, wdeg=3
    xdeg, wdeg = 1, 3 
    nz = get_sparse_nz(xdeg, wdeg)
    
    # Bundle limits for reduction
    # In C++, it's often the bundle fiber range and trace wave range
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
        i = k % (xdeg + 1)
        j = k // (xdeg + 1)
        m.append(mx[i] * mw[j])
    return jnp.stack(m, axis=1) # (Ns, Npoly)

# --- JAX Matrix Filling (Legendre Coeffs) ---

@partial(jit, static_argnums=(9, 10, 11))
def compute_chi2_ab_bundle_coeffs_jnp(flux, psf_coeffs, trace_coeffs, continuum, 
                                     xc_init, yc_init, monomials, 
                                     image_data, weight_data, 
                                     xpix, ypix, degree):
    """
    Fits everything simultaneously using AD on Legendre Coeffs.
    flux: (Ns,)
    psf_coeffs: (Nparams, Npoly)
    trace_coeffs: (2, Npoly)
    monomials: (Ns, Npoly)
    """
    Ns = flux.shape[0]
    Nparams, Npoly = psf_coeffs.shape
    Np = xpix.shape[0]
    
    def model_fn(f, pc, tc, cont):
        # 1. Local params from Legendre
        gh_params = jnp.dot(monomials, pc.T) # (Ns, Nparams)
        dx = jnp.dot(monomials, tc[0])
        dy = jnp.dot(monomials, tc[1])
        xc, yc = xc_init + dx, yc_init + dy
        
        # 2. PSF evaluation
        psf_vals = GaussHermitePSF.pix_value_jnp(xc, yc, xpix, ypix, gh_params, degree)
        return jnp.dot(psf_vals, f) + cont

    signal = model_fn(flux, psf_coeffs, trace_coeffs, continuum)
    res = image_data - signal
    chi2 = jnp.sum(weight_data * res**2)
    
    # Jacobian wrt all coefficients
    # Ntot = Ns (flux) + Nparams*Npoly (psf) + 2*Npoly (trace) + 1 (cont)
    jac_fn = jacfwd(model_fn, argnums=(0, 1, 2, 3))
    df, dpc, dtc, dc = jac_fn(flux, psf_coeffs, trace_coeffs, continuum)
    
    # Reshape and flatten
    H = jnp.concatenate([
        df, 
        dpc.reshape(Np, -1),
        dtc.reshape(Np, -1),
        dc
    ], axis=1)
    
    B = jnp.dot(H.T, weight_data * res)
    A = jnp.dot(H.T * weight_data, H)
    
    return chi2, A, B

# --- Rest of Helpers ---

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
    pixels = np.array(pixels)
    return pixels[:, 0], pixels[:, 1]

# --- Brent's Method (same) ---
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
            v,w,x = w,x,u; fv,fw,fx = fw,fw,fu
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

    def fit(self, image, weight, spots, bundle_id, fit_type='full', max_iter=15):
        """
        Refactored to optimize Legendre coefficients.
        """
        print(f"Starting Python/JAX Legendre fit for bundle {bundle_id}...")
        
        # 1. Setup
        flux = np.array([s['flux'] for s in spots])
        xc_init = jnp.array([s['xc_init'] for s in spots])
        yc_init = jnp.array([s['yc_init'] for s in spots])
        
        # Pre-compute monomials: (Ns, Npoly)
        monomials = get_bundle_monomials_jnp(self.psf, bundle_id, spots)
        Npoly = monomials.shape[1]
        
        # Initialize Legendre coefficients from bundle model
        bundle = self.psf.params_of_bundles[bundle_id]
        # (For this validation, I'll extract the 1D coeffs per fiber as a proxy if 2D is missing)
        # In a real fit, we start from a single 2D model.
        # Placeholder: start from constant defaults
        psf_coeffs = jnp.zeros((55, Npoly))
        # Initial sigma guess (coeffs[0] is the constant term in Legendre)
        psf_coeffs = psf_coeffs.at[0, 0].set(1.1) # GHSIGX
        psf_coeffs = psf_coeffs.at[1, 0].set(1.1) # GHSIGY
        
        trace_coeffs = jnp.zeros((2, Npoly))
        continuum = 0.0
        
        xpix, ypix = get_bundle_footprint(spots)
        xpix, ypix = jnp.array(xpix), jnp.array(ypix)
        image_data, weight_data = jnp.array(image[xpix, ypix]), jnp.array(weight[xpix, ypix])
        
        degree = self.psf.gh_psf.degree
        old_chi2 = 1e30
        
        for i in range(max_iter):
            # Evaluate model
            chi2, A, B = compute_chi2_ab_bundle_jnp(
                jnp.array(flux), psf_coeffs, trace_coeffs, continuum,
                monomials, xc_init, yc_init, image_data, weight_data, xpix, ypix, degree
            )
            
            print(f"Iter {i}: chi2 = {float(chi2):.4f}")
            if abs(old_chi2 - chi2) < self.chi2_precision: break
            old_chi2 = chi2
            
            # Solve updates
            A_reg = A + 1e-4 * jnp.eye(A.shape[0])
            try: delta_P = jnp.array(jax.lax.linalg.cholesky_solve(B, jax.lax.linalg.cholesky(A_reg, upper=False), left_side=True))
            except: delta_P = np.array(jnp.linalg.solve(A_reg, B))
            
            # Apply updates (Flux + PSF + Trace + Cont)
            Ns = len(flux)
            flux = flux + np.array(delta_P[:Ns])
            psf_coeffs = psf_coeffs + delta_P[Ns : Ns + 55*Npoly].reshape(55, Npoly)
            trace_coeffs = trace_coeffs + delta_P[Ns + 55*Npoly : Ns + 55*Npoly + 2*Npoly].reshape(2, Npoly)
            continuum = continuum + float(delta_P[-1])
            
        return old_chi2

@partial(jit, static_argnums=(10, 11))
def compute_chi2_ab_bundle_jnp(flux, psf_coeffs, trace_coeffs, continuum, 
                               monomials, xc_init, yc_init, 
                               image_data, weight_data, xpix, ypix, degree):
    """
    Actual implementation of bundle-wide accumulation.
    """
    def model_fn(f, pc, tc, cont):
        gh_params = jnp.dot(monomials, pc.T)
        dx = jnp.dot(monomials, tc[0]); dy = jnp.dot(monomials, tc[1])
        xc, yc = xc_init + dx, yc_init + dy
        psf_vals = GaussHermitePSF.pix_value_jnp(xc, yc, xpix, ypix, gh_params, degree)
        return jnp.dot(psf_vals, f) + cont

    signal = model_fn(flux, psf_coeffs, trace_coeffs, continuum)
    res = image_data - signal
    chi2 = jnp.sum(weight_data * res**2)
    
    jac_fn = jacfwd(model_fn, argnums=(0, 1, 2, 3))
    df, dpc, dtc, dc = jac_fn(flux, psf_coeffs, trace_coeffs, continuum)
    
    H = jnp.concatenate([
        df, 
        dpc.reshape(xpix.shape[0], -1),
        dtc.reshape(xpix.shape[0], -1),
        dc
    ], axis=1)
    
    B = jnp.dot(H.T, weight_data * res)
    A = jnp.dot(H.T * weight_data, H)
    return chi2, A, B
