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

def get_bundle_spots(psf, fiber_min, fiber_max, lamp_lines, image=None, weight=None, 
                     min_dist_angstrom=0.0, sn_threshold=0.0, wave_min=None, wave_max=None,
                     broken_fibers=None):
    lines_sc = [l for l in lamp_lines if 1 <= l.get('score', 1) <= 4]
    broken_list = []
    if broken_fibers:
        if isinstance(broken_fibers, str): broken_list = [int(f) for f in broken_fibers.split(",") if f.strip()]
        else: broken_list = list(broken_fibers)
    
    nx, ny = (4114, 4128)
    if image is not None: nx, ny = image.shape
    
    spots = []
    for fiber in range(fiber_min, fiber_max + 1):
        if fiber in broken_list: continue
        trace = psf.fiber_traces[fiber]
        wmin_f, wmax_f = trace['X_vs_W'].xmin, trace['X_vs_W'].xmax
        lines_f = [l for l in lines_sc if wmin_f <= l['wave'] <= wmax_f]
        lines_f = sorted(lines_f, key=lambda x: x['wave'])
        
        for i, line in enumerate(lines_f):
            wave = line['wave']
            dist_prev = wave - lines_f[i-1]['wave'] if i > 0 else 1e9
            dist_next = lines_f[i+1]['wave'] - wave if i < len(lines_f)-1 else 1e9
            if min(dist_prev, dist_next) < min_dist_angstrom: continue
            
            xc = psf.x_ccd(fiber, wave); yc = psf.y_ccd(fiber, wave)
            if 0 <= xc < nx and 0 <= yc < ny:
                if image is not None:
                    ic, jc = int(np.floor(xc+0.5)), int(np.floor(yc+0.5))
                    i0, i1 = max(0, ic-1), min(nx, ic+2); j0, j1 = max(0, jc-1), min(ny, jc+2)
                    sig = np.sum(image[i0:i1, j0:j1])
                    var = np.sum(1.0/weight[i0:i1, j0:j1]) if weight is not None else sig
                    if var <= 0 or (sig / np.sqrt(var)) < sn_threshold: continue

                hsize_x, hsize_y = psf.h_size_x, psf.h_size_y
                spot = {
                    'fiber': fiber, 'wave': wave, 'xc_init': xc, 'yc_init': yc, 'flux': 1000.0,
                    'stamp_imin': int(np.floor(xc + 0.5)) - hsize_x, 'stamp_imax': int(np.floor(xc + 0.5)) + hsize_x + 1,
                    'stamp_jmin': int(np.floor(yc + 0.5)) - hsize_y, 'stamp_jmax': int(np.floor(yc + 0.5)) + hsize_y + 1
                }
                spots.append(spot)
    return spots

def apply_dead_column_mask(psf, fiber_min, fiber_max, weight):
    nx, ny = weight.shape
    w_new = weight.copy(); rows_j = np.arange(ny).astype(float); count = 0
    for fib in range(fiber_min, fiber_max + 1):
        try:
            w_vals = psf.fiber_traces[fib]['Y_vs_W'].invert(rows_j)
            x_vals = psf.x_ccd(fib, w_vals); i_centers = np.floor(x_vals + 0.5).astype(int)
            valid = (rows_j >= 0) & (rows_j < ny) & (i_centers >= 0) & (i_centers < nx)
            bad = weight[i_centers[valid], rows_j[valid].astype(int)] == 0
            for j in np.where(valid)[0][bad]:
                ic = i_centers[j]; i_s, i_e = max(0, ic - 4), min(nx, ic + 5)
                w_new[i_s:i_e, j] = 0.0; count += 1
        except: continue
    return w_new

def get_bundle_footprint(psf, spots, fiber_min, fiber_max, weight=None):
    margin = 0; nx, ny = (4114, 4128)
    if weight is not None: nx, ny = weight.shape
    rows_j = np.arange(ny).astype(float)
    w1 = psf.fiber_traces[fiber_min]['Y_vs_W'].invert(rows_j); w2 = psf.fiber_traces[fiber_max]['Y_vs_W'].invert(rows_j)
    x1, x2 = psf.x_ccd(fiber_min, w1), psf.x_ccd(fiber_max, w2)
    xmin_env, xmax_env = np.floor(np.minimum(x1, x2) + 0.5).astype(int) - margin, np.floor(np.maximum(x1, x2) + 0.5).astype(int) + margin + 1
    pixels = set()
    for s in spots:
        j_min, j_max = max(0, s['stamp_jmin']), min(ny, s['stamp_jmax'])
        for j in range(j_min, j_max):
            i_s, i_e = max(s['stamp_imin'], xmin_env[j]), min(s['stamp_imax'], xmax_env[j])
            for i in range(i_s, i_e):
                if 0 <= i < nx and (weight is None or weight[i, j] > 0): pixels.add((i, j))
    pixels = sorted(list(pixels)); pixels = np.array(pixels)
    if not pixels.size: return np.array([]), np.array([]), {}
    return pixels[:, 0], pixels[:, 1], {(p[0], p[1]): i for i, p in enumerate(pixels)}

# --- Standardized Coupled GPU Accumulation ---

@partial(jit, static_argnums=(11, 17, 18))
def accumulate_bundle_gpu_jnp(flux, psf_coeffs, trace_coeffs, continuum_coeffs, 
                             xc_init, yc_init, monomials, 
                             image_data, weight_data, xpix, ypix, degree,
                             stamps_x, stamps_y, stamp_indices,
                             trace_x_vals, trace_w_vals, wmin_cont, wmax_cont):
    Ns = flux.shape[0]; Npoly = monomials.shape[1]; Nparams = psf_coeffs.shape[0]; Ncont = continuum_coeffs.shape[0]
    Np = xpix.shape[0]; Nsh = (Nparams + 2) * Npoly; Ntot = Ns + Nsh + Ncont; stamp_area = stamps_x.shape[1]
    gh_all = jnp.dot(monomials, psf_coeffs.T); dx, dy = jnp.dot(monomials, trace_coeffs[0]), jnp.dot(monomials, trace_coeffs[1])
    xc_all, yc_all = xc_init + dx, yc_init + dy
    
    def get_unit_sig(i): return GaussHermitePSF.single_pix_value_jnp(xc_all[i], yc_all[i], stamps_x[i], stamps_y[i], gh_all[i], degree)
    unit_sigs = vmap(get_unit_sig)(jnp.arange(Ns))
    tsig_p = jnp.zeros(Np + 1).at[stamp_indices.flatten()].add((unit_sigs * flux[:, jnp.newaxis]).flatten())
    
    from .math import legendre_pol_jnp
    rw = 2 * (trace_w_vals - wmin_cont) / (wmax_cont - wmin_cont) - 1
    m_cont = jnp.stack([legendre_pol_jnp(k, rw) for k in range(Ncont)], axis=0)
    f_cont = jnp.tensordot(continuum_coeffs, m_cont, axes=([0], [0]))
    striped_cont = jnp.sum(vmap(lambda fi: f_cont[fi] * jnp.exp(-0.5 * (xpix - trace_x_vals[fi])**2) / jnp.sqrt(2 * jnp.pi))(jnp.arange(25)), axis=0)
    
    total_sig = tsig_p[:Np] + striped_cont; res = image_data - total_sig; chi2 = jnp.sum(weight_data * res**2)
    def get_spot_chi2(i):
        mask = (stamp_indices[i] < Np); idx = jnp.where(mask, stamp_indices[i], 0)
        return jnp.sum(weight_data[idx] * (res[idx] * mask)**2)
    chi2_spots = vmap(get_spot_chi2)(jnp.arange(Ns))

    h_cont = vmap(lambda k: jnp.sum(vmap(lambda fi: m_cont[k, fi] * jnp.exp(-0.5 * (xpix - trace_x_vals[fi])**2) / jnp.sqrt(2 * jnp.pi))(jnp.arange(25)), axis=0))(jnp.arange(Ncont)).T

    batch_size = 100; n_batches = (Ns + batch_size - 1) // batch_size; n_padded = n_batches * batch_size
    gh_p, xc_p, yc_p, f_p, m_p = jnp.pad(gh_all, ((0, n_padded - Ns), (0, 0))), jnp.pad(xc_all, (0, n_padded - Ns)), jnp.pad(yc_all, (0, n_padded - Ns)), jnp.pad(flux, (0, n_padded - Ns)), jnp.pad(monomials, ((0, n_padded - Ns), (0, 0)))
    sx_p, sy_p, idx_p, mask_p = jnp.pad(stamps_x, ((0, n_padded - Ns), (0, 0))), jnp.pad(stamps_y, ((0, n_padded - Ns), (0, 0))), jnp.pad(stamp_indices, ((0, n_padded - Ns), (0, 0)), constant_values=Np), jnp.arange(n_padded) < Ns
    
    def scan_body(carry, b_step):
        A, B = carry; start = b_step * batch_size; idx = start + jnp.arange(batch_size)
        b_xc, b_yc, b_gh, b_f, b_m, b_mask, b_sx, b_sy, b_idx = xc_p[idx], yc_p[idx], gh_p[idx], f_p[idx], m_p[idx], mask_p[idx], sx_p[idx], sy_p[idx], idx_p[idx]
        def spot_jac(x, y, g, fi, mi, s_x, s_y):
            def psf_fn(xi, yi, gi): return GaussHermitePSF.single_pix_value_jnp(xi, yi, s_x, s_y, gi, degree)
            p_s = psf_fn(x, y, g); jx, jy, jg = jacfwd(psf_fn, argnums=(0, 1, 2))(x, y, g)
            return p_s, jnp.concatenate([(fi * jg[:, :, jnp.newaxis] * mi).reshape(stamp_area, -1), (fi * jx[:, jnp.newaxis] * mi), (fi * jy[:, jnp.newaxis] * mi)], axis=1)
        b_psf, b_hsh = vmap(spot_jac)(b_xc, b_yc, b_gh, b_f, b_m, b_sx, b_sy)
        bm = b_mask[:, jnp.newaxis]; b_psf, b_hsh = b_psf * bm, b_hsh * bm[:, jnp.newaxis]
        b_res, b_w = res[jnp.where(b_idx < Np, b_idx, 0)] * (b_idx < Np), weight_data[jnp.where(b_idx < Np, b_idx, 0)] * (b_idx < Np); wr = b_w * b_res
        B = B.at[idx].add(jnp.sum(wr * b_psf, axis=1)).at[Ns:Ns+Nsh].add(jnp.sum(jnp.sum(b_hsh * wr[:, :, jnp.newaxis], axis=1), axis=0))
        A = A.at[idx, idx].add(jnp.sum(b_w * b_psf**2, axis=1)).at[Ns:Ns+Nsh, Ns:Ns+Nsh].add(jnp.sum(vmap(lambda i: jnp.dot(b_hsh[i].T * b_w[i], b_hsh[i]))(jnp.arange(batch_size)), axis=0))
        cr = vmap(lambda i: jnp.dot(b_hsh[i].T, b_w[i] * b_psf[i]))(jnp.arange(batch_size))
        A = A.at[Ns:Ns+Nsh, idx].add(cr.T).at[idx, Ns:Ns+Nsh].add(cr)
        def get_cross_sc(i): return jnp.dot(b_psf[i] * b_w[i], h_cont[jnp.where(b_idx[i] < Np, b_idx[i], 0)] * (b_idx[i] < Np)[:, jnp.newaxis])
        batch_cross_sc = vmap(get_cross_sc)(jnp.arange(batch_size))
        A = A.at[idx, -Ncont:].add(batch_cross_sc).at[-Ncont:, idx].add(batch_cross_sc.T)
        def get_cross_shc(i): return jnp.dot(b_hsh[i].T * b_w[i], h_cont[jnp.where(b_idx[i] < Np, b_idx[i], 0)] * (b_idx[i] < Np)[:, jnp.newaxis])
        batch_cross_shc = jnp.sum(vmap(get_cross_shc)(jnp.arange(batch_size)), axis=0)
        A = A.at[Ns:Ns+Nsh, -Ncont:].add(batch_cross_shc).at[-Ncont:, Ns:Ns+Nsh].add(batch_cross_shc.T)
        return (A, B), None
        
    (A_f, B_f), _ = lax.scan(scan_body, (jnp.zeros((Ntot, Ntot)), jnp.zeros(Ntot)), jnp.arange(n_batches))
    B_f = B_f.at[-Ncont:].set(jnp.dot(h_cont.T, weight_data * res)); A_f = A_f.at[-Ncont:, -Ncont:].set(jnp.dot(h_cont.T * weight_data, h_cont))
    return chi2, A_f, B_f, chi2_spots

class PSF_Fitter:
    def __init__(self, psf):
        self.psf = psf; self.chi2_precision = 1e-4
    def fit(self, image, weight, spots, bundle_id, fit_type='full', max_iter=50):
        print(f"Starting FULLY-COUPLED BATCHED-STRIPED GPU fit for bundle {bundle_id}...")
        weight = apply_dead_column_mask(self.psf, spots[0]['fiber'], spots[-1]['fiber'], weight)
        fmin, fmax = spots[0]['fiber'], spots[-1]['fiber']
        xpix, ypix, pix_idx = get_bundle_footprint(self.psf, spots, fmin, fmax, weight)
        Np = len(xpix); print(f"Footprint: {Np} pixels."); area = (2*self.psf.h_size_x+1)*(2*self.psf.h_size_y+1); Ns = len(spots)
        sx, sy = np.zeros((Ns, area)), np.zeros((Ns, area)); idx_g = np.full((Ns, area), Np, dtype=np.int32)
        for s_i, s in enumerate(spots):
            k = 0
            for j in range(s['stamp_jmin'], s['stamp_jmax']):
                for i in range(s['stamp_imin'], s['stamp_imax']):
                    sx[s_i, k], sy[s_i, k] = i, j
                    if (i, j) in pix_idx: idx_g[s_i, k] = pix_idx[(i, j)]
                    k += 1
        rows_u = np.unique(ypix); row_m = {j: i for i, j in enumerate(rows_u)}
        tx_j, tw_j = np.zeros((25, len(rows_u))), np.zeros((25, len(rows_u)))
        for f_i in range(25):
            fib = fmin + f_i; w_v = self.psf.fiber_traces[fib]['Y_vs_W'].invert(rows_u.astype(float))
            tw_j[f_i], tx_j[f_i] = w_v, self.psf.x_ccd(fib, w_v)
        ix_r = np.array([row_m[j] for j in ypix]); tx_g, tw_g = jnp.array(tx_j[:, ix_r]), jnp.array(tw_j[:, ix_r])
        flux = jnp.array([s['flux'] for s in spots]); xc_init, yc_init = jnp.array([s['xc_init'] for s in spots]), jnp.array([s['yc_init'] for s in spots])
        monomials = get_bundle_monomials_jnp(self.psf, bundle_id, spots); pc = jnp.zeros((55, monomials.shape[1])).at[0, 0].set(1.1).at[1, 0].set(1.1)
        tc = jnp.zeros((2, monomials.shape[1])); Ncont = 4; cc = jnp.zeros(Ncont) 
        img_d, w_d = jnp.array(image[xpix, ypix]), jnp.array(weight[xpix, ypix])
        wmin_c, wmax_c = float(self.psf.fiber_traces[fmin]['X_vs_W'].xmin), float(self.psf.fiber_traces[fmin]['X_vs_W'].xmax); old_chi2 = 1e30
        for i in range(max_iter):
            chi2, A, B, chi2_spots = accumulate_bundle_gpu_jnp(flux, pc, tc, cc, xc_init, yc_init, monomials, img_d, w_d, jnp.array(xpix), jnp.array(ypix), self.psf.gh_psf.degree, jnp.array(sx), jnp.array(sy), jnp.array(idx_g), tx_g, tw_g, wmin_c, wmax_c)
            mode = 'flux' if i < 2 else 'trace' if i < 7 else 'full'
            print(f"Iter {i}: chi2 = {float(chi2):.4f} [Mode: {mode}]")
            Ns_l, Npoly = len(flux), monomials.shape[1]
            if mode == 'flux': idx = jnp.concatenate([jnp.arange(Ns_l), jnp.arange(A.shape[0]-4, A.shape[0])])
            elif mode == 'trace': idx = jnp.concatenate([jnp.arange(Ns_l), jnp.arange(Ns_l+55*Npoly, Ns_l+57*Npoly), jnp.arange(A.shape[0]-4, A.shape[0])])
            else: idx = jnp.arange(A.shape[0])
            A_sub, B_sub = A[jnp.ix_(idx, idx)], B[idx]; diag = jnp.diag(A_sub); S = jnp.sqrt(diag); S = jnp.where(S < 1e-12, 1.0, S)
            A_reg = (A_sub / jnp.outer(S, S)) + 1e-4 * jnp.eye(A_sub.shape[0]); ds = jnp.linalg.solve(A_reg, B_sub / S); d_p = jnp.zeros(A.shape[0]).at[idx].set(ds / S)
            alpha = 0.5 if mode == 'full' else 1.0
            flux = jnp.maximum(flux + alpha * d_p[:Ns_l], 0.0); pc = pc + alpha * d_p[Ns_l : Ns_l+55*Npoly].reshape(55, Npoly); tc = tc + alpha * d_p[Ns_l+55*Npoly : Ns_l+57*Npoly].reshape(2, Npoly); cc = cc + alpha * d_p[-4:]
            if mode == 'full' and jnp.abs(old_chi2 - chi2) < self.chi2_precision: break
            old_chi2 = chi2
        return float(old_chi2), np.array(pc), np.array(tc), np.array(cc)
