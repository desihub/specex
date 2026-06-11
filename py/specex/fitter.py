import os
import numpy as np

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
    import jax.numpy as jnp
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
                     sn_threshold=3.15, min_dist_angstrom=0.0, wave_min=None, wave_max=None,
                     broken_fibers=None):
    """
    Match C++ Specex spot selection with correct CCD boundaries and outlier rejection.
    """
    lines_sc = [l for l in lamp_lines if 1 <= l.get('score', 1) <= 4]
    broken_list = []
    if broken_fibers:
        if isinstance(broken_fibers, str): broken_list = [int(f) for f in broken_fibers.split(",") if f.strip()]
        else: broken_list = list(broken_fibers)
    
    nx, ny = (4114, 4128)
    if image is not None: nx, ny = image.shape
    
    candidates = []
    for fiber in range(fiber_min, fiber_max + 1):
        if fiber in broken_list: continue
        if fiber not in psf.fiber_traces: continue
        for line in lines_sc:
            wave = line['wave']
            xc = psf.x_ccd(fiber, wave); yc = psf.y_ccd(fiber, wave)
            if 0 <= xc < nx and -4 <= yc < ny + 4:
                hsize_x, hsize_y = psf.h_size_x, psf.h_size_y
                candidates.append({
                    'fiber': fiber, 'wave': wave, 'xc_init': xc, 'yc_init': yc,
                    'stamp_imin': int(np.floor(xc + 0.5)) - hsize_x, 'stamp_imax': int(np.floor(xc + 0.5)) + hsize_x + 1,
                    'stamp_jmin': int(np.floor(yc + 0.5)) - hsize_y, 'stamp_jmax': int(np.floor(yc + 0.5)) + hsize_y + 1
                })
    if not candidates: return []
    for spot in candidates:
        if image is None: 
            spot['snr'] = 100.0; spot['flux'] = 1000.0; spot['chi2'] = 0.0
            continue
        i0, i1 = max(0, spot['stamp_imin']), min(nx, spot['stamp_imax'])
        j0, j1 = max(0, spot['stamp_jmin']), min(ny, spot['stamp_jmax'])
        gh = psf.gh_params(spot['fiber'], spot['wave'])
        ix, iy = np.meshgrid(np.arange(i0, i1), np.arange(j0, j1), indexing='ij')
        psf_val = GaussHermitePSF.single_pix_value_np(spot['xc_init'], spot['yc_init'], ix.flatten(), iy.flatten(), gh, psf.gh_psf.degree)
        img_val = image[i0:i1, j0:j1].flatten()
        w_val = weight[i0:i1, j0:j1].flatten() if weight is not None else 1.0
        A = np.sum(w_val * psf_val**2); B = np.sum(w_val * img_val * psf_val)
        if A > 0:
            flux = B/A; spot['flux'] = flux; spot['snr'] = flux/np.sqrt(1.0/A)
            res = img_val - flux * psf_val; spot['chi2'] = np.sum(w_val * res**2)
        else:
            spot['flux'] = 0.0; spot['snr'] = -1.0; spot['chi2'] = 1e10
    nsig = 4.0
    for i, s in enumerate(candidates):
        if s['chi2'] > 1e9: continue
        others = [candidates[j]['chi2'] for j in range(len(candidates)) if i != j and abs(candidates[j]['wave'] - s['wave']) < 1.0]
        if len(others) < 2: continue
        m_c2 = np.mean(others); rms_c2 = np.std(others)
        if s['chi2'] > (m_c2 + nsig * rms_c2): s['snr'] = -2.0
    selected = []
    for i, s in enumerate(candidates):
        if s['snr'] < sn_threshold: continue
        if min_dist_angstrom > 0:
            if any(abs(o['wave'] - s['wave']) < min_dist_angstrom for j, o in enumerate(candidates) if i != j and o['fiber'] == s['fiber']):
                continue
        selected.append(s)
    return selected

def apply_dead_column_mask(psf, fiber_min, fiber_max, weight):
    nx, ny = weight.shape
    w_new = weight.copy(); rows_j = np.arange(ny).astype(float)
    for fib in range(fiber_min, fiber_max + 1):
        try:
            w_vals = psf.fiber_traces[fib]['Y_vs_W'].invert(rows_j)
            x_vals = psf.x_ccd(fib, w_vals); i_centers = np.floor(x_vals + 0.5).astype(int)
            valid = (rows_j >= 0) & (rows_j < ny) & (i_centers >= 0) & (i_centers < nx)
            bad = weight[i_centers[valid], rows_j[valid].astype(int)] == 0
            for j in np.where(valid)[0][bad]:
                ic = i_centers[j]; i_s, i_e = max(0, ic - 4), min(nx, ic + 5)
                w_new[i_s:i_e, j] = 0.0
        except: continue
    return w_new

def get_bundle_footprint(psf, spots, fiber_min, fiber_max, weight=None):
    nx, ny = (4114, 4128)
    if weight is not None: nx, ny = weight.shape
    rows_j = np.arange(ny).astype(float)
    w1 = psf.fiber_traces[fiber_min]['Y_vs_W'].invert(rows_j); w2 = psf.fiber_traces[fiber_max]['Y_vs_W'].invert(rows_j)
    x1, x2 = psf.x_ccd(fiber_min, w1), psf.x_ccd(fiber_max, w2)
    xmin_env, xmax_env = np.floor(np.minimum(x1, x2) + 0.5).astype(int), np.floor(np.maximum(x1, x2) + 0.5).astype(int) + 1
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

# --- High-Performance Analytical GPU Driver ---

def accumulate_bundle_gpu_jnp(flux, psf_coeffs, trace_coeffs, continuum_coeffs, 
                             xc_init, yc_init, monomials, 
                             image_data, weight_data, xpix, ypix, degree,
                             stamps_x, stamps_y, stamp_indices,
                             trace_x_vals, trace_w_vals, wmin_cont, wmax_cont):
    import jax
    import jax.numpy as jnp
    from jax import jit, vmap, lax
    from functools import partial
    
    @partial(jit, static_argnums=(11, 17, 18))
    def _accumulate(flux, psf_coeffs, trace_coeffs, continuum_coeffs, 
                    xc_init, yc_init, monomials, 
                    image_data, weight_data, xpix, ypix, degree,
                    stamps_x, stamps_y, stamp_indices,
                    trace_x_vals, trace_w_vals, wmin_cont, wmax_cont):
        Ns = flux.shape[0]; Npoly = monomials.shape[1]; Nparams = psf_coeffs.shape[0]; Ncont = continuum_coeffs.shape[0]
        Np = xpix.shape[0]; Nsh = (Nparams + 2) * Npoly; Ntot = Ns + Nsh + Ncont; stamp_area = stamps_x.shape[1]
        
        gh_all = jnp.dot(monomials, psf_coeffs.T); dx, dy = jnp.dot(monomials, trace_coeffs[0]), jnp.dot(monomials, trace_coeffs[1])
        xc_all, yc_all = xc_init + dx, yc_init + dy
        
        batch_size = 2000; n_pad = batch_size - Ns
        gh_p = jnp.pad(gh_all, ((0, n_pad), (0, 0))); xc_p, yc_p, f_p = jnp.pad(xc_all, (0, n_pad)), jnp.pad(yc_all, (0, n_pad)), jnp.pad(flux, (0, n_pad))
        m_p = jnp.pad(monomials, ((0, n_pad), (0, 0))); sx_p = jnp.pad(stamps_x, ((0, n_pad), (0, 0))); sy_p = jnp.pad(stamps_y, ((0, n_pad), (0, 0)))
        idx_p = jnp.pad(stamp_indices, ((0, n_pad), (0, 0)), constant_values=Np); mask_p = jnp.arange(batch_size) < Ns

        def pix_fn_scalar(x, y, g, px, py): return GaussHermitePSF.single_pix_value_jnp(x, y, px, py, g, degree)
        val_grad_all = vmap(vmap(jax.value_and_grad(pix_fn_scalar, argnums=(0, 1, 2)), in_axes=(None, None, None, 0, 0)), in_axes=(0, 0, 0, 0, 0))
        
        b_psf, (b_jx, b_jy, b_jg) = val_grad_all(xc_p, yc_p, gh_p, sx_p, sy_p)
        bm = mask_p[:, jnp.newaxis]; b_psf *= bm
        
        j_psf = (f_p[:, jnp.newaxis, jnp.newaxis, jnp.newaxis] * b_jg[:, :, :, jnp.newaxis] * m_p[:, jnp.newaxis, jnp.newaxis, :]).reshape(batch_size, stamp_area, -1)
        j_xc = f_p[:, jnp.newaxis, jnp.newaxis] * b_jx[:, :, jnp.newaxis] * m_p[:, jnp.newaxis, :]
        j_yc = f_p[:, jnp.newaxis, jnp.newaxis] * b_jy[:, :, jnp.newaxis] * m_p[:, jnp.newaxis, :]
        b_jac = jnp.concatenate([j_psf, j_xc, j_yc], axis=2) * bm[:, jnp.newaxis]

        flat_idx = idx_p.flatten(); valid = flat_idx < Np
        total_sig = jnp.zeros(Np + 1).at[flat_idx].add((b_psf * f_p[:, jnp.newaxis]).flatten())
        
        from .math import legendre_pol_jnp
        rw = 2 * (trace_w_vals - wmin_cont) / (wmax_cont - wmin_cont) - 1
        m_cont = jnp.stack([legendre_pol_jnp(k, rw) for k in range(Ncont)], axis=0)
        f_cont = jnp.tensordot(continuum_coeffs, m_cont, axes=([0], [0]))
        striped_cont = jnp.sum(vmap(lambda fi: f_cont[fi] * jnp.exp(-0.5 * (xpix - trace_x_vals[fi])**2) / jnp.sqrt(2 * jnp.pi))(jnp.arange(25)), axis=0)
        res = image_data - (total_sig[:Np] + striped_cont); chi2 = jnp.sum(weight_data * res**2)
        h_cont = vmap(lambda k: jnp.sum(vmap(lambda fi: m_cont[k, fi] * jnp.exp(-0.5 * (xpix - trace_x_vals[fi])**2) / jnp.sqrt(2 * jnp.pi))(jnp.arange(25)), axis=0))(jnp.arange(Ncont)).T

        b_res = jnp.where(valid, res[jnp.where(valid, flat_idx, 0)], 0.0).reshape(batch_size, stamp_area)
        b_w = jnp.where(valid, weight_data[jnp.where(valid, flat_idx, 0)], 0.0).reshape(batch_size, stamp_area)
        wr = b_w * b_res
        
        B = jnp.zeros(Ntot).at[:Ns].set(jnp.sum(wr * b_psf, axis=1)[:Ns]).at[Ns:Ns+Nsh].set(jnp.sum(jnp.sum(b_jac * wr[:, :, jnp.newaxis], axis=1), axis=0)).at[-Ncont:].set(jnp.dot(h_cont.T, weight_data * res))
        A = jnp.zeros((Ntot, Ntot)).at[jnp.arange(Ns), jnp.arange(Ns)].set(jnp.sum(b_w * b_psf**2, axis=1)[:Ns])
        A = A.at[Ns:Ns+Nsh, Ns:Ns+Nsh].set(jnp.einsum('bij,bi,bik->jk', b_jac, b_w, b_jac))
        A_fs = jnp.einsum('bij,bi,bi->bj', b_jac, b_w, b_psf); A = A.at[Ns:Ns+Nsh, :Ns].set(A_fs[:Ns].T).at[:Ns, Ns:Ns+Nsh].set(A_fs[:Ns])
        h_lookup = h_cont[jnp.where(valid, flat_idx, 0)]; b_hcont = jnp.where(valid[:, jnp.newaxis], h_lookup, 0.0).reshape(batch_size, stamp_area, Ncont)
        A_fc = jnp.einsum('bi,bi,bik->bk', b_psf, b_w, b_hcont); A = A.at[:Ns, -Ncont:].set(A_fc[:Ns]).at[-Ncont:, :Ns].set(A_fc[:Ns].T)
        A_sc = jnp.einsum('bij,bi,bik->jk', b_jac, b_w, b_hcont); A = A.at[Ns:Ns+Nsh, -Ncont:].set(A_sc).at[-Ncont:, Ns:Ns+Nsh].set(A_sc.T)
        A = A.at[-Ncont:, -Ncont:].set(jnp.dot(h_cont.T * weight_data, h_cont))
        return chi2, A[:Ns+Nsh+Ncont, :Ns+Nsh+Ncont], B[:Ns+Nsh+Ncont], None

    return _accumulate(flux, psf_coeffs, trace_coeffs, continuum_coeffs, xc_init, yc_init, monomials, image_data, weight_data, xpix, ypix, degree, stamps_x, stamps_y, stamp_indices, trace_x_vals, trace_w_vals, wmin_cont, wmax_cont)

class PSF_Fitter:
    def __init__(self, psf):
        self.psf = psf; self.chi2_precision = 10.0
    def fit(self, image, weight, spots, bundle_id, fit_type='full', max_iter=15):
        import jax.numpy as jnp
        from jax import jit, vmap
        print(f"Starting HIGH-PERFORMANCE OPTIMIZED COUPLED GPU fit for bundle {bundle_id}...")
        weight = apply_dead_column_mask(self.psf, spots[0]['fiber'], spots[-1]['fiber'], weight)
        fmin, fmax = spots[0]['fiber'], spots[-1]['fiber']
        xpix, ypix, pix_idx = get_bundle_footprint(self.psf, spots, fmin, fmax, weight)
        Np = len(xpix); area = (2*self.psf.h_size_x+1)*(2*self.psf.h_size_y+1); Ns = len(spots)
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
        sx_g, sy_g, idx_gg = jnp.array(sx), jnp.array(sy), jnp.array(idx_g)
        @jit
        def get_chi2(f, p, t, c):
            Ns_l = f.shape[0]; gh_all = jnp.dot(monomials, p.T); dx, dy = jnp.dot(monomials, t[0]), jnp.dot(monomials, t[1]); xc_all, yc_all = xc_init + dx, yc_init + dy
            unit_sigs = vmap(lambda i: GaussHermitePSF.single_pix_value_jnp(xc_all[i], yc_all[i], sx_g[i], sy_g[i], gh_all[i], self.psf.gh_psf.degree))(jnp.arange(Ns_l))
            tsig_p = jnp.zeros(Np + 1).at[idx_gg.flatten()].add((unit_sigs * f[:, jnp.newaxis]).flatten())
            from .math import legendre_pol_jnp
            rw = 2 * (tw_g - wmin_c) / (wmax_c - wmin_c) - 1
            m_c = jnp.stack([legendre_pol_jnp(k, rw) for k in range(Ncont)], axis=0); f_c = jnp.tensordot(c, m_c, axes=([0], [0]))
            striped_cont = jnp.sum(vmap(lambda fi: f_c[fi] * jnp.exp(-0.5 * (xpix - tx_g[fi])**2) / jnp.sqrt(2 * jnp.pi))(jnp.arange(25)), axis=0)
            total_sig = tsig_p[:Np] + striped_cont; return jnp.sum(w_d * (img_d - total_sig)**2)
        for i in range(max_iter):
            chi2, A, B, _ = accumulate_bundle_gpu_jnp(flux, pc, tc, cc, xc_init, yc_init, monomials, img_d, w_d, jnp.array(xpix), jnp.array(ypix), self.psf.gh_psf.degree, sx_g, sy_g, idx_gg, tx_g, tw_g, wmin_c, wmax_c)
            mode = 'flux' if i < 2 else 'trace' if i < 5 else 'full'
            print(f"Iter {i}: chi2 = {float(chi2):.4f} [Mode: {mode}]", flush=True)
            Npoly = monomials.shape[1]; Ns_l = len(flux)
            if mode == 'flux': idx = jnp.concatenate([jnp.arange(Ns_l), jnp.arange(A.shape[0]-Ncont, A.shape[0])])
            elif mode == 'trace': idx = jnp.concatenate([jnp.arange(Ns_l), jnp.arange(Ns_l+55*Npoly, Ns_l+57*Npoly), jnp.arange(A.shape[0]-Ncont, A.shape[0])])
            else: idx = jnp.arange(A.shape[0])
            A_sub, B_sub = A[jnp.ix_(idx, idx)], B[idx]; diag = jnp.diag(A_sub); S = jnp.sqrt(diag); S = jnp.where(S < 1e-12, 1.0, S)
            A_reg = (A_sub / jnp.outer(S, S)) + 1e-4 * jnp.eye(A_sub.shape[0]); ds = jnp.linalg.solve(A_reg, B_sub / S); d_p = jnp.zeros(A.shape[0]).at[idx].set(ds / S)
            best_alpha, best_chi2 = 0.0, float(chi2)
            for alpha in [0.2, 0.5, 1.0, 1.5, 2.0]:
                f_try = jnp.maximum(flux + alpha * d_p[:Ns_l], 0.0); p_try = pc + alpha * d_p[Ns_l : Ns_l+55*Npoly].reshape(55, Npoly); t_try = tc + alpha * d_p[Ns_l+55*Npoly : Ns_l+57*Npoly].reshape(2, Npoly); c_try = cc + alpha * d_p[-Ncont:]; c2 = get_chi2(f_try, p_try, t_try, c_try)
                if c2 < best_chi2: best_alpha, best_chi2 = alpha, c2
            if best_alpha == 0: flux = jnp.maximum(flux + 0.1 * d_p[:Ns_l], 0.0); best_alpha = 0.1
            flux = jnp.maximum(flux + best_alpha * d_p[:Ns_l], 0.0); pc = pc + best_alpha * d_p[Ns_l : Ns_l+55*Npoly].reshape(55, Npoly); tc = tc + best_alpha * d_p[Ns_l+55*Npoly : Ns_l+57*Npoly].reshape(2, Npoly); cc = cc + best_alpha * d_p[-Ncont:]
            if mode == 'full' and jnp.abs(old_chi2 - chi2) < self.chi2_precision: break
            old_chi2 = chi2
        return float(old_chi2), np.array(pc), np.array(tc), np.array(cc)
