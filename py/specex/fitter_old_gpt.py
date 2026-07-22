import os
import time
import numpy as np
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

# --- Global JIT Kernels (Standardized signatures) ---

def _predict_bundle_jax(flux, psf_coeffs, trace_coeffs, continuum_coeffs,
                        xc_init, yc_init, monomials, xpix, ypix,
                        sx_g, sy_g, idx_gg, degree,
                        tx_g, tw_g, wmin_c, wmax_c, image_data, weight_data):
    import jax.numpy as jnp
    from jax import vmap
    from .math import legendre_pol_jnp
    
    Ns = flux.shape[0]; Npoly = monomials.shape[1]; Np = xpix.shape[0]
    gh_all = jnp.dot(monomials, psf_coeffs.T)
    dx, dy = jnp.dot(monomials, trace_coeffs[0]), jnp.dot(monomials, trace_coeffs[1])
    xc_all, yc_all = xc_init + dx, yc_init + dy
    
    def spot_sig(i):
        Bx, By = GaussHermitePSF.get_gh_basis(xc_all[i], yc_all[i], sx_g[i], sy_g[i], gh_all[i], degree)
        psf_v = By[0] * Bx[0] 
        nx_p, ny_p = degree + 1, degree + 1; k = 2
        for j_p in range(ny_p):
            bj = By[j_p]; imin_p = 1 if j_p == 0 else 0
            for i_p in range(imin_p, nx_p):
                psf_v += gh_all[i, k] * bj * Bx[i_p]; k += 1
        return psf_v
    
    unit_sigs = vmap(spot_sig)(jnp.arange(Ns))
    tsig_p = jnp.zeros(Np + 1).at[idx_gg.flatten()].add((unit_sigs * flux[:, jnp.newaxis]).flatten())
    
    Ncont = continuum_coeffs.shape[0]
    rw = 2 * (tw_g - wmin_c) / (wmax_c - wmin_c) - 1
    m_c = jnp.stack([legendre_pol_jnp(k, rw) for k in range(Ncont)], axis=0)
    f_c = jnp.tensordot(continuum_coeffs, m_c, axes=([0], [0]))
    striped_cont = jnp.sum(vmap(lambda fi: f_c[fi] * jnp.exp(-0.5 * (xpix - tx_g[fi])**2) / jnp.sqrt(2 * jnp.pi))(jnp.arange(25)), axis=0)
    
    total_sig = tsig_p[:Np] + striped_cont
    return jnp.sum(weight_data * (image_data - total_sig)**2)

from jax import jit
_predict_bundle_jax_jit = jit(_predict_bundle_jax, static_argnums=(12,))

def _accumulate_bundle_jax(flux, psf_coeffs, trace_coeffs, continuum_coeffs, 
                               xc_init, yc_init, monomials, xpix, ypix,
                               sx_g, sy_g, idx_gg, degree,
                               tx_g, tw_g, wmin_c, wmax_c, image_data, weight_data,
                               gain, psf_error, wscale):

    import jax
    import jax.numpy as jnp
    from jax import vmap, lax
    
    Ns = flux.shape[0]; Npoly = monomials.shape[1]; Np = xpix.shape[0]; Ncont = continuum_coeffs.shape[0]
    Nparams = psf_coeffs.shape[0]; Nsh = (Nparams + 2) * Npoly; Ntot = Ns + Nsh + Ncont; stamp_area = sx_g.shape[1]
    
    gh_all = jnp.dot(monomials, psf_coeffs.T); dx, dy = jnp.dot(monomials, trace_coeffs[0]), jnp.dot(monomials, trace_coeffs[1])
    xc_all, yc_all = xc_init + dx, yc_init + dy
    
    batch_size = 2000; n_pad = batch_size - Ns
    gh_p = jnp.pad(gh_all, ((0, n_pad), (0, 0))); xc_p, yc_p, f_p = jnp.pad(xc_all, (0, n_pad)), jnp.pad(yc_all, (0, n_pad)), jnp.pad(flux, (0, n_pad))
    m_p = jnp.pad(monomials, ((0, n_pad), (0, 0))); sx_p = jnp.pad(sx_g, ((0, n_pad), (0, 0))); sy_p = jnp.pad(sy_g, ((0, n_pad), (0, 0)))
    idx_p = jnp.pad(idx_gg, ((0, n_pad), (0, 0)), constant_values=Np); mask_p = jnp.arange(batch_size) < Ns

    def get_all_grads(xi, yi, gi, si_x, si_y):
        sigx = jnp.maximum(gi[0], 0.1); sigy = jnp.maximum(gi[1], 0.1)
        isx = 1.0 / sigx; isy = 1.0 / sigy
        x1 = (jnp.floor(si_x + 0.5) - xi - 0.5) * isx; x2 = (jnp.floor(si_x + 0.5) - xi + 0.5) * isx
        y1 = (jnp.floor(si_y + 0.5) - yi - 0.5) * isy; y2 = (jnp.floor(si_y + 0.5) - yi + 0.5) * isy
        isq2 = 1.0 / jnp.sqrt(2.0); isq2pi = 1.0 / jnp.sqrt(2.0 * jnp.pi)
        gx1 = isq2pi * isx * jnp.exp(-0.5 * x1**2); gx2 = isq2pi * isx * jnp.exp(-0.5 * x2**2)
        gy1 = isq2pi * isy * jnp.exp(-0.5 * y1**2); gy2 = isq2pi * isy * jnp.exp(-0.5 * y2**2)
        ex = 0.5 * (jax.scipy.special.erf(x2 * isq2) - jax.scipy.special.erf(x1 * isq2))
        ey = 0.5 * (jax.scipy.special.erf(y2 * isq2) - jax.scipy.special.erf(y1 * isq2))
        
        from .psf import hermite_pol_jnp
        nx_p = degree + 1; ny_p = degree + 1
        H1_u = jnp.stack([hermite_pol_jnp(n, x1) for n in range(nx_p)], axis=0)
        H2_u = jnp.stack([hermite_pol_jnp(n, x2) for n in range(nx_p)], axis=0)
        H1_v = jnp.stack([hermite_pol_jnp(n, y1) for n in range(ny_p)], axis=0)
        H2_v = jnp.stack([hermite_pol_jnp(n, y2) for n in range(ny_p)], axis=0)
        
        def get_P(n, val, g1, g2, h1, h2, sig):
            return jnp.where(n == 0, val, sig * (g1 * h1[n-1] - g2 * h2[n-1]))
        
        Bx = vmap(lambda n: get_P(n, ex, gx1, gx2, H1_u, H2_u, sigx))(jnp.arange(nx_p))
        By = vmap(lambda n: get_P(n, ey, gy1, gy2, H1_v, H2_v, sigy))(jnp.arange(ny_p))
        dexdsx = (x1 * gx1 - x2 * gx2); deydsy = (y1 * gy1 - y2 * gy2); dexdx = (gx1 - gx2); deydy = (gy1 - gy2)
        
        def get_dPds(n, g1, g2, h1, h2, x1, x2, val, dvalds, sig):
            h1_m1 = h1[jnp.maximum(n-1, 0)]; h2_m1 = h2[jnp.maximum(n-1, 0)]
            h1_m2 = h1[jnp.maximum(n-2, 0)]; h2_m2 = h2[jnp.maximum(n-2, 0)]
            t1 = sig * g1 * h1_m1; t2 = sig * g2 * h2_m1
            h_prime1 = (n-1)*h1_m2; h_prime2 = (n-1)*h2_m2
            res = (-g1*x1*h_prime1 + g2*x2*h_prime2) + (t1*x1*x1/sig - t2*x2*x2/sig)
            return jnp.where(n == 0, dvalds, res)

        def get_dPdx(n, g1, g2, h1, h2, x1, x2, val, dvaldx, sig):
            h1_m1 = h1[jnp.maximum(n-1, 0)]; h2_m1 = h2[jnp.maximum(n-1, 0)]
            h1_m2 = h1[jnp.maximum(n-2, 0)]; h2_m2 = h2[jnp.maximum(n-2, 0)]
            t1 = sig * g1 * h1_m1; t2 = sig * g2 * h2_m1
            h_prime1 = (n-1)*h1_m2; h_prime2 = (n-1)*h2_m2
            res = -1.0/sig * (sig*g1*h_prime1 - sig*g2*h_prime2 - x1*t1 + x2*t2)
            return jnp.where(n == 0, dvaldx, res)

        dBxdsx = vmap(lambda n: get_dPds(n, gx1, gx2, H1_u, H2_u, x1, x2, ex, dexdsx, sigx))(jnp.arange(nx_p))
        dBydsy = vmap(lambda n: get_dPds(n, gy1, gy2, H1_v, H2_v, y1, y2, ey, deydsy, sigy))(jnp.arange(ny_p))
        dBxdx = vmap(lambda n: get_dPdx(n, gx1, gx2, H1_u, H2_u, x1, x2, ex, dexdx, sigx))(jnp.arange(nx_p))
        dBydy = vmap(lambda n: get_dPdx(n, gy1, gy2, H1_v, H2_v, y1, y2, ey, deydy, sigy))(jnp.arange(ny_p))

        psf_v = By[0] * Bx[0]; basis_gh = []; k = 2
        dsigx = dBxdsx[0] * By[0]; dsigy = Bx[0] * dBydsy[0]
        dxc = dBxdx[0] * By[0]; dyc = Bx[0] * dBydy[0]
        for j in range(ny_p):
            bj = By[j]; dbjdsy = dBydsy[j]; dbjdy = dBydy[j]
            imin = 1 if j == 0 else 0
            for i in range(imin, nx_p):
                bi = Bx[i]; dbidsx = dBxdsx[i]; dbidx = dBxdx[i]
                term = bj * bi; psf_v += gi[k] * term; basis_gh.append(term)
                dsigx += gi[k] * bj * dbidsx; dsigy += gi[k] * bi * dbjdsy
                dxc += gi[k] * bj * dbidx; dyc += gi[k] * bi * dbjdy
                k += 1
        return psf_v, jnp.stack(basis_gh, axis=1), dsigx, dsigy, dxc, dyc

    b_psf, b_gh_basis, b_jsigx, b_jsigy, b_jx, b_jy = vmap(get_all_grads)(xc_p, yc_p, gh_p, sx_p, sy_p)
    bm = mask_p[:, jnp.newaxis]; b_psf *= bm
    j_sx = (f_p[:, jnp.newaxis, jnp.newaxis] * b_jsigx[:, :, jnp.newaxis] * m_p[:, jnp.newaxis, :])
    j_sy = (f_p[:, jnp.newaxis, jnp.newaxis] * b_jsigy[:, :, jnp.newaxis] * m_p[:, jnp.newaxis, :])
    j_gh = (f_p[:, jnp.newaxis, jnp.newaxis, jnp.newaxis] * b_gh_basis[:, :, :, jnp.newaxis] * m_p[:, jnp.newaxis, jnp.newaxis, :]).reshape(batch_size, stamp_area, -1)
    j_xc = (f_p[:, jnp.newaxis, jnp.newaxis] * b_jx[:, :, jnp.newaxis] * m_p[:, jnp.newaxis, :])
    j_yc = (f_p[:, jnp.newaxis, jnp.newaxis] * b_jy[:, :, jnp.newaxis] * m_p[:, jnp.newaxis, :])
    b_jac = jnp.concatenate([j_sx, j_sy, j_gh, j_xc, j_yc], axis=2) * bm[:, jnp.newaxis]

    flat_idx = idx_p.flatten(); valid = flat_idx < Np
    total_sig = jnp.zeros(Np + 1).at[flat_idx].add((b_psf * f_p[:, jnp.newaxis]).flatten())
    from .math import legendre_pol_jnp
    rw = 2 * (tw_g - wmin_c) / (wmax_c - wmin_c) - 1
    m_cont = jnp.stack([legendre_pol_jnp(k, rw) for k in range(Ncont)], axis=0)
    f_cont = jnp.tensordot(continuum_coeffs, m_cont, axes=([0], [0]))
    striped_cont = jnp.sum(vmap(lambda fi: f_cont[fi] * jnp.exp(-0.5 * (xpix - tx_g[fi])**2) / jnp.sqrt(2 * jnp.pi))(jnp.arange(25)), axis=0)
    res = image_data - (total_sig[:Np] + striped_cont); chi2 = jnp.sum(weight_data * res**2)
    h_cont = vmap(lambda k: jnp.sum(vmap(lambda fi: m_cont[k, fi] * jnp.exp(-0.5 * (xpix - tx_g[fi])**2) / jnp.sqrt(2 * jnp.pi))(jnp.arange(25)), axis=0))(jnp.arange(Ncont)).T

    b_res = jnp.where(valid, res[jnp.where(valid, flat_idx, 0)], 0.0).reshape(batch_size, stamp_area)
    b_w = jnp.where(valid, weight_data[jnp.where(valid, flat_idx, 0)], 0.0).reshape(batch_size, stamp_area); wr = b_w * b_res
    
    # --- C++ Parity: B-vector Poisson correction (lines 670-673 in specex_psf_fitter.cc) ---
    # bfact = w*res + (1/wscale)*0.5*(w*res)^2 * (1/gain + 2*psf_error^2*signal)
    # We need 'signal' for the correction term. signal = total_sig[:Np]
    signal_vec = total_sig[:Np]
    # we need to map this signal back to the batch/stamp structure
    b_signal = jnp.where(valid, signal_vec[jnp.where(valid, flat_idx, 0)], 0.0).reshape(batch_size, stamp_area)
    
    # Use actual values passed into the JIT function
    gain = gain
    psf_error = psf_error
    wscale = wscale
    
    # The correction only applies if recompute_weight_in_fit is True.
    # In Python, we implement the correction directly into the residual product.
    correction = (1.0 / wscale) * 0.5 * (wr**2) * (1.0/gain + 2.0 * psf_error**2 * b_signal)
    wr_corrected = wr + correction
    
    B = jnp.zeros(Ntot).at[:Ns].set(jnp.sum(wr * b_psf, axis=1)[:Ns]).at[Ns:Ns+Nsh].set(jnp.sum(jnp.sum(b_jac * wr_corrected[:, :, jnp.newaxis], axis=1), axis=0)).at[-Ncont:].set(jnp.dot(h_cont.T, weight_data * res))
    A = jnp.zeros((Ntot, Ntot)).at[jnp.arange(Ns), jnp.arange(Ns)].set(jnp.sum(b_w * b_psf**2, axis=1)[:Ns])
    A = A.at[Ns:Ns+Nsh, Ns:Ns+Nsh].set(jnp.einsum('bij,bi,bik->jk', b_jac, b_w, b_jac))
    A_fs = jnp.einsum('bij,bi,bi->bj', b_jac, b_w, b_psf); A = A.at[Ns:Ns+Nsh, :Ns].set(A_fs[:Ns].T).at[:Ns, Ns:Ns+Nsh].set(A_fs[:Ns])
    h_lookup = h_cont[jnp.where(valid, flat_idx, 0)]; b_hcont = jnp.where(valid[:, jnp.newaxis], h_lookup, 0.0).reshape(batch_size, stamp_area, Ncont)
    A_fc = jnp.einsum('bi,bi,bik->bk', b_psf, b_w, b_hcont); A = A.at[:Ns, -Ncont:].set(A_fc[:Ns]).at[-Ncont:, :Ns].set(A_fc[:Ns].T)
    A_sc = jnp.einsum('bij,bi,bik->jk', b_jac, b_w, b_hcont); A = A.at[Ns:Ns+Nsh, -Ncont:].set(A_sc).at[-Ncont:, Ns:Ns+Nsh].set(A_sc.T)
    A = A.at[-Ncont:, -Ncont:].set(jnp.dot(h_cont.T * weight_data, h_cont))
    return chi2, A[:Ns+Nsh+Ncont, :Ns+Nsh+Ncont], B[:Ns+Nsh+Ncont]

_accumulate_bundle_jax_jit = jit(_accumulate_bundle_jax, static_argnums=(12, 19, 20, 21))

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
    t0 = time.time(); nx, ny = (4114, 4128)
    if weight is not None: nx, ny = weight.shape
    rows_j = np.arange(ny).astype(float)
    w1 = psf.fiber_traces[fiber_min]['Y_vs_W'].invert(rows_j); w2 = psf.fiber_traces[fiber_max]['Y_vs_W'].invert(rows_j)
    x1, x2 = psf.x_ccd(fiber_min, w1), psf.x_ccd(fiber_max, w2)
    xmin_env, xmax_env = np.floor(np.minimum(x1, x2) + 0.5).astype(int), np.floor(np.maximum(x1, x2) + 0.5).astype(int) + 1
    j_min_all = max(0, min(s['stamp_jmin'] for s in spots)); j_max_all = min(ny, max(s['stamp_jmax'] for s in spots))
    i_min_all = max(0, min(s['stamp_imin'] for s in spots)); i_max_all = min(nx, max(s['stamp_imax'] for s in spots))
    active_mask = np.zeros((i_max_all - i_min_all, j_max_all - j_min_all), dtype=bool)
    for s in spots:
        i0, i1 = max(i_min_all, s['stamp_imin']), min(i_max_all, s['stamp_imax'])
        j0, j1 = max(j_min_all, s['stamp_jmin']), min(j_max_all, s['stamp_jmax'])
        active_mask[i0 - i_min_all : i1 - i_min_all, j0 - j_min_all : j1 - j_min_all] = True
    ii, jj = np.where(active_mask); global_i, global_j = ii + i_min_all, jj + j_min_all
    valid_env = (global_i >= xmin_env[global_j]) & (global_i < xmax_env[global_j])
    if weight is not None:
        valid_weight = weight[global_i[valid_env], global_j[valid_env]] > 0
        fi, fj = global_i[valid_env][valid_weight], global_j[valid_env][valid_weight]
    else: fi, fj = global_i[valid_env], global_j[valid_env]
    pixels = np.stack([fi, fj], axis=1); pixels = pixels[np.lexsort((pixels[:, 1], pixels[:, 0]))]
    print(f"  Footprint generation took {time.time() - t0:.2f}s ({len(pixels)} pixels)", flush=True)
    return pixels[:, 0], pixels[:, 1], {(p[0], p[1]): i for i, p in enumerate(pixels)}

def get_bundle_spots(psf, fiber_min, fiber_max, lamp_lines, image=None, weight=None, 
                      sn_threshold=3.0, min_dist_angstrom=0.0, wave_min=None, wave_max=None,
                      broken_fibers=None, max_number_of_lines=100, provided_candidates=None):
    import jax.numpy as jnp
    t0 = time.time(); lines_sc = [l for l in lamp_lines if 1 <= l.get('score', 1) <= 4]
    broken_list = []
    if broken_fibers:
        if isinstance(broken_fibers, str): broken_list = [int(f) for f in broken_fibers.split(",") if f.strip()]
        else: broken_list = list(broken_fibers)
    nx, ny = (4114, 4128)
    if image is not None: nx, ny = image.shape
    
    if provided_candidates is not None:
        candidates = provided_candidates
    else:
        candidates = []
        for fiber in range(fiber_min, fiber_max + 1):
            if fiber in broken_list or fiber not in psf.fiber_traces: continue
            for line in lines_sc:
                wave = line['wave']; xc = psf.x_ccd(fiber, wave); yc = psf.y_ccd(fiber, wave)
                if 0 <= xc < nx and 0 <= yc < ny:
                    candidates.append({'fiber': fiber, 'wave': wave, 'xc_init': xc, 'yc_init': yc})

    
    print(f"PYTHON SPOT SELECTION: Starting first pass. Initial candidates: {len(candidates)}", flush=True)

    if not candidates: return []

    if hasattr(psf, 'output_psf_path') and psf.output_psf_path:
        raw_path = psf.output_psf_path.replace('.fits', '.pyrawspots.txt')
        with open(raw_path, 'w') as f:
            for s in candidates:
                f.write(f"{s['fiber']},{s['wave']:.15f},{s['xc_init']:.15f},{s['yc_init']:.15f}\n")
        print(f"  Written {len(candidates)} raw candidates to {raw_path}", flush=True)

    Ns = len(candidates)
    # Preserve the full candidate list for the distance‑check (C++ loops over all input_spots)
    all_candidates = list(candidates)
    if image is not None:
        c_xc = jnp.array([s['xc_init'] for s in candidates]); c_yc = jnp.array([s['yc_init'] for s in candidates])
        gh_all = jnp.array([psf.gh_params(s['fiber'], s['wave']) for s in candidates])
        # C++ Logic: Fit each spot first to get robust flux/S/N (mirrors FitIndividualSpotFluxes)
        # In our current JAX implementation, _get_spot_stats_jax already performs the a‑priori 
        # flux estimation (B/A). To strictly mirror C++, we also return the weight‑sum A so we can
        # apply the exact eflux ≤ 0 check used by the C++ code.
        
        fluxes, snrs, chi2s, As = _get_spot_stats_jax(jnp.array(image), jnp.array(weight), c_xc, c_yc, gh_all, psf.gh_psf.degree, psf.h_size_x, psf.h_size_y)
        fluxes, snrs, chi2s, As = np.array(fluxes), np.array(snrs), np.array(chi2s), np.array(As)
        waves = np.array([s['wave'] for s in candidates])
        
        # C++ Logic: 5x5 Zero-Weight Check (Strict Flux Measurement)
        # Spots are ignored if > 5 pixels in a 5x5 center window have zero weight.
        weight_np = np.array(weight)
        nx, ny = weight_np.shape
        can_measure_flux = np.ones(Ns, dtype=bool)
        for i in range(Ns):
            s = candidates[i]
            s_imin = int(np.floor(s['xc_init'] + 0.5)) - psf.h_size_x
            s_imax = int(np.floor(s['xc_init'] + 0.5)) + psf.h_size_x + 1
            s_jmin = int(np.floor(s['yc_init'] + 0.5)) - psf.h_size_y
            s_jmax = int(np.floor(s['yc_init'] + 0.5)) + psf.h_size_y + 1
            
            c_x, c_y = int(np.floor(c_xc[i] + 0.5)), int(np.floor(c_yc[i] + 0.5))
            i_start, i_end = max(s_imin, c_x - 2), min(s_imax, c_x + 3)
            j_start, j_end = max(s_jmin, c_y - 2), min(s_jmax, c_y + 3)
            
            i_start, i_end = max(0, i_start), min(nx, i_end)
            j_start, j_end = max(0, j_start), min(ny, j_end)
            
            window = weight_np[i_start:i_end, j_start:j_end]
            if np.sum(window == 0) > 5:
                can_measure_flux[i] = False
        
        # C++ Logic: Chi2 masking per wavelength (compare_spots_chi2_and_mask)
        status = np.ones(Ns, dtype=int)
        status[~can_measure_flux] = 0
        unique_waves = np.unique(waves)
        for uw in unique_waves:
            mask = np.abs(waves - uw) < 1.0
            if np.sum(mask) >= 2:
                w_chi2 = chi2s[mask]
                m_chi2 = np.mean(w_chi2); std_chi2 = np.std(w_chi2)
                if std_chi2 > 0:
                    outliers = mask & (chi2s > (m_chi2 + 4.0 * std_chi2))
                    status[outliers] = 0
        
        # First selection pass (S/N, Image Bounds, Distance)
        first_pass_count = 0
        for i in range(Ns):
            if status[i] == 0: continue
            # C++ parity: reject spots with non‑positive eflux (spot->eflux <= 0).
            # In the C++ code eflux is computed after a dedicated flux fit.  Our JAX helper now
            # returns the weight‑sum A = Σ w·p², which is proportional to 1/eflux².  A non‑positive
            # A corresponds exactly to eflux ≤ 0.  We therefore mask any spot with A <= 0.
            if As[i] <= 0:
                status[i] = 0
                continue
            # After the eflux check we apply the usual S/N cut.
            if snrs[i] < sn_threshold:
                status[i] = 0
                continue
            # Distance check: C++ checks distance against all input_spots for the SAME fiber
            if min_dist_angstrom > 0:
                dist = 1000.0
                for j in range(Ns):
                    if i == j: continue
                    if all_candidates[j]['fiber'] != all_candidates[i]['fiber']: continue
                    # C++ does NOT check if input_spots[j] is masked here, it uses all input_spots
                    dist = min(dist, abs(waves[j] - waves[i]))
                if dist < min_dist_angstrom:
                    status[i] = 0
            if status[i] == 1: first_pass_count += 1
        
        # CHECKPOINT 1: Post-S/N & Bounds
        if hasattr(psf, 'output_psf_path') and psf.output_psf_path:
            cp1_path = psf.output_psf_path.replace('.fits', '.py_cp1.txt')
            with open(cp1_path, 'w') as f:
                for i in range(Ns):
                    if status[i] == 1:
                        s = candidates[i]
                        f.write(f"{s['fiber']},{s['wave']:.15f},{s['xc_init']:.15f},{s['yc_init']:.15f}\n")
            print(f"  Written {first_pass_count} spots to {cp1_path}", flush=True)
        
        # Second selection pass (Complex coverage algorithm)
        if max_number_of_lines > 0:
            selected_waves = {} # waveid -> bool
            nspots_per_wave = {} # waveid -> count
            snr_per_wave = {} # waveid -> avg_snr
            
            def get_id(w): return int(round(w * 10))
            
            for i in range(Ns):
                if status[i] == 0: continue
                wid = get_id(waves[i])
                if wid not in nspots_per_wave:
                    selected_waves[wid] = 1
                    nspots_per_wave[wid] = 1
                    snr_per_wave[wid] = snrs[i]
                else:
                    snr_per_wave[wid] = (snr_per_wave[wid] * nspots_per_wave[wid] + snrs[i]) / (nspots_per_wave[wid] + 1)
                    nspots_per_wave[wid] += 1
            
            sorted_ids = sorted(nspots_per_wave.keys())
            if not sorted_ids:
                first_pass_count = sum(status)
                # Ensure we return something if selection fails early
                # (Though usually we'd just skip to final)
                # But let's make sure we handle empty sorted_ids
                pass 
            else:
                max_fibers = max(nspots_per_wave.values())
                
                begin_id = 0
                end_id = 0
                has_found_max = False
                for wid in sorted_ids:
                    if has_found_max:
                        begin_id = wid
                        break
                    if nspots_per_wave[wid] == max_fibers:
                        has_found_max = True
                for wid in reversed(sorted_ids):
                    if nspots_per_wave[wid] == max_fibers:
                        end_id = wid - 1
                        break
                
            # --- Pass 2a: Remove low-SNR lines while maintaining gap constraints ---
            if sorted_ids:
                while True:
                    # C++: count lines in at least 60% of fibers
                    num_lines = sum(1 for wid in sorted_ids if selected_waves.get(wid, 0) == 1 and nspots_per_wave[wid] >= max_fibers * 0.6)
                    if num_lines <= max_number_of_lines: 
                        break
                    
                    # C++: sort selected flux by increasing order
                    # wave_vs_snr = {snr: waveid}
                    wave_vs_snr = sorted([(snr_per_wave[wid], wid) for wid in sorted_ids if selected_waves.get(wid, 0) == 1])
                    
                    waveid_to_remove = None
                    for snr, wid in wave_vs_snr:
                        if wid <= begin_id or wid >= end_id:
                            continue
                        
                        # Find previous and next selected waveids
                        previous_selected_waveid = begin_id
                        next_selected_waveid = end_id
                        
                        # C++: search through all selected to find tightest bounds
                        for swid in sorted_ids:
                            if selected_waves.get(swid, 0) == 1:
                                if swid < wid and swid > previous_selected_waveid:
                                    previous_selected_waveid = swid
                                if swid > wid and swid < next_selected_waveid:
                                    next_selected_waveid = swid
                        
                        dwave = (next_selected_waveid - previous_selected_waveid) / 10.0
                        if dwave > 300.0:
                            # C++ Logic: if(dwave > max_dwave) continue;
                            continue
                        
                        waveid_to_remove = wid
                        break
                    
                    if waveid_to_remove is None: 
                        break
                        
                    selected_waves[waveid_to_remove] = 0
            
            # --- Pass 2b: Bring back lines that are very close to selected lines ---
            while True:
                brought_back = False
                # C++: Nested loop over all sorted_ids (including those not yet selected)
                # The C++ loop is: for(it = selected.begin... for(jt = selected.begin...))
                # 'selected' is a map of waveid -> bool.
                for wid_s in sorted_ids:
                    if selected_waves.get(wid_s, 0) == 0: continue
                    for wid_j in sorted_ids:
                        if selected_waves.get(wid_j, 0) == 1: continue
                        # C++: if(fabs( it->first/10. - jt->first/10.) < min_dwave)
                        if abs(wid_s - wid_j) / 10.0 < 5.0:
                            selected_waves[wid_j] = 1
                            brought_back = True
                if not brought_back: break
                
                # Apply changes to status array after each full iteration of the "bring back" loop
                for i in range(Ns):
                    if selected_waves.get(get_id(waves[i]), 0) == 0:
                        status[i] = 0

        
        # CHECKPOINT 2: Final Selection
        if hasattr(psf, 'output_psf_path') and psf.output_psf_path:
            cp2_path = psf.output_psf_path.replace('.fits', '.py_cp2.txt')
            with open(cp2_path, 'w') as f:
                for i in range(Ns):
                    if status[i] == 1:
                        s = candidates[i]
                        f.write(f"{s['fiber']},{s['wave']:.15f},{s['xc_init']:.15f},{s['yc_init']:.15f}\n")
            print(f"  Written {sum(status)} spots to {cp2_path}", flush=True)

        selected = []
        for i in range(Ns):
            if status[i] == 1:
                s = candidates[i]; s.update({'flux': float(fluxes[i]), 'snr': float(snrs[i]), 'chi2': float(chi2s[i]), 'stamp_imin': int(np.floor(s['xc_init'] + 0.5)) - psf.h_size_x, 'stamp_imax': int(np.floor(s['xc_init'] + 0.5)) + psf.h_size_x + 1, 'stamp_jmin': int(np.floor(s['yc_init'] + 0.5)) - psf.h_size_y, 'stamp_jmax': int(np.floor(s['yc_init'] + 0.5)) + psf.h_size_y + 1})
                selected.append(s)
        
        # Write spots to sidecar file if requested (via psf.output_psf_path)
        if hasattr(psf, 'output_psf_path') and psf.output_psf_path:
            spots_path = psf.output_psf_path.replace('.fits', '.pyspots.txt')
            with open(spots_path, 'w') as f:
                for s in selected:
                    f.write(f"{s['fiber']},{s['wave']:.15f},{s['xc_init']:.15f},{s['yc_init']:.15f}\n")
            print(f"  Written {len(selected)} spots to {spots_path}", flush=True)

        print(f"  Spot selection took {time.time() - t0:.2f}s ({len(selected)} spots)", flush=True)
        return selected
    return candidates


def _fit_one_spot_jax(image, weight, xc, yc, gh, degree, hsize_x, hsize_y):
    import jax.numpy as jnp
    from jax import jit
    nx, ny = image.shape
    ix_rel, iy_rel = jnp.meshgrid(jnp.arange(2*hsize_x+1), jnp.arange(2*hsize_y+1), indexing='ij')
    dx = ix_rel.flatten() - hsize_x; dy = iy_rel.flatten() - hsize_y
    
    def fit_loop(flux, x_shift, y_shift):
        # Simple 3-parameter fit (flux, dx, dy) for selection purposes
        # This mimics C++ FitOneSpot's initial phase
        im = jnp.floor(xc + 0.5).astype(int); jm = jnp.floor(yc + 0.5).astype(int)
        gix = im + dx + x_shift; giy = jm + dy + y_shift
        valid = (gix >= 0) & (gix < nx) & (giy >= 0) & (giy < ny)
        p_val = GaussHermitePSF.single_pix_value_jnp(xc + x_shift, yc + y_shift, gix, giy, gh, degree)
        p_val = jnp.where(valid, p_val, 0.0)
        d_val = image[jnp.clip(gix, 0, nx-1), jnp.clip(giy, 0, ny-1)]
        w_val = weight[jnp.clip(gix, 0, nx-1), jnp.clip(giy, 0, ny-1)]
        w_val = jnp.where(valid, w_val, 0.0)
        
        A = jnp.sum(w_val * p_val**2)
        B = jnp.sum(w_val * d_val * p_val)
        flux_est = jnp.where(A > 0, B/A, 0.0)
        return flux_est, A

    return jit(fit_loop)(0.0, 0.0, 0.0)

def _get_spot_stats_jax(image, weight, cand_xc, cand_yc, gh_params, degree, hsize_x, hsize_y):
    import jax.numpy as jnp
    from jax import vmap, jit, grad
    nx, ny = image.shape; area = (2*hsize_x+1)*(2*hsize_y+1)
    ix_rel, iy_rel = jnp.meshgrid(jnp.arange(2*hsize_x+1), jnp.arange(2*hsize_y+1), indexing='ij')
    dx = ix_rel.flatten() - hsize_x; dy = iy_rel.flatten() - hsize_y

    def spot_objective(params, xc_init, yc_init, gh):
        flux, dx_s, dy_s = params
        xc = xc_init + dx_s; yc = yc_init + dy_s
        im = jnp.floor(xc + 0.5).astype(int); jm = jnp.floor(yc + 0.5).astype(int)
        gix = im + dx; giy = jm + dy
        valid = (gix >= 0) & (gix < nx) & (giy >= 0) & (giy < ny)
        p_val = GaussHermitePSF.single_pix_value_jnp(xc, yc, gix, giy, gh, degree)
        p_val = jnp.where(valid, p_val, 0.0)
        d_val = image[jnp.clip(gix, 0, nx-1), jnp.clip(giy, 0, ny-1)]
        w_val = weight[jnp.clip(gix, 0, nx-1), jnp.clip(giy, 0, ny-1)]
        w_val = jnp.where(valid, w_val, 0.0)
        return jnp.sum(w_val * (d_val - flux * p_val)**2)

    def fit_spot(xc, yc, gh):
        # Simple flux‑only fit matching C++ FitIndividualSpotFluxes (no centroid shift)
        im = jnp.floor(xc + 0.5).astype(int); jm = jnp.floor(yc + 0.5).astype(int)
        gix = im + dx; giy = jm + dy
        valid = (gix >= 0) & (gix < nx) & (giy >= 0) & (giy < ny)
        p_val = GaussHermitePSF.single_pix_value_jnp(xc, yc, gix, giy, gh, degree)
        p_val = jnp.where(valid, p_val, 0.0)
        d_val = image[jnp.clip(gix, 0, nx-1), jnp.clip(giy, 0, ny-1)]
        w_val = weight[jnp.clip(gix, 0, nx-1), jnp.clip(giy, 0, ny-1)]
        w_val = jnp.where(valid, w_val, 0.0)
        A = jnp.sum(w_val * p_val**2)
        B = jnp.sum(w_val * d_val * p_val)
        flux = jnp.where(A > 0, B/A, 0.0)
        snr = jnp.where(A > 0, flux * jnp.sqrt(A), -1.0)
        chi2 = jnp.sum(w_val * (d_val - flux * p_val)**2)
        return flux, snr, chi2, A

    # Return flux, snr, chi2, and the weight‑sum A for each spot (required for eflux parity)
    return jit(vmap(fit_spot))(cand_xc, cand_yc, gh_params)

class PSF_Fitter:
    def __init__(self, psf):
        self.psf = psf; self.chi2_precision = 10.0
    def fit(self, image, weight, spots, bundle_id, fit_type='full', max_iter=20):
        import jax.numpy as jnp
        print(f"Starting HIGH-PERFORMANCE OPTIMIZED fit for bundle {bundle_id}...")
        fmin, fmax = spots[0]['fiber'], spots[-1]['fiber']
        weight = apply_dead_column_mask(self.psf, fmin, fmax, weight)
        xpix, ypix, pix_idx = get_bundle_footprint(self.psf, spots, fmin, fmax, weight)
        Np = len(xpix); area = (2*self.psf.h_size_x+1)*(2*self.psf.h_size_y+1); Ns = len(spots)
        t0 = time.time(); sx, sy = np.zeros((Ns, area)), np.zeros((Ns, area)); idx_g = np.full((Ns, area), Np, dtype=np.int32)
        nx, ny = image.shape; idx_map = np.full((nx, ny), -1, dtype=np.int32); idx_map[xpix, ypix] = np.arange(Np)
        for s_i, s in enumerate(spots):
            ix, iy = np.meshgrid(np.arange(s['stamp_imin'], s['stamp_imax']), np.arange(s['stamp_jmin'], s['stamp_jmax']), indexing='ij')
            sx[s_i], sy[s_i] = ix.flatten(), iy.flatten(); st_idx = idx_map[ix.astype(int), iy.astype(int)].flatten(); mask = st_idx >= 0; idx_g[s_i, mask] = st_idx[mask]
        print(f"  Stamp indexing took {time.time() - t0:.2f}s", flush=True)
        rows_u = np.unique(ypix); row_m = {j: i for i, j in enumerate(rows_u)}
        tx_j, tw_j = np.zeros((25, len(rows_u))), np.zeros((25, len(rows_u)))
        for f_i in range(25):
            fib = fmin + f_i; w_v = self.psf.fiber_traces[fib]['Y_vs_W'].invert(rows_u.astype(float)); tw_j[f_i], tx_j[f_i] = w_v, self.psf.x_ccd(fib, w_v)
        ix_r = np.array([row_m[j] for j in ypix]); tx_g, tw_g = jnp.array(tx_j[:, ix_r]), jnp.array(tw_j[:, ix_r])
        flux = jnp.array([s['flux'] for s in spots]); xc_init, yc_init = jnp.array([s['xc_init'] for s in spots]), jnp.array([s['yc_init'] for s in spots])
        monomials = get_bundle_monomials_jnp(self.psf, bundle_id, spots)
        gh_deg = self.psf.gh_psf.degree; n_gh = (gh_deg + 1) * (gh_deg + 1) - 1
        pc = jnp.zeros((n_gh + 2, monomials.shape[1])).at[0, 0].set(1.1).at[1, 0].set(1.1)
        tc = jnp.zeros((2, monomials.shape[1])); Ncont = 4; cc = jnp.zeros(Ncont) 
        img_d, w_d = jnp.array(image[xpix, ypix]), jnp.array(weight[xpix, ypix])
        wmin_c, wmax_c = float(self.psf.fiber_traces[fmin]['X_vs_W'].xmin), float(self.psf.fiber_traces[fmin]['X_vs_W'].xmax); old_chi2 = 1e30
        
        best_chi2 = 1e30
        best_tc = tc.copy()
        best_pc = pc.copy()
        best_cc = cc.copy()
        best_flux = flux.copy()
        
        sx_g, sy_g, idx_gg = jnp.array(sx), jnp.array(sy), jnp.array(idx_g)
        for i in range(max_iter):
            chi2, A, B = _accumulate_bundle_jax_jit(flux, pc, tc, cc, xc_init, yc_init, monomials, jnp.array(xpix), jnp.array(ypix), sx_g, sy_g, idx_gg, gh_deg, tx_g, tw_g, wmin_c, wmax_c, img_d, w_d, self.psf.gain, self.psf.psf_error, 1.0)
            
            if chi2 < best_chi2:
                best_chi2 = chi2
                best_tc = tc.copy()
                best_pc = pc.copy()
                best_cc = cc.copy()
                best_flux = flux.copy()
                
            mode = 'flux' if i < 2 else 'trace' if i < 5 else 'full'
            print(f"Iter {i}: chi2 = {float(chi2):.4f} [Mode: {mode}]", flush=True)
            Npoly = monomials.shape[1]; Ns_l = len(flux); n_psf_tot = (n_gh + 2) * Npoly
            if mode == 'flux': idx = jnp.concatenate([jnp.arange(Ns_l), jnp.arange(A.shape[0]-Ncont, A.shape[0])])
            elif mode == 'trace': idx = jnp.concatenate([jnp.arange(Ns_l), jnp.arange(Ns_l + n_psf_tot, Ns_l + n_psf_tot + 2*Npoly), jnp.arange(A.shape[0]-Ncont, A.shape[0])])
            else: idx = jnp.arange(A.shape[0])
            A_sub, B_sub = A[jnp.ix_(idx, idx)], B[idx]; diag = jnp.diag(A_sub); S = jnp.sqrt(diag); S = jnp.where(S < 1e-12, 1.0, S)
            A_reg = (A_sub / jnp.outer(S, S)) + 1e-8 * jnp.eye(A_sub.shape[0])
            try: ds = jnp.linalg.solve(A_reg, B_sub / S); d_p = jnp.zeros(A.shape[0]).at[idx].set(ds / S)
            except: d_p = jnp.zeros(A.shape[0])
            best_alpha, best_chi2 = 0.0, float(chi2)
            if jnp.any(d_p != 0):
                for alpha in [0.2, 0.5, 1.0]:
                    f_try = jnp.maximum(flux + alpha * d_p[:Ns_l], 0.0); p_try = pc + alpha * d_p[Ns_l : Ns_l + n_psf_tot].reshape(n_gh + 2, Npoly); t_try = tc + alpha * d_p[Ns_l + n_psf_tot : Ns_l + n_psf_tot + 2*Npoly].reshape(2, Npoly); c_try = cc + alpha * d_p[-Ncont:]; c2 = _predict_bundle_jax_jit(f_try, p_try, t_try, c_try, xc_init, yc_init, monomials, jnp.array(xpix), jnp.array(ypix), sx_g, sy_g, idx_gg, gh_deg, tx_g, tw_g, wmin_c, wmax_c, img_d, w_d)
                    if c2 < best_chi2: best_alpha, best_chi2 = alpha, c2
            if best_alpha == 0 and i > 5: break
            if best_alpha == 0: best_alpha = 0.1
            flux = jnp.maximum(flux + best_alpha * d_p[:Ns_l], 0.0); pc = pc + best_alpha * d_p[Ns_l : Ns_l + n_psf_tot].reshape(n_gh + 2, Npoly); tc = tc + best_alpha * d_p[Ns_l + n_psf_tot : Ns_l + n_psf_tot + 2*Npoly].reshape(2, Npoly); cc = cc + best_alpha * d_p[-Ncont:]
            
            # --- Iterative Snapping: Update xc_init/yc_init to the current model prediction ---
            # We use a staged approach to prevent oscillation.
            # Flux mode: No snapping. Trace mode: Full snapping. Full mode: Increased frequency.
            if mode == 'trace':
                # Removed iterative snapping to align with C++ logic.
                # xc_init and yc_init must remain constant anchors.
                pass
            
            # REMOVED: All iterative snapping in 'full' mode. 
            # We keep xc_init fixed from the end of trace mode.
            # tc will now naturally accumulate the total shift from this anchor.
            
            if mode == 'full' and jnp.abs(old_chi2 - chi2) < self.chi2_precision: break
            old_chi2 = chi2
        
        # --- C++ Parity: Snap centroids to the final optimized model ---
        # Use the best coefficients found during the optimization process
        import jax.numpy as jnp
        dx_final = jnp.dot(monomials, best_tc[0])
        dy_final = jnp.dot(monomials, best_tc[1])
        
        # DEBUG: Check if the shifts are actually non-zero
        print(f"  DEBUG: dx_final mean={np.mean(np.abs(dx_final)):.6f}, dy_final mean={np.mean(np.abs(dy_final)):.6f}", flush=True)
        
        # The final position is the anchor (xc_init) plus the optimized shift
        # xc_init remained constant throughout the fit, mirroring C++ logic
        xc_final = np.array(xc_init + dx_final)
        yc_final = np.array(yc_init + dy_final)
        
        return float(best_chi2), np.array(best_pc), np.array(best_tc), np.array(best_cc), np.array(best_flux), xc_final, yc_final
