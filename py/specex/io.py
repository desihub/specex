import os
import numpy as np
import fitsio
from fitsio import FITS, FITSHDR
from datetime import datetime
from .psf import PSF, PSF_Params
from .math import SparseLegendre2DPol, Legendre1DPol, legendre_pol_jnp
from .fitter import get_sparse_nz

def meta2header(meta):
    """Convert a plain Python dict of FITS header metadata into the C++ extension's MapStringString header type (specex._libspecex), quoting strings and stripping trailing '.0' from float-looking values the way C++'s own header formatting does.

    Args:
        meta (dict): header key -> value (bool/str/other), typically a
            fitsio-read header converted to a dict.

    Returns:
        specex._libspecex.MapStringString: the C++-side header object.

    Status: ACTIVE (C++-wrapper path) -- only called by read_preproc_cpp,
    itself only used by specex.run_specex(). Real production code for
    the C++ codepath, not unused legacy (see run_specex()'s own
    docstring); just a different entry point than the GPU-native path.
    """
    import specex._libspecex as spx
    header = spx.MapStringString()
    for key in meta:
        mkey = meta[key]
        if type(mkey) == bool:
            header[key]='T' if mkey else 'F'
        elif type(mkey) == str: 
            mstr = mkey
            if len(mstr) < 8: mstr = mstr.ljust(8, ' ')
            header[key]="\'"+mstr+"\'"
        else:
            mstr = str(mkey)
            if len(mstr) > 1 and mstr[-2:]=='.0': mstr=mstr[:-1]
            header[key]=mstr
    return header

def load_python_psf(filename, opts):
    """Read an input (shifted) PSF FITS file into a PSF object: per-fiber XTRACE/YTRACE Legendre trace polynomials, and (if a PSF extension is present) per-bundle Gauss-Hermite shape parameter models.

    Args:
        filename (str): path to the input PSF FITS file.
        opts: unused (accepted for interface symmetry with C++-wrapper-style
            callers).

    Returns:
        psf.PSF: populated with fiber_traces (all fibers in [FIBERMIN,
        FIBERMAX]) and, if a PSF extension exists, params_of_bundles (one
        psf.PSF_Params per bundle found via the header's `B<NN>NDATA` keys).

    Status: ACTIVE (production default path).
    """
    f = fitsio.FITS(filename)
    psf = PSF()
    xt_hdr = f['XTRACE'].read_header()
    psf.fiber_min = xt_hdr['FIBERMIN']
    psf.fiber_max = xt_hdr['FIBERMAX']
    xtrace = f['XTRACE'].read().astype(np.float64)
    ytrace = f['YTRACE'].read().astype(np.float64)
    for fib in range(psf.fiber_min, psf.fiber_max + 1):
        idx = fib - psf.fiber_min
        psf.fiber_traces[fib] = {
            'X_vs_W': Legendre1DPol(deg=xtrace.shape[1]-1, xmin=xt_hdr['WAVEMIN'], xmax=xt_hdr['WAVEMAX'], coeff=xtrace[idx]),
            'Y_vs_W': Legendre1DPol(deg=ytrace.shape[1]-1, xmin=xt_hdr['WAVEMIN'], xmax=xt_hdr['WAVEMAX'], coeff=ytrace[idx])
        }
    if 'PSF' in f:
        table = f['PSF'].read()
        hdr = f['PSF'].read_header()
        psf.gh_psf.degree = hdr['GHDEGX']
        psf.h_size_x = hdr['HSIZEX']; psf.h_size_y = hdr['HSIZEY']
        psf.gain = hdr['GAIN']; psf.readout_noise = hdr['READNOIS']
        bundle_ids = [int(key[1:3]) for key in hdr.keys() if key.startswith('B') and key.endswith('NDATA')]
        coeffs_all = table['COEFF'].astype(np.float64)
        param_names = [p.strip() for p in table['PARAM']]
        for bid in bundle_ids:
            b_fmin, b_fmax = bid * 25, (bid + 1) * 25 - 1
            bundle = PSF_Params(bid, b_fmin, b_fmax)
            bundle.param_names = param_names; bundle.param_models = {}
            for i, name in enumerate(param_names):
                bundle.param_models[name] = [
                    Legendre1DPol(deg=coeffs_all.shape[2]-1, xmin=hdr['WAVEMIN'], xmax=hdr['WAVEMAX'], coeff=coeffs_all[i][fib])
                    for fib in range(500)
                ]
            psf.params_of_bundles[bid] = bundle
    return psf

def write_python_psf(filename, bundle_results, input_template):
    """Merge all fitted bundles' results into a single output PSF FITS file: apply each bundle's trace correction to the input template's XTRACE/YTRACE, write its fitted Gauss-Hermite coefficients into the PSF table, restore never-fit fibers (broken/masked-amp) to their input values with STATUS=-1, run an inline trace-crossing QA pass (STATUS=4), and copy through any other input extensions unchanged.

    Args:
        filename (str): output PSF FITS file path (overwritten if it exists).
        bundle_results (dict[int, dict]): per-bundle result dicts as returned
            by specex.fit_bundle_task (or `{'skip_bundle': True, ...}` for a
            fully-excluded bundle).
        input_template (str): input (shifted) PSF FITS file path -- supplies
            the starting-guess trace/PSF-shape values every correction is
            applied on top of, and every extension not itself recomputed here.

    Returns:
        None. Writes `filename` as a side effect.

    Status: ACTIVE (production default path).
    """
    import fitsio
    fin = fitsio.FITS(input_template)
    xtrace_out = fin['XTRACE'].read().astype(np.float64)
    ytrace_out = fin['YTRACE'].read().astype(np.float64)
    xtrace_input_orig = xtrace_out.copy()
    ytrace_input_orig = ytrace_out.copy()
    psf_table = fin['PSF'].read()
    psf_hdr = fin['PSF'].read_header()
    param_names = [p.strip() for p in psf_table['PARAM']]
    name_to_idx = {name: i for i, name in enumerate(param_names)}
    status_idx = name_to_idx.get('STATUS')
    # Pristine copy of every PSF-table param (GH shape, tails, STATUS, ...)
    # for the never-fit-fiber restoration below -- mirrors
    # xtrace_input_orig/ytrace_input_orig, needed because the per-bundle
    # loop below zeroes psf_table['COEFF'] in place before we get a chance
    # to restore a specific fiber's original value.
    psf_coeff_input_orig = psf_table['COEFF'].copy()
    xdeg_b = 1
    for bid, res in bundle_results.items():
        fmin, fmax = bid * 25, (bid + 1) * 25 - 1
        # A bundle where EVERY fiber was excluded (fully inside a masked
        # amp, or all explicitly broken) never got fit at all -- skip the
        # normal per-bundle correction write entirely; the never-fit-fiber
        # restoration below (which runs unconditionally) handles it.
        if not res.get('skip_bundle'):
            pc = res['psf_coeffs']; tc = res['trace_coeffs']
            # wdeg varies by band (3 for z, 1 otherwise -- see fit_ccd_native's
            # auto-detection); carried through bundle_results since the writer
            # has no other way to know what basis pc/tc were fit in. trace_wdeg
            # is a separate, independently-sized basis for tc (defaults to
            # wdeg when not set -- see docs/python-port/porting-notes.md's r2@20250109
            # investigation for why pc and tc can now differ here).
            wdeg_b = res.get('wdeg', 3); nz_b = get_sparse_nz(xdeg_b, wdeg_b)
            rf = 2 * (np.arange(fmin, fmax + 1) - fmin) / (fmax - fmin) - 1
            poly_f = np.stack([legendre_pol_jnp(i, rf) for i in range(xdeg_b + 1)], axis=0)
            trace_per_fiber_deg = res.get('trace_per_fiber_deg')
            if trace_per_fiber_deg is not None:
                # Stage 1 of the full per-fiber redesign (see docs/python-port/porting-notes.md):
                # tc[0]/tc[1] are (n_fibers*(deg+1),) block-diagonal-by-fiber
                # coefficient vectors, not a shared basis -- reshape to
                # (n_fibers, deg+1) and add each fiber's own coefficients
                # straight into its own XTRACE/YTRACE row, no fiber-position
                # broadcast basis (poly_f) involved at all.
                n_fibers = fmax - fmin + 1; n_coefs = trace_per_fiber_deg + 1
                tc_x_pf = tc[0].reshape(n_fibers, n_coefs); tc_y_pf = tc[1].reshape(n_fibers, n_coefs)
                n_write = min(n_coefs, xtrace_out.shape[1])
                xtrace_out[fmin:fmax+1, :n_write] += tc_x_pf[:, :n_write]
                ytrace_out[fmin:fmax+1, :n_write] += tc_y_pf[:, :n_write]
            else:
                trace_wdeg_b = res.get('trace_wdeg', wdeg_b); nz_trace_b = get_sparse_nz(xdeg_b, trace_wdeg_b)
                for k_nz, k_lin in enumerate(nz_trace_b):
                    i_p, j_p = k_lin % 2, k_lin // 2
                    if j_p < xtrace_out.shape[1]:
                        xtrace_out[fmin:fmax+1, j_p] += tc[0, k_nz] * poly_f[i_p]
                        ytrace_out[fmin:fmax+1, j_p] += tc[1, k_nz] * poly_f[i_p]
            # Zero the trace row for any fiber with zero selected spots, matching
            # C++'s own output for such a fiber (specex_psf_proc.cc:49,58 --
            # Trace::resize(0) empties that fiber's coeff array, so it's never
            # copied into the zero-initialized output buffer). Without this,
            # Python instead writes a smoothed/interpolated position inherited
            # from neighboring fibers -- plausible-looking but never actually
            # constrained by any real data for that fiber, unlike C++'s
            # unambiguous all-zero "don't trust this" signal. Overridden last so
            # it wins regardless of which branch above ran.
            for fib in res.get('zero_spot_fibers', []):
                if fmin <= fib <= fmax:
                    xtrace_out[fib, :] = 0.0
                    ytrace_out[fib, :] = 0.0
            for row in range(len(param_names)):
                psf_table['COEFF'][row, fmin:fmax+1, :] = 0.0
                if param_names[row] == 'GH-0-0': psf_table['COEFF'][row, fmin:fmax+1, 0] = 1.0
                # TAILXSCA/TAILYSCA/TAILCORE are never fit (tail amplitude is
                # always 0 here), but C++ still writes 1.0 for these three,
                # not 0.0. Matching that convention is functionally a no-op
                # (TAILAMP=0 zeroes the whole tail term regardless), but a 0
                # here makes downstream tools that evaluate the tail formula
                # unconditionally (e.g. specter's gausshermite.py, used by
                # Julien's plot_psf_comparison_using_specter.py) divide by
                # zero -- see docs/python-port/porting-notes.md 2026-08-23.
                if param_names[row] in ('TAILXSCA', 'TAILYSCA', 'TAILCORE'): psf_table['COEFF'][row, fmin:fmax+1, 0] = 1.0
            # Name mapping for GH terms. Must match the fitter's pc row order
            # exactly: the full (deg+1)^2-1 grid excluding only (0,0), same
            # convention as PSF.canonical_param_names() and the inner loop of
            # _accumulate_bundle_jax (fitter.py). No i+j<=deg triangular filter -
            # that filter (present here previously) desynchronized the enumerate
            # index from the pc rows starting at GH-6-1, scrambling all shape
            # coefficients written after that point.
            gh_deg = psf_hdr['GHDEGX']
            param_mapping = ['GHSIGX', 'GHSIGY']
            for j_gh in range(gh_deg + 1):
                for i_gh in range(gh_deg + 1):
                    if i_gh == 0 and j_gh == 0: continue
                    param_mapping.append(f'GH-{i_gh}-{j_gh}')
            if len(param_mapping) != pc.shape[0]:
                raise ValueError(f"PSF param mapping length {len(param_mapping)} != fitted coeff rows {pc.shape[0]}")

            for i_par, pname in enumerate(param_mapping):
                idx = name_to_idx.get(pname)
                if idx is not None:
                    for k_nz, k_lin in enumerate(nz_b):
                        i_p, j_p = k_lin % 2, k_lin // 2
                        if j_p < psf_table['COEFF'].shape[2]:
                            psf_table['COEFF'][idx, fmin:fmax+1, j_p] += pc[i_par, k_nz] * poly_f[i_p]
            psf_hdr[f'B{bid:02d}RCHI2'] = res['chi2'] / (120000.0)

        # Never-fit fibers -- explicitly-listed --broken-fibers, dynamically-
        # detected masked-amp fibers, or an entire skip_bundle bundle --
        # are excluded from the fit entirely and must be left completely
        # untouched at the input template's own values (trace AND PSF
        # shape), with STATUS flagged -1. Confirmed against real C++
        # production output (bit-for-bit identical to the input PSF, both
        # trace and GH-shape rows, STATUS=-1) for both the explicit-broken
        # case (z8@20260401 fibers 473/474, z3@20260401 fiber 368,
        # b8@20221121 fibers 348/473/474) and the masked-amp case
        # (r8@20211028 fibers 0-254) -- see docs/python-port/porting-notes.md's 2026-08-14
        # writeup. Runs after the normal per-bundle write above (when it
        # ran at all) so it always wins for these specific fibers,
        # regardless of what the broadcast correction wrote elsewhere in
        # the bundle -- the broadcast is a no-op to undo since
        # *_input_orig are the pristine pre-any-bundle values.
        never_fit_fibers = set(res.get('explicitly_broken_fibers', [])) | set(res.get('masked_amp_fibers', []))
        for fib in never_fit_fibers:
            if fmin <= fib <= fmax:
                xtrace_out[fib, :] = xtrace_input_orig[fib, :]
                ytrace_out[fib, :] = ytrace_input_orig[fib, :]
                psf_table['COEFF'][:, fib, :] = psf_coeff_input_orig[:, fib, :]
                if status_idx is not None:
                    psf_table['COEFF'][status_idx, fib, :] = 0.0
                    psf_table['COEFF'][status_idx, fib, 0] = -1.0

    # Overlapping-trace QA (adapted from specex#91's trace_psf_qa, py/specex/
    # qa.py -- present in this repo but dead code, imported only by the
    # legacy C++-wrapper path and never actually called). Detects
    # neighboring fiber pairs whose FITTED traces cross, and flags the pair
    # PLUS their immediate neighbors (one fiber above/below) as STATUS=4
    # ("overlapping traces (QA)", matching specex#91's own repurposing of
    # that code) -- the neighbor-flagging is a deliberate difference from
    # specex#91's own version (which only flags the crossing pair itself)
    # and from this project's own prior behavior (flag the whole bundle,
    # or crash the whole camera) -- see docs/python-port/porting-notes.md's 2026-08-14
    # writeup. This is the NON-fatal-accuracy case: trace values are left
    # exactly as fitted, only STATUS changes. Never-fit fibers (STATUS=-1
    # already, from the block above) are excluded on both sides of the
    # comparison -- their trace is a pass-through input value, not a real
    # fit, so an apparent "crossing" against a real neighbor is a
    # comparison artifact, not a genuine overlap, and their STATUS=-1
    # must not be downgraded to a 4.
    if status_idx is not None:
        no_data_fibers = set()
        for res in bundle_results.values():
            no_data_fibers.update(res.get('explicitly_broken_fibers', []))
            no_data_fibers.update(res.get('masked_amp_fibers', []))
        xt_hdr = fin['XTRACE'].read_header()
        wmin, wmax = xt_hdr['WAVEMIN'], xt_hdr['WAVEMAX']
        ww = np.linspace(wmin, wmax, 200)
        n_fibers_total = xtrace_out.shape[0]
        x_cache = {}
        def _xval(fib):
            """Evaluate (and memoize) a fiber's fitted X-trace position over the shared 200-point wavelength grid, for the trace-crossing QA pass.

            Args:
                fib (int): absolute fiber index.

            Returns:
                np.ndarray: shape (200,), X position at each grid wavelength.

            Status: ACTIVE (production default path) -- internal helper closure of
            write_python_psf's QA pass.
            """
            if fib not in x_cache:
                x_cache[fib] = np.array(Legendre1DPol(deg=xtrace_out.shape[1]-1, xmin=wmin, xmax=wmax, coeff=xtrace_out[fib]).value(ww))
            return x_cache[fib]
        crossing_fibers = set()
        for f0 in range(n_fibers_total - 1):
            f1 = f0 + 1
            if f0 in no_data_fibers or f1 in no_data_fibers:
                continue
            if not np.all(_xval(f1) > _xval(f0)):
                for nb in (f0 - 1, f0, f1, f1 + 1):
                    if 0 <= nb < n_fibers_total and nb not in no_data_fibers:
                        crossing_fibers.add(nb)
        if crossing_fibers:
            print(f"QA: {len(crossing_fibers)} fiber(s) flagged STATUS=4 for overlapping traces "
                  f"(crossing pair + immediate neighbors): {sorted(crossing_fibers)}", flush=True)
            for fib in crossing_fibers:
                psf_table['COEFF'][status_idx, fib, :] = 0.0
                psf_table['COEFF'][status_idx, fib, 0] = 4.0

    if os.path.exists(filename): os.remove(filename)
    fout = fitsio.FITS(filename, 'rw')
    fout.write(xtrace_out, header=fin['XTRACE'].read_header(), extname='XTRACE')
    fout.write(ytrace_out, header=fin['YTRACE'].read_header(), extname='YTRACE')
    fout.write(psf_table, header=psf_hdr, extname='PSF')

    # Pass through any other extensions from the input template unchanged
    # (e.g. EXTOFF/INTOFF wavelength-offset tables). Real desi_compute_psf
    # merges never modify these -- merge_psf() only touches XTRACE/YTRACE/
    # PSF, so anything else in the input PSF file is inherited byte-for-byte
    # from an earlier upstream calibration step (desi_compute_trace_shifts),
    # not recomputed here.
    known = {'XTRACE', 'YTRACE', 'PSF'}
    for hdu in fin:
        extname = hdu.get_extname()
        if extname and extname not in known:
            fout.write(hdu.read(), header=hdu.read_header(), extname=extname)

    fout.close()

def write_psf(pyps, opts, pyio):
    """Write a C++-side fitted PSF (PyPSF object) out to a PSF FITS file, via the specex._libspecex extension's own trace/table accessors.

    Args:
        pyps: a specex._libspecex.PyPSF object (already fit).
        opts: a specex._libspecex.PyOptions object; `opts.output_fits_filename`
            is the output path (overwritten if it exists).
        pyio: a specex._libspecex.PyIO object.

    Returns:
        None. Writes `opts.output_fits_filename` as a side effect.

    Status: ACTIVE (C++-wrapper path) -- only called by
    specex.run_specex(). Real production code for the C++ codepath,
    not unused legacy; just a different entry point than the GPU-native
    path.
    """
    import specex._libspecex as spx
    pyio.load_psf(opts, pyps); spx.tablewrite_init(pyps)
    xtrace = spx.get_trace(pyps, 'x'); ytrace = spx.get_trace(pyps, 'y')
    xtrace = np.reshape(xtrace, (pyps.nfibers, pyps.trace_ncoeff))
    ytrace = np.reshape(ytrace, (pyps.nfibers, pyps.trace_ncoeff))
    table_col0 = spx.VectorString(); table_col1 = spx.VectorDouble(); table_col2 = spx.VectorInt(); table_col3 = spx.VectorInt()
    table_bundle_id = spx.VectorInt(); table_bundle_ndata = spx.VectorInt(); table_bundle_nparams = spx.VectorInt(); table_bundle_chi2pdf = spx.VectorDouble()
    spx.get_table(pyps, table_col0, table_col1, table_col2, table_col3, table_bundle_id, table_bundle_ndata, table_bundle_nparams, table_bundle_chi2pdf)
    col0 = np.zeros(pyps.table_nrows, dtype='S8')
    col1 = np.zeros((pyps.table_nrows, pyps.nfibers, pyps.ncoeff))
    col2 = np.zeros(pyps.table_nrows, dtype='i4')
    col3 = np.zeros(pyps.table_nrows, dtype='i4')
    i = 0
    for r in range(pyps.table_nrows):
        col0[r] = table_col0[r]
        for t2 in range(pyps.nfibers):
            for t1 in range(pyps.ncoeff):
                col1[r, t2, t1] = table_col1[i]; i += 1
        col2[r] = table_col2[r]; col3[r] = table_col3[r]
    data = np.zeros(pyps.table_nrows, dtype=[('PARAM', 'S8'), ('COEFF', 'f8', (pyps.nfibers, pyps.ncoeff)), ('LEGDEGX', 'i4'), ('LEGDEGW', 'i4')])
    data['PARAM'] = col0; data['COEFF'] = col1; data['LEGDEGX'] = col2; data['LEGDEGW'] = col3
    header = FITSHDR()
    header['GHDEGX'] = pyps.GHDEGX; header['GHDEGY'] = pyps.GHDEGY; header['MJD'] = pyps.mjd; header['PLATEID'] = pyps.plate_id; header['CAMERA'] = pyps.camera_id
    header['ARCEXP'] = pyps.arc_exposure_id; header['NPIX_X'] = pyps.NPIX_X; header['NPIX_Y'] = pyps.NPIX_Y; header['HSIZEX'] = pyps.hSizeX; header['HSIZEY'] = pyps.hSizeY
    header['FIBERMIN'] = pyps.FIBERMIN; header['FIBERMAX'] = pyps.FIBERMAX; header['WAVEMIN'] = pyps.table_WAVEMIN; header['WAVEMAX'] = pyps.table_WAVEMAX; header['LEGDEG'] = pyps.LEGDEG
    for i in range(len(table_bundle_id)):
        header[f'B{table_bundle_id[i]:02d}NDATA'] = table_bundle_ndata[i]; header[f'B{table_bundle_id[i]:02d}NPAR'] = table_bundle_nparams[i]; header[f'B{table_bundle_id[i]:02d}RCHI2'] = table_bundle_chi2pdf[i]
    if os.path.exists(opts.output_fits_filename): os.remove(opts.output_fits_filename)
    fout = fitsio.FITS(opts.output_fits_filename, 'rw')
    fout.write(xtrace, extname='XTRACE'); fout.write(ytrace, extname='YTRACE'); fout.write(data, header=header, extname='PSF')
    fout.close()

def read_image(filename):
    """Read a preproc FITS file's IMAGE/IVAR/MASK extensions and IMAGE header.

    Args:
        filename (str): path to the preproc FITS file.

    Returns:
        dict: {'image': np.ndarray (float64), 'ivar': np.ndarray (float64),
        'mask': np.ndarray (int32), 'meta': fitsio header of the IMAGE
        extension}.

    Status: ACTIVE (production default path) -- called by both read_preproc
    (the active Python pipeline) and read_preproc_cpp (LEGACY, C++-wrapper
    path).
    """
    f = fitsio.FITS(filename)
    image = f['IMAGE'].read().astype(np.float64); ivar = f['IVAR'].read().astype(np.float64); mask = f['MASK'].read().astype(np.int32); meta = f['IMAGE'].read_header()
    return {'image': image, 'ivar': ivar, 'mask': mask, 'meta': meta}

def read_preproc(opts):
    """Read a preproc arc image for the Python/JAX pipeline: read_image plus zeroing ivar at masked pixels and synthesizing a per-pixel readout-noise array from the header's RDNOISE.

    Args:
        opts: an object with an `arc_image_filename` (str) attribute (e.g.
            specex.fit_bundle_task's local Opts class).

    Returns:
        dict: read_image's dict, with 'ivar' zeroed where 'mask' != 0 and an
        added 'rdnoise' key (np.ndarray, same shape as 'image', filled with
        the header's RDNOISE value, default 0.0 if absent).

    Status: ACTIVE (production default path).
    """
    ddata = read_image(opts.arc_image_filename)
    ddata['ivar'][ddata['mask'] != 0] = 0.0
    rdnoise_meta = ddata['meta'].get('RDNOISE', 0.0)
    ddata['rdnoise'] = np.full_like(ddata['image'], float(rdnoise_meta))
    return ddata

def read_psf(opts, pyps):
    """Read an input PSF FITS file's trace and (if present) PSF-shape data directly into a C++-side PyPSF object, via the specex._libspecex extension.

    Args:
        opts: a specex._libspecex.PyOptions object (supplies
            `input_psf_filename`, `trace_deg_x`, `trace_deg_wave`).
        pyps: a specex._libspecex.PyPSF object, populated in place.

    Returns:
        None. Mutates `pyps` in place.

    Status: ACTIVE (C++-wrapper path) -- only called by
    specex.run_specex(). Real production code for the C++ codepath,
    not unused legacy; just a different entry point than the GPU-native
    path.
    """
    import specex._libspecex as spx
    pyps.init_traces(opts)
    fitsfilename = opts.input_psf_filename
    fitsfile     = FITS(fitsfilename,'r')
    xt_hdr = fitsfile['XTRACE'].read_header()
    pyps.trace_ncoeff  = xt_hdr['NAXIS1']
    pyps.nfibers       = xt_hdr['NAXIS2']
    pyps.trace_WAVEMIN = xt_hdr['WAVEMIN']
    pyps.trace_WAVEMAX = xt_hdr['WAVEMAX']
    pyps.TRDEGW        = pyps.trace_ncoeff - 1
    xtrace = fitsfile['XTRACE'].read()
    ytrace = fitsfile['YTRACE'].read()
    pyps.set_trace(xtrace,opts.trace_deg_x,1)
    pyps.set_trace(ytrace,opts.trace_deg_wave,0)
    pyps.synchronize_traces() 
    if 'PSF' in fitsfile:
        psf_hdr = fitsfile['PSF'].read_header()
        pyps.GHDEGX          = psf_hdr['GHDEGX']
        pyps.GHDEGY          = psf_hdr['GHDEGY']
        pyps.mjd             = psf_hdr['MJD']
        pyps.plate_id        = psf_hdr['PLATEID']
        pyps.camera_id       = psf_hdr['CAMERA']
        pyps.arc_exposure_id = psf_hdr['ARCEXP']
        pyps.NPIX_X          = psf_hdr['NPIX_X']
        pyps.NPIX_Y          = psf_hdr['NPIX_Y']
        pyps.hSizeX          = psf_hdr['HSIZEX']
        pyps.hSizeY          = psf_hdr['HSIZEY']
        pyps.FIBERMIN        = psf_hdr['FIBERMIN']
        pyps.FIBERMAX        = psf_hdr['FIBERMAX']
        pyps.table_WAVEMIN   = psf_hdr['WAVEMIN']
        pyps.table_WAVEMAX   = psf_hdr['WAVEMAX']
        pyps.LEGDEG          = psf_hdr['LEGDEG']
        table_col0 = spx.VectorString(); table_col1 = spx.VectorDouble(); table_col2 = spx.VectorInt(); table_col3 = spx.VectorInt()
        col0 = fitsfile['PSF']['PARAM'][:]; col1 = fitsfile['PSF']['COEFF'][:]; col2 = fitsfile['PSF']['LEGDEGX'][:]; col3 = fitsfile['PSF']['LEGDEGW'][:]
        for t in col0: table_col0.append(t)
        for t in col2: table_col2.append(t)
        for t in col3: table_col3.append(t)
        pyps.table_nrows = np.shape(col1)[0]; pyps.nfibers = np.shape(col1)[1]; pyps.ncoeff = np.shape(col1)[2]
        for r in np.arange(pyps.table_nrows):
            for f in np.arange(pyps.nfibers):
                for c in np.arange(pyps.ncoeff): table_col1.append(col1[r,f,c])
        pyps.set_psf(table_col0,table_col1,table_col2,table_col3)

def read_preproc_cpp(opts):
    """Read a preproc arc image (same processing as read_preproc) and wrap it as a C++-side PyImage object for the C++ fitter.

    Args:
        opts: a specex._libspecex.PyOptions object (supplies
            `arc_image_filename`).

    Returns:
        specex._libspecex.PyImage: the image/ivar/mask/rdnoise/header data,
        C++-extension-ready.

    Status: ACTIVE (C++-wrapper path) -- only called by
    specex.run_specex(). Real production code for the C++ codepath,
    not unused legacy; just a different entry point than the GPU-native
    path.
    """
    import specex._libspecex as spx
    ddata = read_image(opts.arc_image_filename)
    ddata['ivar'][ddata['mask'] != 0] = 0.0
    rdnoise_meta = ddata['meta'].get('RDNOISE', 0.0)
    ddata['rdnoise'] = np.full_like(ddata['image'], float(rdnoise_meta))
    hdr = meta2header(ddata['meta'])
    return spx.PyImage(ddata['image'], ddata['ivar'], ddata['mask'], ddata['rdnoise'], hdr)

def get_desi_linelist_file():
    """Resolve the default DESI arc-lamp line list file's path, without needing the caller to know where it's installed.

    Same resolution `desispec.scripts.specex` already uses for the C++
    path: an explicit `$SPECEXDATA` override first, else the file installed
    alongside this package's own `data/` directory (`importlib.resources`,
    works whether this package is a `pip install`, an editable install, or
    a plain `PYTHONPATH`-prepended checkout like this branch's own
    `env_setup.sh` sets up).

    Returns:
        str: absolute path to `specex_linelist_desi.txt`. Not guaranteed to
        exist (e.g. a corrupted install) -- callers that need to fail fast
        on a missing file should check `os.path.exists()` themselves, as
        `py/specex/test/test_specex.py` does.

    Status: ACTIVE -- was previously duplicated ad hoc inline in `main()`
    (`specex.py`) via a hand-rolled `os.path.dirname` climb; added here
    2026-09-13 as a single reusable helper, matching the more robust
    `importlib.resources` approach `desispec.scripts.specex` already used
    independently. `main()` still has its own inline version (not switched
    over in this pass, to avoid touching the CLI's argument-parsing flow
    for an unrelated fix) -- prefer this function for any new code.
    """
    if "SPECEXDATA" in os.environ:
        specexdata = os.environ["SPECEXDATA"]
    else:
        from importlib import resources
        specexdata = resources.files("specex").joinpath("data")
    return os.path.join(str(specexdata), "specex_linelist_desi.txt")

def read_lamp_lines(filename):
    """Parse a lamp line list file (whitespace-separated `name wave score` per line, '#'-prefixed comments and blank/malformed lines skipped) into candidate-generation input.

    Args:
        filename (str): path to the lamp line list file.

    Returns:
        list[dict]: one {'wave': float, 'name': str, 'score': int} per valid
        line, in file order.

    Status: ACTIVE (production default path).
    """
    lines = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'): continue
            parts = line.split()
            if len(parts) < 3: continue
            try:
                wave = float(parts[1]); name = parts[0]; score = int(parts[2])
                lines.append({'wave': wave, 'name': name, 'score': score})
            except: continue
    return lines
