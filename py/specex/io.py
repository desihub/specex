import os
import numpy as np
import fitsio
from fitsio import FITS, FITSHDR
from datetime import datetime
from .psf import PSF, PSF_Params
from .math import SparseLegendre2DPol, Legendre1DPol, legendre_pol_jnp

def meta2header(meta):
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

def get_sparse_nz(xdeg, ydeg):
    nz = []
    for j in range(ydeg + 1):
        for i in range(xdeg + 1):
            if i == 0: nz.append(i + j*(xdeg + 1))
            elif i == 1 and j < 2: nz.append(i + j*(xdeg + 1))
            elif i > 1 and j == 0: nz.append(i + j*(xdeg + 1))
    return nz

def load_python_psf(filename, opts):
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
    import fitsio
    fin = fitsio.FITS(input_template)
    xtrace_out = fin['XTRACE'].read().astype(np.float64)
    ytrace_out = fin['YTRACE'].read().astype(np.float64)
    psf_table = fin['PSF'].read()
    psf_hdr = fin['PSF'].read_header()
    param_names = [p.strip() for p in psf_table['PARAM']]
    name_to_idx = {name: i for i, name in enumerate(param_names)}
    xdeg_b = 1; wdeg_b = 3; nz_b = get_sparse_nz(xdeg_b, wdeg_b)
    for bid, res in bundle_results.items():
        fmin, fmax = bid * 25, (bid + 1) * 25 - 1
        pc = res['psf_coeffs']; tc = res['trace_coeffs']
        rf = 2 * (np.arange(fmin, fmax + 1) - fmin) / (fmax - fmin) - 1
        poly_f = np.stack([legendre_pol_jnp(i, rf) for i in range(xdeg_b + 1)], axis=0)
        for k_nz, k_lin in enumerate(nz_b):
            i_p, j_p = k_lin % 2, k_lin // 2
            if j_p < xtrace_out.shape[1]:
                xtrace_out[fmin:fmax+1, j_p] += tc[0, k_nz] * poly_f[i_p]
                ytrace_out[fmin:fmax+1, j_p] += tc[1, k_nz] * poly_f[i_p]
        for row in range(len(param_names)):
            psf_table['COEFF'][row, fmin:fmax+1, :] = 0.0
            if param_names[row] == 'GH-0-0': psf_table['COEFF'][row, fmin:fmax+1, 0] = 1.0
        for i_par in range(55):
            if i_par == 0: pname = 'GHSIGX'
            elif i_par == 1: pname = 'GHSIGY'
            elif 2 <= i_par <= 49:
                idx_gh = i_par - 2; gh_i = (idx_gh + 1) % 7; gh_j = (idx_gh + 1) // 7; pname = f'GH-{gh_i}-{gh_j}'
            else: pname = ['TAILAMP', 'TAILCORE', 'TAILXSCA', 'TAILYSCA', 'TAILINDE'][i_par - 50]
            idx = name_to_idx.get(pname)
            if idx is not None:
                for k_nz, k_lin in enumerate(nz_b):
                    i_p, j_p = k_lin % 2, k_lin // 2
                    if j_p < psf_table['COEFF'].shape[2]:
                        psf_table['COEFF'][idx, fmin:fmax+1, j_p] += pc[i_par, k_nz] * poly_f[i_p]
        psf_hdr[f'B{bid:02d}RCHI2'] = res['chi2'] / (120000.0)
    if os.path.exists(filename): os.remove(filename)
    fout = fitsio.FITS(filename, 'rw')
    fout.write(xtrace_out, header=fin['XTRACE'].read_header(), extname='XTRACE')
    fout.write(ytrace_out, header=fin['YTRACE'].read_header(), extname='YTRACE')
    fout.write(psf_table, header=psf_hdr, extname='PSF')
    fout.close()

def write_psf(pyps, opts, pyio):
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
    f = fitsio.FITS(filename)
    image = f['IMAGE'].read().astype(np.float64); ivar = f['IVAR'].read().astype(np.float64); mask = f['MASK'].read().astype(np.int32); meta = f['IMAGE'].read_header()
    return {'image': image, 'ivar': ivar, 'mask': mask, 'meta': meta}

def read_preproc(opts):
    ddata = read_image(opts.arc_image_filename)
    ddata['ivar'][ddata['mask'] != 0] = 0.0
    rdnoise_meta = ddata['meta'].get('RDNOISE', 0.0)
    ddata['rdnoise'] = np.full_like(ddata['image'], float(rdnoise_meta))
    return ddata

def read_psf(opts, pyps):
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
    import specex._libspecex as spx
    ddata = read_image(opts.arc_image_filename)
    ddata['ivar'][ddata['mask'] != 0] = 0.0
    rdnoise_meta = ddata['meta'].get('RDNOISE', 0.0)
    ddata['rdnoise'] = np.full_like(ddata['image'], float(rdnoise_meta))
    hdr = meta2header(ddata['meta'])
    return spx.PyImage(ddata['image'], ddata['ivar'], ddata['mask'], ddata['rdnoise'], hdr)

def read_lamp_lines(filename):
    lines = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            parts = line.split()
            if len(parts) < 3: continue
            try:
                # C++ istringstream >> ion >> wave >> score >> intensity
                # If line is "#ArI 7725.887 1 15000", ion is "#ArI", wave is 7725.887
                wave = float(parts[1]); name = parts[0]; score = int(parts[2])
                lines.append({'wave': wave, 'name': name, 'score': score})
            except: continue
    return lines
