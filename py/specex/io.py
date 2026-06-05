import os
import numpy as np
import fitsio
from fitsio import FITS, FITSHDR
from datetime import datetime
from .psf import PSF, PSF_Params
from .math import SparseLegendre2DPol, Legendre1DPol

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
    """
    Pure Python loader for specex PSF files.
    Ensures native endianness for JAX compatibility.
    """
    f = fitsio.FITS(filename)
    
    psf = PSF()
    xt_hdr = f['XTRACE'].read_header()
    psf.fiber_min = xt_hdr['FIBERMIN']
    psf.fiber_max = xt_hdr['FIBERMAX']
    
    # 1. Traces (1D Legendre per fiber)
    xtrace = f['XTRACE'].read().astype(np.float64)
    ytrace = f['YTRACE'].read().astype(np.float64)
    for fib in range(psf.fiber_min, psf.fiber_max + 1):
        idx = fib - psf.fiber_min
        psf.fiber_traces[fib] = {
            'X_vs_W': Legendre1DPol(deg=xt_hdr['NAXIS1']-1, xmin=xt_hdr['WAVEMIN'], xmax=xt_hdr['WAVEMAX'], coeff=xtrace[idx]),
            'Y_vs_W': Legendre1DPol(deg=xt_hdr['NAXIS1']-1, xmin=xt_hdr['WAVEMIN'], xmax=xt_hdr['WAVEMAX'], coeff=ytrace[idx])
        }

    # 2. PSF Parameters
    if 'PSF' in f:
        table = f['PSF'].read()
        hdr = f['PSF'].read_header()
        
        psf.gh_psf.degree = hdr['GHDEGX']
        psf.h_size_x = hdr['HSIZEX']
        psf.h_size_y = hdr['HSIZEY']
        psf.gain = hdr['GAIN']
        psf.readout_noise = hdr['READNOIS']
        
        bundle_ids = []
        for key in hdr.keys():
            if key.startswith('B') and key.endswith('NDATA'):
                bundle_ids.append(int(key[1:3]))
        
        coeffs_all = table['COEFF'].astype(np.float64)
        param_names = [p.strip() for p in table['PARAM']]

        for bid in bundle_ids:
            b_fmin, b_fmax = bid * 25, (bid + 1) * 25 - 1
            bundle = PSF_Params(bid, b_fmin, b_fmax)
            bundle.param_names = param_names
            bundle.param_models = {}
            for i, name in enumerate(bundle.param_names):
                bundle.param_models[name] = [
                    Legendre1DPol(deg=hdr['LEGDEG'], xmin=hdr['WAVEMIN'], xmax=hdr['WAVEMAX'], coeff=coeffs_all[i][fib])
                    for fib in range(500)
                ]
            psf.params_of_bundles[bid] = bundle
    
    return psf

def read_image(filename):
    f = fitsio.FITS(filename)
    image = f['IMAGE'].read().astype(np.float64)
    ivar = f['IVAR'].read().astype(np.float64)
    mask = f['MASK'].read().astype(np.int32)
    meta = f['IMAGE'].read_header()
    return {'image': image, 'ivar': ivar, 'mask': mask, 'meta': meta}

def read_preproc(opts):
    """
    Reads and pre-processes image data, returning a clean dictionary of NumPy arrays.
    """
    ddata = read_image(opts.arc_image_filename)
    # Masking logic
    ddata['ivar'][ddata['mask'] != 0] = 0.0
    
    # Metadata for readout noise
    rdnoise_meta = ddata['meta'].get('RDNOISE', 0.0)
    ddata['rdnoise'] = np.full_like(ddata['image'], float(rdnoise_meta))
    
    return ddata

def create_cpp_image(ddata):
    """
    Bridge helper to create a legacy C++ PyImage if needed for baseline comparison.
    """
    import specex._libspecex as spx
    hdr = meta2header(ddata['meta'])
    return spx.PyImage(ddata['image'], ddata['ivar'], ddata['mask'], ddata['rdnoise'], hdr)

def read_psf(opts, pyps):
    """
    Original behavior for C++ baseline fits.
    """
    pyps.init_traces(opts)
    fitsfilename = opts.input_psf_filename
    fitsfile     = FITS(fitsfilename,'r')
    xtrace_header = fitsfile['XTRACE'].read_header()
    ytrace_header = fitsfile['YTRACE'].read_header()
    pyps.trace_ncoeff  = xtrace_header['NAXIS1']
    pyps.nfibers       = xtrace_header['NAXIS2']
    pyps.trace_WAVEMIN = xtrace_header['WAVEMIN']
    pyps.trace_WAVEMAX = xtrace_header['WAVEMAX']
    pyps.TRDEGW        = pyps.trace_ncoeff - 1
    xtrace = fitsfile['XTRACE'].read()
    ytrace = fitsfile['YTRACE'].read()
    pyps.set_trace(xtrace,opts.trace_deg_x   ,1)
    pyps.set_trace(ytrace,opts.trace_deg_wave,0)
    pyps.synchronize_traces() 
    if 'PSF' in fitsfile:
        psf_header = fitsfile['PSF'].read_header()
        pyps.GHDEGX          = psf_header['GHDEGX']
        pyps.GHDEGY          = psf_header['GHDEGY']
        pyps.mjd             = psf_header['MJD']
        pyps.plate_id        = psf_header['PLATEID']
        pyps.camera_id       = psf_header['CAMERA']
        pyps.arc_exposure_id = psf_header['ARCEXP']
        pyps.NPIX_X          = psf_header['NPIX_X']
        pyps.NPIX_Y          = psf_header['NPIX_Y']
        pyps.hSizeX          = psf_header['HSIZEX']
        pyps.hSizeY          = psf_header['HSIZEY']
        pyps.FIBERMIN        = psf_header['FIBERMIN']
        pyps.FIBERMAX        = psf_header['FIBERMAX']
        pyps.table_WAVEMIN   = psf_header['WAVEMIN']
        pyps.table_WAVEMAX   = psf_header['WAVEMAX']
        pyps.LEGDEG          = psf_header['LEGDEG']
        
        # Table data handling (requires spx imports internally)
        import specex._libspecex as spx
        table_col0 = spx.VectorString()
        table_col1 = spx.VectorDouble()
        table_col2 = spx.VectorInt()
        table_col3 = spx.VectorInt()
        col0 = fitsfile['PSF']['PARAM'][:]
        col1 = fitsfile['PSF']['COEFF'][:]
        col2 = fitsfile['PSF']['LEGDEGX'][:]
        col3 = fitsfile['PSF']['LEGDEGW'][:]
        for t in col0: table_col0.append(t)
        for t in col2: table_col2.append(t)
        for t in col3: table_col3.append(t)
        pyps.table_nrows = np.shape(col1)[0]
        pyps.nfibers     = np.shape(col1)[1]
        pyps.ncoeff      = np.shape(col1)[2]
        for r in np.arange(pyps.table_nrows):
            for f in np.arange(pyps.nfibers):
                for c in np.arange(pyps.ncoeff):
                    table_col1.append(col1[r,f,c])
        pyps.set_psf(table_col0,table_col1,table_col2,table_col3)
    return

def get_desi_linelist_file():
    specexdata = os.environ.get('SPECEXDATA', '')
    if not specexdata:
        from importlib import resources
        specexdata = resources.files('specex').joinpath('data')
    return os.path.join(specexdata,'specex_linelist_desi.txt')

def read_lamp_lines(filename):
    """
    Reads lamp line ASCII files.
    """
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
