import os
import numpy as np
from .psf import PSF, PSF_Params
from .math import SparseLegendre2DPol, Legendre1DPol
import fitsio

def load_python_psf(filename):
    """
    Pure Python loader for specex PSF files.
    """
    f = fitsio.FITS(filename)
    
    # 1. Traces
    xtrace = f['XTRACE'].read()
    ytrace = f['YTRACE'].read()
    xt_hdr = f['XTRACE'].read_header()
    
    psf = PSF()
    psf.fiber_min = xt_hdr['FIBERMIN']
    psf.fiber_max = xt_hdr['FIBERMAX']
    
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
        
        # Initializing Bundle 5 (example)
        bundle_id = 5
        b_fmin, b_fmax = bundle_id * 25, (bundle_id + 1) * 25 - 1
        bundle = PSF_Params(bundle_id, b_fmin, b_fmax)
        bundle.param_names = [p.strip() for p in table['PARAM']]
        bundle.param_models = {}
        
        # Each row in table['COEFF'] is (Nfibers, Ncoeff_wave)
        for i, name in enumerate(bundle.param_names):
            bundle.param_models[name] = [
                Legendre1DPol(deg=hdr['LEGDEG'], xmin=hdr['WAVEMIN'], xmax=hdr['WAVEMAX'], coeff=table['COEFF'][i][fib])
                for fib in range(500)
            ]
        psf.params_of_bundles[bundle_id] = bundle
    
    return psf

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

def read_image(filename):
    f = fitsio.FITS(filename)
    image = f['IMAGE'].read().astype(np.float64)
    ivar = f['IVAR'].read().astype(np.float64)
    mask = f['MASK'].read().astype(np.int32)
    meta = f['IMAGE'].read_header()
    return {'image': image, 'ivar': ivar, 'mask': mask, 'meta': meta}

def read_preproc(opts):
    dsmg = read_image(opts.arc_image_filename)
    dsmg['ivar'][dsmg['mask']!=0] = 0.0
    import specex._libspecex as spx
    hdr = meta2header(dsmg['meta'])
    pymg = spx.PyImage(dsmg['image'], dsmg['ivar'], dsmg['mask'], 0.0, hdr) # dummy rdnoise
    return pymg

def read_psf(opts, pyps):
    # Header only loading into pyps for dimensions
    f = fitsio.FITS(opts.input_psf_filename)
    xt_hdr = f['XTRACE'].read_header()
    pyps.trace_ncoeff  = xt_hdr['NAXIS1']
    pyps.nfibers       = xt_hdr['NAXIS2']
    pyps.trace_WAVEMIN = xt_hdr['WAVEMIN']
    pyps.trace_WAVEMAX = xt_hdr['WAVEMAX']
    pyps.FIBERMIN      = xt_hdr['FIBERMIN']
    pyps.FIBERMAX      = xt_hdr['FIBERMAX']
    if 'PSF' in f:
        p_hdr = f['PSF'].read_header()
        pyps.GHDEGX = p_hdr['GHDEGX']
        pyps.GHDEGY = p_hdr['GHDEGY']
        pyps.hSizeX = p_hdr['HSIZEX']
        pyps.hSizeY = p_hdr['HSIZEY']
        pyps.table_nrows = p_hdr['NPARAMS']
        pyps.ncoeff = p_hdr['LEGDEG'] + 1
        pyps.table_WAVEMIN = p_hdr['WAVEMIN']
        pyps.table_WAVEMAX = p_hdr['WAVEMAX']
    return
