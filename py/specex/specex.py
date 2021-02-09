from astropy.io import fits
import fitsio
from fitsio import FITS,FITSHDR
from datetime import datetime

import numpy as np
from specex._libspecex import (PyOptions,PyIO,PyPrior,PyImage,PyPSF,PyFitting,VectorString)
from specex.io import (read_desi_ppimage, write_psf)

def run_specex(com):
    
    # instantiate specex c++ objects exposed to python        
    opts = PyOptions() # input options
    pyio = PyIO()      # IO options and methods
    pypr = PyPrior()   # Gaussian priors
    pyps = PyPSF()     # psf data
    pyft = PyFitting() # psf fitting
    
    # copy com to opaque pybind VectorString object args
    spxargs = VectorString()
    for strs in com:
        spxargs.append(strs)
        
    opts.parse(spxargs)         # parse args
    pyio.check_input_psf(opts)  # set input psf bools
    pypr.deal_with_priors(opts) # set Gaussian priors
        
    pymg = read_desi_ppimage(opts) # read preproc images (desispec)        
    pyio.read_psf_data(opts,pyps)      # read psf (specex)        
    
    pyft.fit_psf(opts,pyio,pypr,pymg,pyps) # fit psf (specex)
    
    pyio.prepare_psf(opts,pyps) # prepare psf (specex)
    write_psf(pyps,opts)        # write psf (fitsio)
    pyio.write_spots(opts,pyps) # write spots

    return 0
