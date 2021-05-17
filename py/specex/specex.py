from astropy.io import fits
import fitsio
from fitsio import FITS,FITSHDR
from datetime import datetime

import numpy as np
from specex._libspecex import (PyOptions,PyIO,PyPrior,PyImage,PyPSF,PyFitting,VectorString)
from specex.io import (read_preproc, write_psf)

def run_specex(com):

    # instantiate specex c++ objects exposed to python        
    opts = PyOptions() # input options
    pyio = PyIO()      # IO options and methods
    pypr = PyPrior()   # Gaussian priors
    pyps = PyPSF()     # psf data
    pyft = PyFitting() # psf fitting
    
    # copy com to opaque pybind VectorString object spxargs
    spxargs = VectorString()
    for strs in com:
        spxargs.append(strs)
        
    opts.parse(spxargs)     # parse args
    pyio.set_inputpsf(opts) # set input psf bools
    pypr.set_priors(opts)   # set Gaussian priors

    pymg = read_preproc(opts) # read preproc image 
    pyio.read_psf(opts,pyps)  # read psf 
    
    pyft.fit_psf(opts,pyio,pypr,pymg,pyps) # fit psf 
    
    pyio.load_psf(opts,pyps)    # load psf
    write_psf(pyps,opts)        # write psf 

    return 0
