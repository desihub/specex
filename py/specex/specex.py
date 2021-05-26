from astropy.io import fits
import fitsio
from fitsio import FITS,FITSHDR
from datetime import datetime

import numpy as np
from specex._libspecex import (PyOptions,PyIO,PyPrior,PyImage,PyPSF,PyFitting,VectorString)
from specex.io import (read_preproc, write_psf, read_psf)

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

    # parse args
    opts.parse(spxargs)

    # set input psf bools
    pyio.set_inputpsf(opts)

     # set Gaussian priors
    pypr.set_priors(opts)    

    # read psf 
    read_psf(opts,pyio,pyps)  

    # read preproc
    pymg = read_preproc(opts) 
    
    # fit psf 
    pyft.fit_psf(opts,pyio,pypr,pymg,pyps) 
    
    # write psf 
    write_psf(pyps,opts,pyio)        

    return 0
