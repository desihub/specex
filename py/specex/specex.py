import numpy as np
from specex._libspecex import (PyOptions,PyIO,PyPrior,PyPSF,PyFitting,VectorString)
from specex.io import (read_preproc, write_psf, read_psf)

def run_specex(com):

    # instantiate specex C++ objects exposed to python        
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
    retval = opts.parse(spxargs)
    if retval != 0: return retval

    # read psf
    read_psf(opts,pyps)

    # set input psf bools
    pyio.set_inputpsf(opts,pyps)

     # set Gaussian priors
    pypr.set_priors(opts)

    # read preproc
    pymg = read_preproc(opts) 
    
    # fit psf 
    retval = pyft.fit_psf(opts,pyio,pypr,pymg,pyps) 
    
    # write psf 
    write_psf(pyps,opts,pyio)        

    return retval
