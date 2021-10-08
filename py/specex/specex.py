import numpy as np
from specex._libspecex import (PyOptions,PyIO,PyPrior,PyPSF,PyFitting,VectorString)
from specex.io import (read_preproc, write_psf, read_psf)

class timer:
    def __init__(self, **kwargs):
        self.tinit = time.time()
    def report(self,message):
        dt = time.time() - self.tinit
        print('specex-timer: ',message.ljust(30), dt,' seconds',flush=True);

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
    opts.parse(spxargs)

    # read psf
    read_psf(opts,pyps)

    # set input psf bools
    pyio.set_inputpsf(opts,pyps)

     # set Gaussian priors
    pypr.set_priors(opts)

    # read preproc
    pymg = read_preproc(opts) 

    # fit psf 
    fit_timer=timer()
    pyft.fit_psf(opts,pyio,pypr,pymg,pyps) 
    fit_timer.report('PSF fit for '+opts.arc_image_filename[-16:-14]+' bundles '
                     +str(opts.first_fiber_bundle)+' to '+str(opts.last_fiber_bundle))
    
    # write psf 
    write_psf(pyps,opts,pyio)        

    return 0
