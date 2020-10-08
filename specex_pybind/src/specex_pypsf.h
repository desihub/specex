#ifndef SPECEX_PYPSF__H
#define SPECEX_PYPSF__H

#include <boost/program_options.hpp>
#include <vector>
#include <string>

#include <harp.hpp>

#include <specex_psf.h>
#include <specex_gauss_hermite_psf.h>

#include <specex_pyoptions.h>

namespace specex {
  
  class PyPSF : public std::enable_shared_from_this <PyPSF> {

  public :

    typedef std::shared_ptr <PyPSF> pshr;

    specex::PSF_p psf;
    
    PyPSF(){}
    
  };
  
}

#endif
