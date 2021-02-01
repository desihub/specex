#ifndef SPECEX_PYPRIOR__H
#define SPECEX_PYPRIOR__H

#include <vector>
#include <string>

#include <harp.hpp>

#include <specex_pyoptions.h>
#include <specex_psf.h>

namespace specex {
  
  class PyPrior : public std::enable_shared_from_this <PyPrior> {

  public :

    typedef std::shared_ptr <PyPrior> pshr;

    map<string,Prior*> priors;
    
    int deal_with_priors( specex::PyOptions );

    PyPrior(){}
    
  };
  
}

#endif
