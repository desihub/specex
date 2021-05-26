#ifndef SPECEX_PYPRIOR__H
#define SPECEX_PYPRIOR__H

#include <vector>
#include <string>

#include <unbls.h>

#include <specex_pyoptions.h>
#include <specex_psf.h>

namespace specex {
  
  class PyPrior : public std::enable_shared_from_this <PyPrior> {

  public :

    typedef std::shared_ptr <PyPrior> pshr;

    map<string,Prior*> priors;
    
    int set_priors( specex::PyOptions );

    PyPrior(){}
    
  };
  
}

#endif
