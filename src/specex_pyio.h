#ifndef SPECEX_PYIO__H
#define SPECEX_PYIO__H

#include <vector>
#include <string>

#include <specex_unbls.h>

#include <specex_pyoptions.h>
#include <specex_pyimage.h>
#include <specex_pypsf.h>
#include <specex_psf_fitter.h>

namespace specex {
  
  class PyIO : public std::enable_shared_from_this <PyIO> {

  public :

    typedef std::shared_ptr <PyIO> pshr;

    bool use_input_specex_psf;
    bool psf_change_req;

    int load_psf(      specex::PyOptions, specex::PyPSF&  );
    int write_spots(   specex::PyOptions, specex::PyPSF&  );
    int set_inputpsf(  specex::PyOptions, specex::PyPSF&  );
    
    PyIO()
      : use_input_specex_psf(false)
      , psf_change_req(false)
      {}
    
  };
  
}

#endif
