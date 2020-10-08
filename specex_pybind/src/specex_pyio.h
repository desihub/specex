#ifndef SPECEX_PYIO__H
#define SPECEX_PYIO__H

#include <boost/program_options.hpp>
#include <vector>
#include <string>

#include <harp.hpp>

#include <specex_options.h>

namespace specex {
  
  class PyIO : public std::enable_shared_from_this <PyIO> {

  public :

    typedef std::shared_ptr <PyIO> pshr;

    bool use_input_specex_psf;
    bool psf_change_req;

    int check_input_psf(specex::Options opts);

    PyIO()
      : use_input_specex_psf(false)
      , psf_change_req(false)
      {}
    
  };
  
}

#endif
