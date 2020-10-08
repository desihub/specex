#ifndef SPECEX_PYIMAGE__H
#define SPECEX_PYIMAGE__H

#include <boost/program_options.hpp>
#include <vector>
#include <string>

#include <harp.hpp>

#include <specex_pyoptions.h>
#include <specex_image_data.h>

namespace specex {
  
  class PyImage : public std::enable_shared_from_this <PyImage> {

  public :

    typedef std::shared_ptr <PyImage> pshr;

    std::map<std::string,std::string> header;
    image_data image,weight,mask,rdnoise;

    PyImage(){}
    
  };
  
}

#endif
