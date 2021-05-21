#ifndef SPECEX_PYIMAGE__H
#define SPECEX_PYIMAGE__H

#include <boost/program_options.hpp>
#include <vector>
#include <string>

#include <specex_unbls.h>

#include <specex_pyoptions.h>
#include <specex_image_data.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

namespace specex {
  
  class PyImage : public std::enable_shared_from_this <PyImage> {

  public :

    typedef std::shared_ptr <PyImage> pshr;

    std::map<std::string,std::string> header;
    image_data image,weight,mask,rdnoise;

    PyImage(){}
    PyImage(py::array,py::array,py::array,py::array,
	    std::map<std::string,std::string>);
    
    std::vector<double> get_data(std::string);
    std::map<std::string,std::string> get_header();    
    
  };
  
}

#endif
