#ifndef SPECEX_PSFPY__H
#define SPECEX_PSFPY__H

#include <boost/program_options.hpp>
#include <vector>
#include <string>

#include <harp.hpp>
#include "specex_fits.h"
#include "specex_image_data.h"

namespace specex {
  
  class PSFPy : public std::enable_shared_from_this <PSFPy> {

  public :

    typedef std::shared_ptr <PSFPy> pshr;
    
    specex::image_data coeff2d_x, coeff2d_y;
    int    FIBERMIN, FIBERMAX;
    double  WAVEMIN,  WAVEMAX;

    specex::FitsTable table;

    std::map<std::string, std::string> tablekeys_comment;
    std::map<std::string, std::string> tablekeys_string;
    std::map<std::string, int>         tablekeys_int;
    std::map<std::string, double>      tablekeys_double;
    
    specex::image_data GetCoeff2d(bool);  
    
    void SetCoeff2d(specex::image_data, bool);    
    
    PSFPy(){}

  };
  
}

#endif
