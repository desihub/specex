#include <string>
#include <iostream>
#include <specex_pypsf.h>
using namespace std;

specex::image_data specex::PyPSF::get_trace(std::string axis){

  ///specex::PSF_p psf = this->psf;
  
  if(axis == "x"){
    return this->psf->pydata.coeff2d_x;
  }
  else{
    return this->psf->pydata.coeff2d_y;
  }
  
}

