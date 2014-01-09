#ifndef SPECEX_PSF_IO__H
#define SPECEX_PSF_IO__H

#include <string>

//
// need to be included before main() :
//   #included <specex_psf_io.h> 
//   #included <specex_serialisation.h> 
//

namespace specex {

  void read_psf_xml(specex::PSF_p& psf, const std::string& filename);

  void write_psf_xml(const specex::PSF_p psf, const std::string& filename);
  void write_psf_fits_image(const specex::PSF_p psf, const string& filename, const int fiber, const double& wavelength, int bundle, int oversampling=1);
  
  
  

}

#endif
