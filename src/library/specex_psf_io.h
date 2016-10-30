#ifndef SPECEX_PSF_IO__H
#define SPECEX_PSF_IO__H

#include <string>



namespace specex {

  void read_psf_xml(specex::PSF_p& psf, const std::string& filename);
  void write_psf_xml(const specex::PSF_p psf, const std::string& filename);
  
  void write_psf_fits_image(const specex::PSF_p psf, const string& filename, 
			    const int fiber, const double& wavelength, int oversampling=1);

  void write_xtrace_fits_hdu(const specex::PSF& psf, fitsfile *fp, int hdu);
  void write_ytrace_fits_hdu(const specex::PSF& psf, fitsfile *fp, int hdu);
  
  void write_psf_fits(const specex::PSF_p psf, fitsfile* fp, int first_hdu=1);
  void write_psf_fits(const specex::PSF_p psf, const string& filename);
  
  void read_psf_fits(specex::PSF_p& psf,  fitsfile* fp, int first_hdu=1);
  void read_psf_fits(specex::PSF_p& psf, const string& filename);


  

  void write_spots_xml(const std::vector<specex::Spot_p>& spots, const std::string& filename);
  
}

#endif
