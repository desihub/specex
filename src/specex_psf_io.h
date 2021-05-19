#ifndef SPECEX_PSF_IO__H
#define SPECEX_PSF_IO__H

#include <string>



namespace specex {

  // to write an image of the psf
  void write_psf_fits_image(const specex::PSF_p psf, const string& filename, 
			    const int fiber, const double& wavelength, int oversampling=1);

  // to write the spots coordinates, fluxes ... etc, this is for debugging
  void write_spots_xml(const std::vector<specex::Spot_p>& spots, const std::string& filename);
  
  void write_psf_xml(const specex::PSF_p psf, const std::string& filename);
  void write_psf_fits_dummy(const string& filename);

  // read routines
  void read_xtrace_fits_hdu(specex::PSF_p psf, fitsfile *fp, int hdu, int requested_deg=0);
  void read_ytrace_fits_hdu(specex::PSF_p psf, fitsfile *fp, int hdu, int requested_deg=0);
  void synchronize_traces(specex::PSF_p psf);
  void read_traceset_fits(specex::PSF_p psf, fitsfile * fp, int degx=0, int degy=0);
  void read_traceset_fits(specex::PSF_p psf, const string& filename, int degx=0, int degy=0);
  
  void read_psf_xml(specex::PSF_p& psf, const std::string& filename);
  void read_psf_fits(specex::PSF_p& psf, const string& filename);
  void read_psf_gen(specex::PSF_p& psf, const string& filename);
  
  
}

#endif
