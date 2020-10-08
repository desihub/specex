#include <string>
#include <iostream>
#include <specex_message.h>
#include <specex_pyio.h>

using namespace std;

int specex::PyIO::check_input_psf(specex::Options opts){
  
  // check input PSF type, can be from boot calib with only the traces
  if (opts.input_psf_filename.find(".xml") != std::string::npos) { // xml file
    SPECEX_INFO("Input PSF file is xml");
    use_input_specex_psf = true;
  }else{ // fits file, look at header  
    try {
      fitsfile * fp;
      harp::fits::open_read(fp,opts.input_psf_filename);
      int status = 0;
      int first_hdu = 1;
      fits_movabs_hdu ( fp, first_hdu, NULL, &status ); harp::fits::check ( status );
      string psftype; harp::fits::key_read(fp,"PSFTYPE",psftype);
      SPECEX_INFO("Input PSF type = " << psftype);
      use_input_specex_psf = (psftype=="GAUSS-HERMITE"); 
    } catch (harp::exception) {
      SPECEX_WARNING("Could not read PSF type in " << opts.input_psf_filename);
    }
  }

  psf_change_req |= (! opts.vm["half-size-x"].defaulted());
  psf_change_req |= (! opts.vm["half-size-y"].defaulted());
  psf_change_req |= (! opts.vm["gauss-hermite-deg"].defaulted());
  psf_change_req |= (! opts.vm["gauss-hermite-deg2"].defaulted());
  psf_change_req |= (! opts.vm["legendre-deg-wave"].defaulted());
  psf_change_req |= (! opts.vm["legendre-deg-x"].defaulted());
  psf_change_req |= (! opts.vm["trace-deg-wave"].defaulted());
  psf_change_req |= (! opts.vm["trace-deg-x"].defaulted());
  
  if(psf_change_req && use_input_specex_psf) {
    SPECEX_WARNING("option(s) were given to specify the psf properties, so we cannot use the input PSF parameters as a starting point (except for the trace coordinates)");
    use_input_specex_psf = false;
  }
  if(use_input_specex_psf) {
    SPECEX_INFO("will start from input psf parameters")
  }

  return EXIT_SUCCESS;
  
}

