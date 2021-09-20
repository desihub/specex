#include <string>
#include <iostream>
#include <specex_message.h>
#include <specex_spot.h>

#include <specex_pyio.h>
#include <specex_pyimage.h>
#include <specex_psf_proc.h>

using namespace std;

int specex::PyIO::load_psf(specex::PyOptions opts, specex::PyPSF& pypsf){

  load_psf_work(pypsf.psf);  
  
  return EXIT_SUCCESS;

}

int specex::PyIO::write_spots(specex::PyOptions opts, specex::PyPSF& pypsf){

  vector <Spot_p> fitted_spots = pypsf.fitted_spots;

  // future location of spot writing
  
  return EXIT_SUCCESS;

}

int specex::PyIO::set_inputpsf(specex::PyOptions opts){
  
  use_input_specex_psf = true;

  psf_change_req |= (! opts.half_size_x_def);
  psf_change_req |= (! opts.half_size_y_def);
  psf_change_req |= (! opts.gauss_hermite_deg_def);
  psf_change_req |= (! opts.gauss_hermite_deg2_def);
  //psf_change_req |= (! opts.legendre_deg_wave_def);
  psf_change_req |= (! opts.legendre_deg_x_def);
  psf_change_req |= (! opts.trace_deg_wave_def);
  psf_change_req |= (! opts.trace_deg_x_def);
  
  if(psf_change_req && use_input_specex_psf) {
    SPECEX_WARNING("option(s) were given to specify the psf properties, so we cannot use the input PSF parameters as a starting point (except for the trace coordinates)");
    use_input_specex_psf = false;
  }
  if(use_input_specex_psf) {
    SPECEX_INFO("will start from input psf parameters")
  }

  return EXIT_SUCCESS;
  
}

