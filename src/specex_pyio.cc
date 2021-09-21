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

int specex::PyIO::set_inputpsf(specex::PyOptions opts, specex::PyPSF& pypsf){

  // use input PSF by default
  use_input_specex_psf = true;

  // check if PSF parameters (whether specified on command line or default value)
  // are different from those in the input PSF; if so, don'tuse PSF as starting point
  psf_change_req |= (pypsf.hSizeX != opts.half_size_x);
  psf_change_req |= (pypsf.hSizeY != opts.half_size_y);
  psf_change_req |= (pypsf.GHDEGX != opts.gauss_hermite_deg);
  psf_change_req |= (pypsf.GHDEGY != opts.gauss_hermite_deg);
  psf_change_req |= (pypsf.TRDEGW != opts.trace_deg_wave);
  psf_change_req |= (pypsf.LEGDEG != opts.legendre_deg_wave);

  printf("IODEBUG: legdeg = %d and legendre_deg_wave = %d for %s\n",pypsf.LEGDEG,opts.legendre_deg_wave,opts.input_psf_filename);
  fflush(stdout);
  if(psf_change_req) {
    SPECEX_WARNING("specified and/or default psf properties differ from those of the input PSF, so we cannot use the input PSF parameters as a starting point (except for the trace coordinates)");
    use_input_specex_psf = false;
  }

  // these parameters are not typically specified and/or deprecated
  // if they are specified (param_def = False) don't use PSF as starting point
  psf_change_req |= (! opts.gauss_hermite_deg2_def);
  psf_change_req |= (! opts.legendre_deg_x_def);
  psf_change_req |= (! opts.trace_deg_x_def);

  if(psf_change_req) {
    SPECEX_WARNING("parameters not typically changed and/or deprecated have been specified, so we cannot use the input PSF parameters as a starting point (except for the trace coordinates)");
    use_input_specex_psf = false;
  }

  if(use_input_specex_psf) {
    SPECEX_INFO("will start from input psf parameters")
  }

  return EXIT_SUCCESS;

}
