#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <map>

#include <boost/program_options.hpp>

#include <harp.hpp>

#include <specex_message.h>

#include <specex_serialization.h>

#include <specex_trace.h>
#include <specex_spot.h>
#include <specex_spot_array.h>
//#include <specex_spectrograph.h>
#include <specex_lamp_lines_utils.h>
#include <specex_fits.h>
#include <specex_desi_io.h>
#include <specex_psf_fitter.h>

#include <specex_psf_io.h>

using namespace std;
using namespace specex;

//namespace popts = boost::program_options;

#define _GNU_SOURCE 1
#ifndef __USE_GNU
#define __USE_GNU
#endif
#include <portable_fenv.h>

/*
  input format of fits file
  HDU 0 : blank
  HDU 1 ELECTRONS : 2D image of noisy electrons
  - RDNOISE is in header
  HDU 2 IVAR : Inverse variance [1/electrons^2]
  - ivar = 1/(pix.clip(0) + rdnoise**2)
  HDU 3 MASK : 0=good.  Currently all 0.
  HDU 4 XCOEFF : Legendre coefficients for mapping wavelength -> x
  - WAVEMIN, WAVEMAX : domain for mapping to [-1,1] for Legendre polynomials
  - image is coefficients for each fiber
  HDU 5 YCOEFF : Legendre coefficients for mapping wavelength -> y
  HDU 5 TRUE_ELECTRONS : original noiseless image in electrons
  ELECTRONS = poisson(TRUE_ELECTRONS) + gaussian(rdnoise)
 */

int specex_desi_psf_fit_main ( int argc, char *argv[] ) {

  return EXIT_SUCCESS;
  
}

//#include <specex_serialization_implement.h>


