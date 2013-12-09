#include <iostream>
#include <fstream>

#include "harp.hpp"

#include "specex_spot.h"
#include "specex_trace.h"
#include "specex_spectrograph.h"
#include "specex_lamp_lines_utils.h"

using namespace std;


void specex::allocate_spots_of_bundle(vector<specex::Spot>& spots, const specex::Spectrograph & spectro, const string& lamp_lines_filename, const specex::TraceSet& traceset, 
				     int fiber_bundle, int ymin, int ymax, 
				      const double& min_wavelength, const double& max_wavelength) {
  
  if(fiber_bundle<0 || fiber_bundle >= spectro.number_of_fiber_bundles_per_ccd) {
    HARP_THROW("inproper fiber bundle id");
  }

  spots.clear();
  
  string word;
  ifstream is(lamp_lines_filename.c_str());
  while(is >> word) {  
    if(word=="arclineid") {
      double wave;
      if(!(is >> wave)) continue;
      
      // for test to go faster
      if(wave<min_wavelength || wave>max_wavelength) continue;
      
      for(int fiber=spectro.number_of_fibers_per_bundle*fiber_bundle; fiber<spectro.number_of_fibers_per_bundle*(fiber_bundle+1); fiber++) {
	
	specex::Spot spot;
	spot.log10_wavelength = log10(wave);
	spot.fiber = fiber;
	spot.xc = traceset[fiber].X_vs_lW.Value(spot.log10_wavelength);
	spot.yc = traceset[fiber].Y_vs_lW.Value(spot.log10_wavelength);
	
	if(spot.yc<ymin) continue;
	if(spot.yc>ymax) continue;
	

	// apply 
	  
	spot.flux = 0;
	spot.eflux = 99;
	
	spot.initial_xc = spot.xc;
	spot.initial_yc = spot.yc;
	spot.initial_flux = spot.flux;
	
	//spot.PSFname
	//spot.PSFParams

	spots.push_back(spot);

      }
    }
  }
  is.close();
}


