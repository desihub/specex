#include <iostream>
#include <fstream>

#include <harp.hpp>

#include "specex_spot.h"
#include "specex_trace.h"
#include "specex_spectrograph.h"
#include "specex_lamp_lines_utils.h"
#include "specex_message.h"

using namespace std;


void specex::allocate_spots_of_bundle(vector<specex::Spot_p>& spots, const specex::Spectrograph & spectro, const string& lamp_lines_filename, const specex::TraceSet& traceset, 
				     int fiber_bundle, int fiber_min, int fiber_max, int ymin, int ymax, 
				      const double& min_wavelength, const double& max_wavelength) {
  
  if(fiber_bundle<0 || fiber_bundle >= spectro.number_of_fiber_bundles_per_ccd) {
    SPECEX_ERROR("inproper fiber bundle id");
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
      
      for(int fiber=max(spectro.number_of_fibers_per_bundle*fiber_bundle,fiber_min); fiber<min(spectro.number_of_fibers_per_bundle*(fiber_bundle+1),fiber_max+1); fiber++) {
	
	const specex::Trace& trace = traceset[fiber];
	if(trace.Off()) {
	  //SPECEX_WARNING("Ignore spot in fiber " << fiber << " because mask=" << trace.mask);
	  continue;
	}
	
	if(wave<trace.X_vs_W.xmin || wave>trace.X_vs_W.xmax) 
	  continue;
	if(wave<trace.Y_vs_W.xmin || wave>trace.Y_vs_W.xmax) 
	  continue;
	
	specex::Spot_p spot(new specex::Spot());
	spot->wavelength = wave;
	spot->fiber = fiber;
	spot->fiber_bundle = fiber_bundle;
	spot->xc = trace.X_vs_W.Value(spot->wavelength);
	spot->yc = trace.Y_vs_W.Value(spot->wavelength);
	
	if(spot->yc<ymin) continue;
	if(spot->yc>ymax) continue;
	

	// apply 
	  
	spot->flux = 0;
	spot->eflux = 99;
	
	spot->initial_xc = spot->xc;
	spot->initial_yc = spot->yc;
	spot->initial_flux = spot->flux;
	
	spots.push_back(spot);

      }
    }
  }
  is.close();
}


