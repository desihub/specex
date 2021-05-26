#include <iostream>
#include <fstream>

#include <specex_unbls.h>

#include "specex_spot.h"
#include "specex_trace.h"
//#include "specex_spectrograph.h"
#include "specex_lamp_lines_utils.h"
#include "specex_message.h"

using namespace std;


void specex::allocate_spots_of_bundle(vector<specex::Spot_p>& spots, const string& lamp_lines_filename, const specex::TraceSet& traceset, 
				      int fiber_bundle, int fiber_min, int fiber_max, int ymin, int ymax, 
				      const double& min_wavelength, const double& max_wavelength) {
  
  /*
    if(fiber_bundle<0 || fiber_bundle >= spectro.number_of_fiber_bundles_per_ccd) {
    SPECEX_ERROR("incorrect fiber bundle id : " << fiber_bundle);
    }
  */

  spots.clear();
  
  string ion;
  double wave;
  int score;
  double intensity;
  ifstream is(lamp_lines_filename.c_str());
  if ( ! is.good()) {
    SPECEX_ERROR("cannot open file " << lamp_lines_filename);
  }
  SPECEX_INFO("reading " << lamp_lines_filename);
  SPECEX_DEBUG("allocating spots in wavelength range " << min_wavelength << " " << max_wavelength);
  
  int nlines=0;
  string line;
  while (std::getline(is, line)) {
    std::istringstream iss(line);
    if( !( iss >> ion >> wave >> score >> intensity) ) continue;
    
    SPECEX_DEBUG(ion << " " << wave << " " << score << " " << intensity);
    
    if(wave<min_wavelength || wave>max_wavelength) continue;
    if(score<1 or score>4) {
      SPECEX_WARNING("ignore emission line " << wave << " with score = " << score);
      continue;
    }
    nlines ++;
    //for(int fiber=max(spectro.number_of_fibers_per_bundle*fiber_bundle,fiber_min); fiber<min(spectro.number_of_fibers_per_bundle*(fiber_bundle+1),fiber_max+1); fiber++) {
    for(int fiber=fiber_min; fiber<=fiber_max; fiber++) {
      
      const specex::Trace& trace = traceset.find(fiber)->second; //traceset[fiber];
      if(trace.Off()) {
	//SPECEX_WARNING("Ignore spot in fiber " << fiber << " because mask=" << trace.mask);
	continue;
      }
      
      if(wave<trace.X_vs_W.xmin || wave>trace.X_vs_W.xmax) {
	//SPECEX_INFO("ignore wave " << wave << " because outside X_vs_W range " <<  trace.X_vs_W.xmin << " " << trace.X_vs_W.xmax);
	continue;
      }
      if(wave<trace.Y_vs_W.xmin || wave>trace.Y_vs_W.xmax) {
	//SPECEX_INFO("ignore wave " << wave << " because outside Y_vs_W range " <<  trace.X_vs_W.xmin << " " << trace.X_vs_W.xmax);	  
	continue;
      }
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
  is.close();
  
  if (nlines==0) {
    SPECEX_ERROR("loaded " << nlines << " lines");
  }else{
    SPECEX_INFO("loaded " << nlines << " lines");
  }
}


