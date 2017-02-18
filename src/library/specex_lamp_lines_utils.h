#ifndef SPECEX_LAMP_LINES_UTILS__H
#define SPECEX_LAMP_LINES_UTILS__H

#include <string>
#include <vector>

namespace specex {
  
  class Spot;
  
  void allocate_spots_of_bundle(vector<specex::Spot_p>& spots,
				const std::string& lamp_lines_filename, const TraceSet& traceset, 
				int fiber_bundle, int fiber_min, int fiber_max, int ymin=0, int ymax=10000, 
				const double& min_wavelength=0, const double& max_wavelength=1e6);

};


#endif
