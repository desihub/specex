#ifndef SPECEX_LAMP_LINES_UTILS__H
#define SPECEX_LAMP_LINES_UTILS__H

#include <string>
#include <vector>

namespace specex {
  
  class Spot;
  class TraceSet;
  
  void allocate_spots_of_bundle(vector<specex::Spot>& spots, const specex::Spectrograph & spectro, 
				const std::string& lamp_lines_filename, const TraceSet& traceset, 
				int fiber_bundle, int ymin=0, int ymax=10000, 
				const double& min_wavelength=0, const double& max_wavelength=1e6);

};


#endif
