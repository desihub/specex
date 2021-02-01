#ifndef SPECEX_SPOT_ARRAY__H
#define SPECEX_SPOT_ARRAY__H

#include <vector>
#include <string>

#include "specex_spot.h"
#include "specex_psf.h"

namespace specex {
  class SpotArray : public std::vector<Spot_p> {
  public :
    double wavelength;
    int fiber_bundle;
    SpotArray() {
      wavelength=0;
      fiber_bundle=-1;
    }
  };
  
  std::vector<SpotArray> find_spot_arrays( std::vector<Spot_p> &spots, bool check_status=false);
  std::vector<SpotArray> isolated_spot_arrays(  std::vector<SpotArray>& other, const double& delta_lambda=10 );
  std::vector<SpotArray> valid_spot_arrays(  std::vector<SpotArray>& other, int required_status=1, const double& min_valid_fraction=0.5);
  SpotArray         merge_spot_arrays(  std::vector<SpotArray>& other);
  
  void write_spots_list(std::vector<Spot_p> spots, const std::string& filename);
  void write_spots_list(std::vector<Spot_p> spots, const specex::PSF_p psf, const std::string& filename);
  
}

#endif
