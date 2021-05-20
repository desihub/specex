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
  
  void write_spots_list(std::vector<Spot_p> spots, const std::string& filename);
  void write_spots_list(std::vector<Spot_p> spots, const specex::PSF_p psf, const std::string& filename);
  
}

#endif
