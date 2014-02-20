#ifndef SPECEX_MODEL_IMAGE__H
#define SPECEX_MODEL_IMAGE__H

#include <specex_psf.h>
#include <specex_spot.h>
#include <specex_stamp.h>
#include <specex_image_data.h>

namespace specex {
  
  Stamp compute_stamp(const image_data& model_image, const PSF_p psf, const std::vector<specex::Spot_p>& spots);
			   
  void compute_model_image(image_data& model_image, const specex::image_data& weight, const PSF_p psf, const std::vector<specex::Spot_p>& spots, const int first_fiber, const int last_fiber, bool only_on_spots, bool only_psf_core, bool only_positive);
  
  
};

#endif
