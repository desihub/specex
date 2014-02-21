#ifndef SPECEX_MODEL_IMAGE__H
#define SPECEX_MODEL_IMAGE__H

#include <specex_psf.h>
#include <specex_spot.h>
#include <specex_stamp.h>
#include <specex_image_data.h>

namespace specex {
  
  Stamp compute_stamp(const image_data& model_image, const PSF_p psf, const std::vector<specex::Spot_p>& spots, int x_margin, int y_margin, int only_this_bundle=-1);
			   
  void compute_model_image(image_data& model_image, const specex::image_data& weight, const PSF_p psf, const std::vector<specex::Spot_p>& spots, bool only_on_spots, bool only_psf_core, bool only_positive, int begin_j, int end_j, int x_margin, int y_margin, int only_this_bundle=-1);
  
  void parallelized_compute_model_image(image_data& model_image, const specex::image_data& weight, const PSF_p psf, const std::vector<specex::Spot_p>& spots, bool only_on_spots, bool only_psf_core, bool only_positive, int x_margin, int y_margin, int only_this_bundle=-1);
  
  
};

#endif
