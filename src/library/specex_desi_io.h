#ifndef SPECEX_BOSS_IO__H
#define SPECEX_BOSS_IO__H

#include <string>
#include <specex_psf.h>

namespace specex {

  class TraceSet;
  
  void read_DESI_traceset_in_fits(
				  TraceSet& traceset,
				  const std::string& x_vs_wave_filename, 
				  int x_vs_wave_hdu_number, 
				  const std::string& y_vs_wave_filename, 
				  int y_vs_wave_hdu_number,
				  int required_x_vs_wave_degree = 0,
				  int required_y_vs_wave_degree = 0
				  );
  void read_DESI_preprocessed_image(const std::string& arc_image_filename, image_data &image, image_data &weight, image_data &mask, image_data &read_noise, std::map<std::string,std::string>& header);
};

#endif
