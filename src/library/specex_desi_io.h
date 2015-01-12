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
				  int y_vs_wave_hdu_number
				  );
  void read_DESI_keywords(
			  const std::string& arc_image_filename, 
			  std::map<std::string,std::string>& infos
			  );
};

#endif
