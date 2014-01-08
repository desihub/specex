#ifndef SPECEX_BOSS_IO__H
#define SPECEX_BOSS_IO__H

#include <string>
#include <specex_psf.h>

namespace specex {

  class TraceSet;

  std::vector<specex::Legendre1DPol> read_BOSS_singleset_in_fits(const std::string& filename, int hdu_number, bool verbose=false);
  
  
  void read_BOSS_traceset_in_fits(
				  TraceSet& traceset,
				  const std::string& wave_vs_y_filename, int wave_vs_y_hdu_number,
				  const std::string& x_vs_y_filename, int x_vs_y_hdu_number,
				  bool verbose=false
				  );
  void read_BOSS_keywords(PSF_p psf, const std::string& arc_image_filename, std::map<std::string,std::string>& infos);
};

#endif
