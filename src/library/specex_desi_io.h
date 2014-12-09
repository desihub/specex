#ifndef SPECEX_BOSS_IO__H
#define SPECEX_BOSS_IO__H

#include <string>
#include <specex_psf.h>

namespace specex {

  class TraceSet;

  
  void read_DESI_traceset_in_fits(
				  TraceSet& traceset,
				  const std::string& arc_filename, 
				  int x_vs_wave_hdu_number, 
				  int y_vs_wave_hdu_number
				  );
};

#endif
