#ifndef SPECEX_SPECTROGRAPH__H
#define SPECEX_SPECTROGRAPH__H


#include <string>

namespace specex {

  class Spectrograph {
    
  public :
    std::string name;
    int number_of_fiber_bundles_per_ccd;
    int number_of_fibers_per_bundle;
    Spectrograph() {
      name = "NoName";
      number_of_fiber_bundles_per_ccd = 0;
      number_of_fibers_per_bundle = 0;
    }
    void AutoConfigure(const TraceSet& traceset);
  };
  
  class BOSS_Spectrograph : public Spectrograph {
  public :
    BOSS_Spectrograph() {
      name = "BOSS";
      number_of_fiber_bundles_per_ccd = 25;
      number_of_fibers_per_bundle = 20;
    };
  };
  
  class DESI_Spectrograph : public Spectrograph {
  public :
    DESI_Spectrograph() {
      name = "DESI";
      number_of_fiber_bundles_per_ccd = 0; // read in file
      number_of_fibers_per_bundle = 0; // read in file
    };
  };

}


#endif
