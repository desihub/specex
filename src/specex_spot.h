#ifndef SPECEX_SPOT__H
#define SPECEX_SPOT__H

#include <vector>
#include <string>
#include <memory>

namespace specex {
  
  class Spot {

//    friend class boost::serialization::access;
    
  public :
    
    double wavelength; 
    
    int fiber; // fiber id, for which the psf is smoothly changing with xy (but need to account for jumps)
    int fiber_bundle; // id of a list of fibers for which the psf is smoothly changing with xy
    
    double xc; // coordinate of center of spot in CCD
    double yc; // coordinate of center of spot in CCD
    double flux;
    
    double initial_xc; // 
    double initial_yc; // 
    double initial_flux; // 
    
    double eflux; // 
        
    double chi2;
    int status;
    Spot() {
      wavelength=0;
      fiber=0;
      fiber_bundle=0;
      xc=0;
      yc=0;
      flux=0;
      chi2=1e20;
      status=0; //1 for successful fit
      
      initial_xc=0;
      initial_yc=0;
      initial_flux=0;
      
    }
    
    
    void write_list_header(std::ostream& os) const;
    void write_list_entry(std::ostream& os) const;
    
    
private :

  };

  typedef std::shared_ptr < specex::Spot > Spot_p;
  typedef std::weak_ptr   < specex::Spot > Spot_wp;
}


#endif
