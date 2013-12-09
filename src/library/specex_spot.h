#ifndef SPECEX_SPOT__H
#define SPECEX_SPOT__H

#include <vector>
#include <string>

#include "harp.hpp"

//#include "spec2dpsf.h"

namespace specex {
  
  class Spot {
    
  public :
    
    //double wavelength; // don't want to do mistakes
    double log10_wavelength;
    int fiber; // fiber id, for which the psf is smoothly changing with xy (but need to account for jumps)
    int fiber_bundle; // id of a list of fibers for which the psf is smoothly changing with xy
    
    double xc; // coordinate of center of spot in CCD
    double yc; // coordinate of center of spot in CCD
    double flux;
    
    double initial_xc; // 
    double initial_yc; // 
    double initial_flux; // 
    
    std::string PSFname;
    harp::vector_double PSFParams;
    
    double eflux; // same info as in xyflux_CovMat(2,2), but might be useful
    harp::matrix_double fxy_CovMat; // covmat of flux, x, y
    harp::matrix_double PSFParams_CovMat;
    harp::matrix_double PSFParams_WeightMat; // inverse of covariance of PSF parameters
    
    harp::vector_double GlobalPSFParams;
    
    double chi2;
    int status;
    Spot() {
      log10_wavelength=0;
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
      
      PSFname="noname";
      
    }
    double wavelength() const {return pow(10,log10_wavelength);}
    
    
    void write_list_header(std::ostream& os) const;
    void write_list_entry(std::ostream& os) const;
    
  };
}


#endif
