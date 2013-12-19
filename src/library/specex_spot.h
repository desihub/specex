#ifndef SPECEX_SPOT__H
#define SPECEX_SPOT__H

#include <vector>
#include <string>

#include "harp.hpp"

//#include "spec2dpsf.h"

namespace specex {
  
  class Spot {

    friend class boost::serialization::access;
    
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
    
    
private :

    template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
      ar & BOOST_SERIALIZATION_NVP(log10_wavelength);
      ar & BOOST_SERIALIZATION_NVP(fiber);
      ar & BOOST_SERIALIZATION_NVP(fiber_bundle);
      ar & BOOST_SERIALIZATION_NVP(xc);
      ar & BOOST_SERIALIZATION_NVP(yc);
      ar & BOOST_SERIALIZATION_NVP(flux);
      ar & BOOST_SERIALIZATION_NVP(initial_xc);
      ar & BOOST_SERIALIZATION_NVP(initial_yc);
      ar & BOOST_SERIALIZATION_NVP(initial_flux);
      ar & BOOST_SERIALIZATION_NVP(PSFname);
      ar & BOOST_SERIALIZATION_NVP(PSFParams);
      ar & BOOST_SERIALIZATION_NVP(eflux);     
      ar & BOOST_SERIALIZATION_NVP(fxy_CovMat);
      ar & BOOST_SERIALIZATION_NVP(chi2);
      ar & BOOST_SERIALIZATION_NVP(status);
      
      return;
    }

  };

  BOOST_SERIALIZATION_SHARED_PTR(Spot)  
    typedef boost::shared_ptr < specex::Spot > Spot_p;
  typedef boost::weak_ptr < specex::Spot > Spot_wp;
}


#endif
