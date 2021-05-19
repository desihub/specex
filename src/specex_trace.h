#ifndef SPECEX_TRACE__H
#define SPECEX_TRACE__H



#include <vector>
#include <string>

#include "specex_legendre.h"
#include "specex_spot.h"
#include <boost/serialization/map.hpp>

namespace specex {
  
  

#define SPECEX_TRACE_DEFAULT_LEGENDRE_POL_DEGREE 6
  
  class Trace {
    
    friend class boost::serialization::access;
    
  protected :
    
    
  public :
    
    int fiber; // this trace is for this fiber
    int mask;
    Legendre1DPol X_vs_W; // maps x ccd coordinate as a function of wavelength
    Legendre1DPol Y_vs_W; // maps y ccd coordinate as a function of wavelength
    
    Legendre1DPol W_vs_Y; // as saved is SDSS fits file
    Legendre1DPol X_vs_Y; // as saved is SDSS fits file
    
    // example in /clusterfs/riemann/raid006/bosswork/boss/spectro/redux/current/4097/spArc-r1-00121688.fits.gz[2]
    double yjumplo; // KEY  XJUMPLO in fits table
    double yjumphi; // KEY  XJUMPHI in fits table
    double yjumpval; // KEY  XJUMPVAL in fits table
    
    bool synchronized; // meaning X_vs_W Y_vs_W  W_vs_Y X_vs_Y are consistent
    
    Trace(int i_fiber=-1);
    
    void resize(int nparams);
    
    bool Fit(std::vector<Spot_p> spots, bool set_xy_range = true);
    
    bool Off() const;

    private :

  };
  
  // SDSS IO
  // must be replaced by reading dimension in traceset
#define NUMBER_OF_FIBERS_PER_CCD 500
#define NUMBER_OF_FIBERS_PER_BUNDLE 20
  
  enum TraceSetType {WY = 0, XY = 1, YW = 2, XW = 3};

  typedef std::map<int, Trace> TraceSet;

  int eval_bundle_size(const TraceSet& traceset);

}
#endif
