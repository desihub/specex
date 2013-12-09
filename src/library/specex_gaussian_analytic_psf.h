#ifndef SPECEX_GAUSSIAN_ANALYTIC_PSF__H
#define SPECEX_GAUSSIAN_ANALYTIC_PSF__H

#include "specex_base_analytic_psf.h"

namespace specex {

//! gaussian PSF with 3 parameters :wxx, wyy , wxy.
  class GaussPSF : public AnalyticPSF {
  
  public :
    
    GaussPSF();
    virtual ~GaussPSF(){};
    
    std::string Name() const;
    
    size_t NPar() const;
    
    double Profile(const double &X, const double &Y,
		   const harp::vector_double &Params,
		   harp::vector_double *PosDer = 0,
		   harp::vector_double *ParamGradient = 0) const;
    
    void InitParams(const double &sigma, harp::vector_double &Params);
    
    bool CheckParams(const harp::vector_double &Params) const;


  

  };
}
  
#endif
