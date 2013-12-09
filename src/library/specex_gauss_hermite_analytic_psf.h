#ifndef SPECEX_GAUSS_HERMITE_ANALYTIC_PSF__H
#define SPECEX_GAUSS_HERMITE_ANALYTIC_PSF__H

#include "specex_base_analytic_psf.h"

#include <string>
#include <vector>

//#define ADD_Y_TAILS_TO_GAUSS_HERMITE
//#define ADD_2D_TAILS_TO_GAUSS_HERMITE

// to go faster, less parameters
//#define LORENTZIAN_TAILS

namespace specex {

  class GaussHermitePSF : public AnalyticPSF {
    
  protected :
    int degree;
    
    // buffers to go faster
    harp::vector_double Hx,Hy,dHx,dHy;

    int tail_norm_index; 

#ifndef LORENTZIAN_TAILS
    int tail_power_index;
    int tail_x_scale_plus_index;
    int tail_x_scale_minus_index;
    int tail_y_scale_minus_index;
#endif  
    
    int y_tail_norm_index; 
    
  public :
  
    double sigma; 
    
    GaussHermitePSF(int ideg=3);
    virtual ~GaussHermitePSF(){}; 
    
    void SetDegree(const int ideg);
  
    size_t NPar() const;
    
    double Degree() const {
      return degree;
    }
    
    std::string Name() const {return "GAUSSHERMITE";}
    
    
    
    double Profile(const double &X, const double &Y,
		   const harp::vector_double &Params,
		   harp::vector_double *PosDer = 0,
		   harp::vector_double *ParamGradient = 0) const;
    void InitParams(const double &sigma, harp::vector_double &Params);
    
    
    bool CheckParams(const harp::vector_double &Params) const 
    { return true;}
	
    
  };
  
}
#endif
