#ifndef SPECEX_GAUSS_HERMITE_PSF__H
#define SPECEX_GAUSS_HERMITE_PSF__H

#include "specex_psf.h"

#include <string>
#include <vector>


namespace specex {

  class GaussHermitePSF : public PSF {

    friend class boost::serialization::access;
    
  protected :
    int degree;
    
  public :

    typedef std::shared_ptr <GaussHermitePSF> pshr;
    
    GaussHermitePSF(int ideg=3);
    virtual ~GaussHermitePSF(){}; 
    
    void SetDegree(const int ideg);
  
    int LocalNAllPar() const;

    double Degree() const {
      return degree;
    }
    
    double Profile(const double &X, const double &Y,
		   const unhrp::vector_double &Params,
		   unhrp::vector_double *PosDer = 0,
		   unhrp::vector_double *ParamGradient = 0) const;
    
    // needed for analytic integration
    double PixValue(const double &Xc, const double &Yc,
				     const double &XPix, const double &YPix,
				     const unhrp::vector_double &Params,
				     unhrp::vector_double *PosDer,
				 unhrp::vector_double *ParamDer) const;
    
    unhrp::vector_double DefaultParams() const;
    std::vector<std::string> DefaultParamNames() const;
    
    bool CheckParams(const unhrp::vector_double &Params) const 
    { return true;}
    
    void Append(const specex::PSF_p other);
    
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(PSF);
      ar & BOOST_SERIALIZATION_NVP(degree);
    }
  };
  
  BOOST_SERIALIZATION_SHARED_PTR(GaussHermitePSF)
    
}



#endif
