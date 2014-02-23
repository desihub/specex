#ifndef SPECEX_HAT_MOFFAT_PSF__H
#define SPECEX_HAT_MOFFAT_PSF__H

#include "specex_psf.h"

#include <string>
#include <vector>


namespace specex {

  class HatMoffatPSF : public PSF {

    friend class boost::serialization::access;
    
  protected :
    int degree;
    
  public :
    
    HatMoffatPSF(int ideg=3);
    virtual ~HatMoffatPSF(){}; 
    
    void SetDegree(const int ideg);
  
    int LocalNAllPar() const;
    
    double Degree() const {
      return degree;
    }
    
    double Profile(const double &X, const double &Y,
		   const harp::vector_double &Params,
		   harp::vector_double *PosDer = 0,
		   harp::vector_double *ParamGradient = 0) const;
    
    harp::vector_double DefaultParams() const;
    std::vector<std::string> DefaultParamNames() const;
    
    bool CheckParams(const harp::vector_double &Params) const 
    { return true;}

    void Append(const specex::PSF_p other);

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(PSF);
      ar & BOOST_SERIALIZATION_NVP(degree);
    }
  };
  
  BOOST_SERIALIZATION_SHARED_PTR(HatMoffatPSF)
    
}



#endif
