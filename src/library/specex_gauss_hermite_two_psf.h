#ifndef SPECEX_GAUSS_HERMITE_TWO_PSF__H
#define SPECEX_GAUSS_HERMITE_TWO_PSF__H

#include "specex_psf.h"

#include <string>
#include <vector>


namespace specex {

  class GaussHermite2PSF : public PSF {

    friend class boost::serialization::access;
    
  public :
    int core_degree, second_degree;
    
  public :
    
    GaussHermite2PSF(int i_core_degree=3, int i_second_degree=3);
    virtual ~GaussHermite2PSF(){}; 
    
    void SetDegree(const int i_core_degree, const int i_second_degree);
  
    int LocalNAllPar() const;
    
    double Profile(const double &X, const double &Y,
		   const harp::vector_double &Params,
		   harp::vector_double *PosDer = 0,
		   harp::vector_double *ParamGradient = 0) const;
    
    harp::vector_double DefaultParams() const;
    std::vector<std::string> DefaultParamNames() const;
    
    bool CheckParams(const harp::vector_double &Params) const 
    { return true;}
    
    void Append(const specex::PSF_p other);
    
    // will be optimized. for now a hack for the tail at the boundary of core radius
    virtual double PSFValueWithParamsXY(const double& X, const double &Y, 
					const int IPix, const int JPix,
					const harp::vector_double &Params,
					harp::vector_double *PosDer, harp::vector_double *ParamDer,
					bool with_core=true, bool with_tail=true) const;
    

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(PSF);
      ar & BOOST_SERIALIZATION_NVP(core_degree);
      ar & BOOST_SERIALIZATION_NVP(second_degree);
    }
  };
  
  BOOST_SERIALIZATION_SHARED_PTR(GaussHermite2PSF)
    
}



#endif
