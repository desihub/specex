#ifndef SPECEX_GAUSS_HERMITE_ANALYTIC_PSF__H
#define SPECEX_GAUSS_HERMITE_ANALYTIC_PSF__H

#include "specex_psf.h"

#include <string>
#include <vector>


namespace specex {

  class GaussHermitePSF : public PSF {

    friend class boost::serialization::access;
    
  protected :
    int degree;
    
  public :
    
    GaussHermitePSF(int ideg=3);
    virtual ~GaussHermitePSF(){}; 
    
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
    
    bool CheckParams(const harp::vector_double &Params) const 
    { return true;}

    void WriteFits(fitsfile* fp, int first_hdu=1) const;
    void ReadFits(fitsfile* fp, int first_hdu=1);
    
  private :

    void WriteFits_v0(fitsfile* fp, int first_hdu=1) const;
    void ReadFits_v0(fitsfile* fp, int first_hdu=1);
    
    
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
