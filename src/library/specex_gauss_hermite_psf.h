#ifndef SPECEX_GAUSS_HERMITE_ANALYTIC_PSF__H
#define SPECEX_GAUSS_HERMITE_ANALYTIC_PSF__H

#include "specex_psf.h"

#include <string>
#include <vector>

//#define ADD_Y_TAILS_TO_GAUSS_HERMITE
//#define ADD_2D_TAILS_TO_GAUSS_HERMITE

// to go faster, less parameters
//#define LORENTZIAN_TAILS

namespace specex {

  class GaussHermitePSF : public PSF {

    friend class boost::serialization::access;
    
  protected :
    int degree;
    
    // buffers to go faster
    harp::vector_double Hx,Hy,dHx,dHy;
    bool need_to_resize_buffer;
    
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
    

    void ResizeBuffer();

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        // serialize base class information
      //ar & BOOST_SERIALIZATION_NVP(boost::serialization::base_object<specex::PSF>(*this));
      ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(PSF);
      ar & BOOST_SERIALIZATION_NVP(degree);
      ar & BOOST_SERIALIZATION_NVP(sigma);

      need_to_resize_buffer = true;
    }
  };
  
  BOOST_SERIALIZATION_SHARED_PTR(GaussHermitePSF)
  
  

  
}



#endif
