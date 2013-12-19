#ifndef SPECEX_LEGENDRE__H
#define SPECEX_LEGENDRE__H

#include <string>

#include <harp.hpp>

using namespace std;
 
namespace specex {
class Legendre1DPol
{

friend class boost::serialization::access;
 
 public :
  harp::vector_double coeff;
  int deg;
  double xmin,xmax;
  
  Legendre1DPol(int i_deg=0, const double& i_xmin=0, const double& i_xmax=0);
  
  harp::vector_double Monomials(const double &x) const;
  double Value(const double &x) const;
  
  bool Fit(const harp::vector_double& x, const harp::vector_double& y, const harp::vector_double* ey=0, bool set_range = true);
  Legendre1DPol Invert(int add_degree=1) const;
  
  private :

    template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
      ar & BOOST_SERIALIZATION_NVP(coeff);
      ar & BOOST_SERIALIZATION_NVP(deg);
      ar & BOOST_SERIALIZATION_NVP(xmin);
      ar & BOOST_SERIALIZATION_NVP(xmax);
      return;
    }

};
  
Legendre1DPol composed_pol(const Legendre1DPol& pol1, const Legendre1DPol& pol2);
 

class Legendre2DPol
{

  friend class boost::serialization::access;

 public :
  harp::vector_double coeff;
  int xdeg,ydeg;
  double xmin,xmax,ymin,ymax;
  
 Legendre2DPol(int i_xdeg=0, const double& i_xmin=0, const double& i_xmax=0, 
	       int i_ydeg=0, const double& i_ymin=0, const double& i_ymax=0);
  
  harp::vector_double Monomials(const double &x,const double &y) const;
  double Value(const double &x,const double &y) const;
  
  private :

    template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
      ar & BOOST_SERIALIZATION_NVP(coeff);
      ar & BOOST_SERIALIZATION_NVP(xdeg);
      ar & BOOST_SERIALIZATION_NVP(ydeg);
      ar & BOOST_SERIALIZATION_NVP(xmin);
      ar & BOOST_SERIALIZATION_NVP(xmax);
      ar & BOOST_SERIALIZATION_NVP(ymin);
      ar & BOOST_SERIALIZATION_NVP(ymax);
      return;
    }

};
}
  
#endif /* SPECEX_LEGENDRE__H */
