#ifndef SPECEX_LEGENDRE__H
#define SPECEX_LEGENDRE__H

#include <string>

#include <unbls.h>

using namespace std;
 
namespace specex {
class Legendre1DPol
{

 public :
 string name;
  unbls::vector_double coeff;
  int deg;
  double xmin,xmax;
  
  Legendre1DPol(int i_deg=0, const double& i_xmin=0, const double& i_xmax=0);
  
  unbls::vector_double Monomials(const double &x) const;
  double Value(const double &x) const;
  
  bool Fit(const unbls::vector_double& x, const unbls::vector_double& y, const unbls::vector_double* ey=0, bool set_range = true);
  Legendre1DPol Invert(int add_degree=0) const;

};

typedef std::shared_ptr < specex::Legendre1DPol > Legendre1DPol_p;
  
Legendre1DPol composed_pol(const Legendre1DPol& pol1, const Legendre1DPol& pol2); 

class Legendre2DPol
{

 public :
  string name;
  unbls::vector_double coeff;
  int xdeg,ydeg;
  double xmin,xmax,ymin,ymax;
  
 Legendre2DPol(int i_xdeg=0, const double& i_xmin=0, const double& i_xmax=0, 
	       int i_ydeg=0, const double& i_ymin=0, const double& i_ymax=0);
 
  unbls::vector_double Monomials(const double &x,const double &y) const;
  double Value(const double &x,const double &y) const;
  void Fill();

};

class SparseLegendre2DPol
{

 protected :
  std::vector<int> non_zero_indices;
  
 public :
  string name;
  unbls::vector_double coeff;
  int xdeg,ydeg;
  double xmin,xmax,ymin,ymax;
  int Npar() const { return non_zero_indices.size();}
  void Add(int i,int j);
  void Fill(bool sparse = true); // this is equivalent to a std Legendre2DPol is sparse=false
  void Clear(); // reset

  SparseLegendre2DPol(int i_xdeg=0, const double& i_xmin=0, const double& i_xmax=0, 
		      int i_ydeg=0, const double& i_ymin=0, const double& i_ymax=0);
  
  unbls::vector_double Monomials(const double &x,const double &y) const;
  double Value(const double &x,const double &y) const;
 
};

// this if for convenience in the rest of the code
typedef specex::SparseLegendre2DPol Pol;
typedef std::shared_ptr < specex::SparseLegendre2DPol > Pol_p;

}
  
#endif /* SPECEX_LEGENDRE__H */
