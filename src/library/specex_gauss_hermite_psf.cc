
#include <cmath>
#include <assert.h>
#include <ctime>



#include "harp.hpp"

#include "specex_hermite.h"
#include "specex_gauss_hermite_psf.h"
#include "specex_linalg.h"
#include "specex_message.h"
#include "specex_fits.h"

using namespace std;



specex::GaussHermitePSF::GaussHermitePSF(int ideg) {
  name = "GaussHermitePSF";
  SetDegree(ideg);
  
}

#define TWO_GAUSSIANS

void specex::GaussHermitePSF::SetDegree(const int ideg) {
  degree = ideg;
}

int specex::GaussHermitePSF::LocalNAllPar() const {
    
  int npar = 2; // sigma_x and sigma_y
  
#ifdef TWO_GAUSSIANS
  npar += 3;
#endif

  // npar += (degree+1)*(degree+1)-4;// skip (0,0)(1,0)(0,1)(1,1) : -1 because normalized, -3 because centered 
  npar += (degree+1)*(degree+1)-1;// skip (0,0) : -1 because normalized
  

  return npar;
}

double specex::GaussHermitePSF::Profile(const double &input_X, const double &input_Y,
			  const harp::vector_double &Params,
			  harp::vector_double *PosDer,
			  harp::vector_double *ParamDer) const
{

  double sigma_x_inv = 1./Params(0);
  double sigma_y_inv = 1./Params(1);
  double x = input_X*sigma_x_inv;
  double y = input_Y*sigma_y_inv;
#ifdef TWO_GAUSSIANS
  double sigma_x_2_inv = 1./Params(2);
  double sigma_y_2_inv = 1./Params(3);
  double scale_2       = Params(4);
  double scale_1       = 1-scale_2;
  double x_2 = input_X*sigma_x_2_inv;
  double y_2 = input_Y*sigma_y_2_inv;
#endif
  

 
  int nx=(degree+1);
  int ny=(degree+1);
  //int nc=nx*ny-4; // skip (0,0)(1,0)(0,1)(1,1)
  int nc=nx*ny-1; // skip (0,0)
  
  // precompute to go faster
  harp::vector_double Monomials;
  harp::vector_double Monomials_dx;
  harp::vector_double Monomials_dy;
  double prefactor=1;
#ifndef TWO_GAUSSIANS
  int first_hermite_param_index = 2; // first 2 params are sigmas
  double expfact=1./(2*M_PI)*sigma_x_inv*sigma_y_inv*exp(-0.5*(x*x+y*y));
#else
  int first_hermite_param_index = 5; // first 5 params are sigmas
  double expfact_1_part=1/(2*M_PI)*sigma_x_inv*sigma_y_inv*exp(-0.5*(x*x+y*y));
  double expfact_2_part=1/(2*M_PI)*sigma_x_2_inv*sigma_y_2_inv*exp(-0.5*(x_2*x_2+y_2*y_2));
  double expfact_1=scale_1*expfact_1_part;
  double expfact_2=scale_2*expfact_2_part;
  double expfact=expfact_1+expfact_2;
#endif


  if(PosDer==0 && ParamDer==0) {
    int param_index=first_hermite_param_index;
    for(int j=0;j<ny;j++) {
      double Hyj=HermitePol(j,y);
      //int imin=0; if(j<2) {imin=2;} // skip (0,0)(1,0)(0,1)(1,1)
      int imin=0; if(j==0) imin=1; // skip (0,0)
      for(int i=imin;i<nx;i++,param_index++) {
	prefactor+=Params(param_index)*Hyj*HermitePol(i,x);
      }
    }
    return expfact*prefactor;
    
  } else if(PosDer==0 && ParamDer) {
    Monomials.resize(nc);
    int index=0;
    for(int j=0;j<ny;j++) {
      double Hyj=HermitePol(j,y);
      //int imin=0; if(j<2) {imin=2;} // skip (0,0)(1,0)(0,1)(1,1)
      int imin=0; if(j==0) imin=1; // skip (0,0)
      for(int i=imin;i<nx;i++,index++) {
	Monomials[index]=Hyj*HermitePol(i,x);
      }
    }
    prefactor += specex::dot(Monomials,ublas::project(Params,ublas::range(first_hermite_param_index,first_hermite_param_index+nc)));
  
  } else if(PosDer && ParamDer==0) {
    Monomials_dx.resize(nc);
    Monomials_dy.resize(nc);
    int param_index=first_hermite_param_index;
    int index=0;
    for(int j=0;j<ny;j++) {
      double Hyj=HermitePol(j,y);
      double dHyj=HermitePolDerivative(j,y);
      
      //int imin=0; if(j<2) {imin=2;} // skip (0,0)(1,0)(0,1)(1,1)
      int imin=0; if(j==0) imin=1; // skip (0,0)
      for(int i=imin;i<nx;i++,index++,param_index++) {
	
	double Hxi=HermitePol(i,x);
	prefactor+=Params(param_index)*Hxi*Hyj;
	Monomials_dx[index]=Hyj*HermitePolDerivative(i,x);
	Monomials_dy[index]=dHyj*Hxi;
      }
    }
  } else if(PosDer && ParamDer) {
    Monomials.resize(nc);
    Monomials_dx.resize(nc);
    Monomials_dy.resize(nc);
    int index=0;
    for(int j=0;j<ny;j++) {
      double Hyj=HermitePol(j,y);
      double dHyj=HermitePolDerivative(j,y);
      
      //int imin=0; if(j<2) {imin=2;} // skip (0,0)(1,0)(0,1)(1,1)
      int imin=0; if(j==0) imin=1; // skip (0,0)
      for(int i=imin;i<nx;i++,index++) {
	double Hxi=HermitePol(i,x);
	Monomials[index]=Hyj*Hxi;
	Monomials_dx[index]=Hyj*HermitePolDerivative(i,x);
	Monomials_dy[index]=dHyj*Hxi;
      }
    }
    prefactor += specex::dot(Monomials,ublas::project(Params,ublas::range(first_hermite_param_index,first_hermite_param_index+nc)));
  }
  
  
  double psf_val = expfact*prefactor;
  
  if(ParamDer) {
#ifndef TWO_GAUSSIANS
    (*ParamDer)[0] = (x*x-1)*sigma_x_inv*psf_val; // exact ONLY if all gauss-hermite terms except zeroth order  = 0
    (*ParamDer)[1] = (y*y-1)*sigma_y_inv*psf_val; // exact ONLY if all gauss-hermite terms except zeroth order  = 0
#else
    double psf_val_1 = expfact_1*prefactor;
    double psf_val_2 = expfact_2*prefactor;
    
    (*ParamDer)[0] = (x*x-1)*sigma_x_inv*psf_val_1;
    (*ParamDer)[1] = (y*y-1)*sigma_y_inv*psf_val_1;
    (*ParamDer)[2] = (x_2*x_2-1)*sigma_x_2_inv*psf_val_2;
    (*ParamDer)[3] = (y_2*y_2-1)*sigma_y_2_inv*psf_val_2;
    (*ParamDer)[4] = (-expfact_1_part+expfact_2_part)*prefactor;
#endif
    ublas::project(*ParamDer,ublas::range(first_hermite_param_index,first_hermite_param_index+nc)) = expfact*Monomials;
  }
  
  if(PosDer) { // wrong sign on purpose (derivatives w.r.t -X)
    double& dvdx=(*PosDer)(0);
    double& dvdy=(*PosDer)(1);  
   
#ifndef TWO_GAUSSIANS 
    dvdx=x*sigma_x_inv*psf_val;
    dvdy=y*sigma_y_inv*psf_val;
    
    double d_poly_dx = specex::dot(Monomials_dx,ublas::project(Params,ublas::range(first_hermite_param_index,first_hermite_param_index+nc)));
    double d_poly_dy = specex::dot(Monomials_dy,ublas::project(Params,ublas::range(first_hermite_param_index,first_hermite_param_index+nc)));
    dvdx -= d_poly_dx*expfact*sigma_x_inv; // minus sign cause derivative wrt -x
    dvdy -= d_poly_dy*expfact*sigma_y_inv;
#else 
    
    double psf_val_1 = expfact_1*prefactor;
    double psf_val_2 = expfact_2*prefactor;
    
    dvdx=x*sigma_x_inv*psf_val_1;
    dvdy=y*sigma_y_inv*psf_val_1;
    dvdx += x_2*sigma_x_2_inv*psf_val_2;
    dvdy += y_2*sigma_y_2_inv*psf_val_2;
    
    double d_poly_dx = specex::dot(Monomials_dx,ublas::project(Params,ublas::range(first_hermite_param_index,first_hermite_param_index+nc)));
    double d_poly_dy = specex::dot(Monomials_dy,ublas::project(Params,ublas::range(first_hermite_param_index,first_hermite_param_index+nc)));
    dvdx -= d_poly_dx*expfact*sigma_x_inv; // minus sign cause derivative wrt -x
    dvdy -= d_poly_dy*expfact*sigma_y_inv;



#endif
}
  
  return psf_val;
}

  
harp::vector_double specex::GaussHermitePSF::DefaultParams() const 
{
  
  harp::vector_double Params(LocalNAllPar());
  Params.clear(); // all = zero at beginning = a pure gaussian
  Params(0) = 1.1; // this is sigma_x
  Params(1) = 1.1; // this is sigma_y
#ifdef TWO_GAUSSIANS
  Params(2) = 3.; // this is sigma_x_2
  Params(3) = 3.; // this is sigma_y_2
  Params(4) = 0.1; // this is the amplitude of the second gaussian
#endif
  return Params;
}

std::vector<std::string> specex::GaussHermitePSF::DefaultParamNames() const
{
  std::vector<std::string> paramNames;
  paramNames.push_back("GHSIGX");
  paramNames.push_back("GHSIGY");
    
#ifdef TWO_GAUSSIANS
  paramNames.push_back("GHSIGX2");
  paramNames.push_back("GHSIGY2");
  paramNames.push_back("GHSCAL2");
#endif 
    
  char n[10];
  for(int j=0;j<degree+1;j++) {
    for(int i=0;i<degree+1;i++) {
      if(i==0 && j==0) continue;
      //if(i==1 && j==0) continue;
      //if(i==0 && j==1) continue;
      //if(i==1 && j==1) continue;
      
      sprintf(n,"GH-%d-%d",i,j);
      paramNames.push_back(n);
    }
  }
  return paramNames;
}


