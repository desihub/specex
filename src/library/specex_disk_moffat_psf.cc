
#include <cmath>
#include <assert.h>
#include <ctime>

#include <harp.hpp>

#include "specex_hermite.h"
#include "specex_disk_moffat_psf.h"
#include "specex_linalg.h"
#include "specex_message.h"
#include "specex_fits.h"
#include "specex_psf.h"

using namespace std;



specex::DiskMoffatPSF::DiskMoffatPSF(int ideg) : PSF() {
  name = "DiskMoffatPSF";
  SetDegree(ideg);
  
}

//#define DISK_AND_GAUSSIAN

void specex::DiskMoffatPSF::SetDegree(const int ideg) {
  degree = ideg;
}



double specex::DiskMoffatPSF::Profile(const double &X, const double &Y,
			  const harp::vector_double &Params,
			  harp::vector_double *PosDer,
			  harp::vector_double *ParamDer) const
{
  
  const double &R    = Params[0];  // fiber radius in pixel units
  const double & sx  =  Params[1]; // gauss-hermit sigma x
  const double & sy  =  Params[2]; // gauss-hermit sigma x
  int first_hermite_param_index = 3; // 
  
  double sxi = 1./sx;
  double syi = 1./sy;
  double x = X*sxi;
  double y = Y*syi;
  
  
  double s2  = sx*sy;
  double rs2 = (x*x+y*y); // ((X/sx)**2+(Y/sy)**2)
  
  
  // reduced fiber radius
  double R2 = R*R;
  double Rs2 = R2/s2; // = R2/r2*rs2 = R2/(X**2+Y**2)*((X/sx)**2+(Y/sy)**2)
  // internal param
  double As = Rs2-rs2-1;
  double As2 = As*As;
  double disk_moffat_val = 1./(2*M_PI*R2)*(1+As/sqrt(As2+4*Rs2));
  
  double disk_moffat_amp = (1-Params[first_hermite_param_index]); // so total integral = 1
  

  if(ParamDer) {
    double d_dmof_d_Rs2_s = 2*(Rs2+rs2+1)/(2*M_PI*R2)/pow(As2+4*Rs2,3/2.);
    double d_dmof_d_Rs2_r = -disk_moffat_val/Rs2 + d_dmof_d_Rs2_s;
    double d_dmof_d_rs2 = -4*Rs2/(2*M_PI*R2)/pow(As2+4*Rs2,3/2.);
    double d_Rs2_d_R    = 2*R/s2;
    
    double d_rs2_d_sx = 0;
    double d_Rs2_d_sx = 0;
    double d_rs2_d_sy = 0;
    double d_Rs2_d_sy = 0;
    if(rs2>0) {
      d_rs2_d_sx = -2*X*X/(sx*sx*sx);
      d_Rs2_d_sx = -Rs2/sx;
      d_rs2_d_sy = -2*Y*Y/(sy*sy*sy);
      d_Rs2_d_sy = -Rs2/sy;
    }      
    // VALIDATED WITH specex_test_derivatives    
    (*ParamDer)(0) = disk_moffat_amp * d_Rs2_d_R * d_dmof_d_Rs2_r; // dvaldR 
    (*ParamDer)(1) = disk_moffat_amp * ( d_rs2_d_sx * d_dmof_d_rs2 + d_Rs2_d_sx *  d_dmof_d_Rs2_s ); // dvaldsx 
    (*ParamDer)(2) = disk_moffat_amp * ( d_rs2_d_sy * d_dmof_d_rs2 + d_Rs2_d_sy *  d_dmof_d_Rs2_s ); // dvaldsy   
  }
  
  if(PosDer) {
    // do something here
    SPECEX_WARNING("not implemented");
  }
  
  // add gauss-hermite psf
  int nx=(degree+1);
  int ny=(degree+1);
  int nc=nx*ny;
  harp::vector_double Monomials;
  harp::vector_double Monomials_dx;
  harp::vector_double Monomials_dy;
  double prefactor=0;
  
  double expfact=1./(2*M_PI)*sxi*syi*exp(-0.5*(x*x+y*y));
  if(PosDer==0 && ParamDer==0) {
    int param_index=first_hermite_param_index;
    for(int j=0;j<ny;j++) {
      double Hyj=HermitePol(j,y);
      for(int i=0;i<nx;i++,param_index++) {
	prefactor+=Params(param_index)*Hyj*HermitePol(i,x);
      }
    }
  } else if(ParamDer) {
    Monomials.resize(nc);
    Monomials_dx.resize(nc); // need this for derivative of sigma
    Monomials_dy.resize(nc);
    int index=0;
    for(int j=0;j<ny;j++) {
      double Hyj=HermitePol(j,y);
      double dHyj=HermitePolDerivative(j,y);
      for(int i=0;i<nx;i++,index++) {
	double Hxi=HermitePol(i,x);
	Monomials[index]=Hyj*Hxi;
	Monomials_dx[index]=Hyj*HermitePolDerivative(i,x);
	Monomials_dy[index]=dHyj*Hxi;
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
      for(int i=0;i<nx;i++,index++,param_index++) {
	double Hxi=HermitePol(i,x);
	prefactor+=Params(param_index)*Hxi*Hyj;
	Monomials_dx[index]=Hyj*HermitePolDerivative(i,x);
	Monomials_dy[index]=dHyj*Hxi;
      }
    }
  }
  
  double gauss_hermite_val = prefactor*expfact;
  
  if(ParamDer) {

    (*ParamDer)[1] += (x*x-1)*sxi*gauss_hermite_val;
    (*ParamDer)[2] += (y*y-1)*syi*gauss_hermite_val;
    
    (*ParamDer)[1] -= (x*sxi*expfact)*specex::dot(Monomials_dx,ublas::project(Params,ublas::range(first_hermite_param_index,first_hermite_param_index+nc)));
    (*ParamDer)[2] -= (y*syi*expfact)*specex::dot(Monomials_dy,ublas::project(Params,ublas::range(first_hermite_param_index,first_hermite_param_index+nc)));
    
    ublas::project(*ParamDer,ublas::range(first_hermite_param_index,first_hermite_param_index+nc)) = expfact*Monomials;
    
    // add change of disk_moffat_amp
    (*ParamDer)[first_hermite_param_index] -= disk_moffat_val;
  }

  return disk_moffat_amp*disk_moffat_val+gauss_hermite_val;
}

 
int specex::DiskMoffatPSF::LocalNAllPar() const 
{
  int npar = 1; // disk-moffat radius
  npar += 2; // gauss-hermite sigmax sigmay
  npar += (degree+1)*(degree+1); // gauss-hermite
  
#ifdef EXTERNAL_TAIL
  npar += 5;
#endif
  
  return npar;
}
 
harp::vector_double specex::DiskMoffatPSF::DefaultParams() const 
{
  
  harp::vector_double Params(LocalNAllPar());
  Params.clear(); // all = zero at beginning 
  int index=0;
  Params(index++) = 107./3.75*1.70/15./2.; // fiber radius in pixels, tuned for DESI
  Params(index++) = 0.6; // this is sigma_x
  Params(index++) = 0.6; // this is sigma_y  
  
  index += (degree+1)*(degree+1); // gauss-hermite coefficients
  
#ifdef EXTERNAL_TAIL
  Params(index++) = 0.; // tail amplitude
  Params(index++) = 1.; // tail core size
  Params(index++) = 1.; // tail x scale
  Params(index++) = 1.; // tail y scale
  Params(index++) = 2.; // tail power law index
#endif  
  return Params;
}

std::vector<std::string> specex::DiskMoffatPSF::DefaultParamNames() const
{
  std::vector<std::string> paramNames;
  paramNames.push_back("RADIUS");
  paramNames.push_back("GHSIGX");
  paramNames.push_back("GHSIGY");
  char n[10];
  for(int j=0;j<degree+1;j++) {
    for(int i=0;i<degree+1;i++) {
      sprintf(n,"GH-%d-%d",i,j);
      paramNames.push_back(n);
    }
  }
#ifdef EXTERNAL_TAIL
  paramNames.push_back("TAILAMP");
  paramNames.push_back("TAILCORE");
  paramNames.push_back("TAILXSCA");
  paramNames.push_back("TAILYSCA");
  paramNames.push_back("TAILINDE");
#endif  
  return paramNames;
}

void specex::DiskMoffatPSF::Append(const specex::PSF_p other) {
  SPECEX_ERROR("specex::DiskMoffatPSF::Append not implemented");
}
