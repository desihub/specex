
#include <cmath>
#include <assert.h>
#include <ctime>

#include <harp.hpp>

#include "specex_hermite.h"
#include "specex_hat_moffat_psf.h"
#include "specex_linalg.h"
#include "specex_message.h"
#include "specex_fits.h"
#include "specex_psf.h"

using namespace std;



specex::HatMoffatPSF::HatMoffatPSF(int ideg) : PSF() {
  name = "HatMoffatPSF";
  SetDegree(ideg);
  
}

//#define HAT_AND_GAUSSIAN

void specex::HatMoffatPSF::SetDegree(const int ideg) {
  degree = ideg;
}



double specex::HatMoffatPSF::Profile(const double &input_X, const double &input_Y,
			  const harp::vector_double &Params,
			  harp::vector_double *PosDer,
			  harp::vector_double *ParamDer) const
{
  /*
    about BOSS fibers diameters in CCD images
    
    * Smee et al. (2013) Section 3.2.3 Slithead :
    "The 120 μm diameter fibers are mounted in
    groups of twenty in v-groove blocks, with 25 v-groove
    blocks being glued to each slithead. The center-to-center
    spacing between fibers on adjacent v-groove blocks is
    624 μm, compared to 260 μm between fibers within a
    v-groove block. The total length of the arc is 138.6 mm
    from outside edge to outside edge of the first and last
    fibers, slightly taller than the original SDSS design."
    
    * Taking into account change of focal :
    - Pixel size = 15 um
    - Measured distance between 2 fibers in b-channel = 6.65 pixels
    - Measured distance between 2 fibers in r-channel = 6.7 pixels (need recheck, should be the same)
    
    - Focal ratio b-channel = (6.65*15)/260
    - Fiber diameter b-channel in pixels = (6.65*15)/260*120./15 = 3.07 pixels
        
  */
  
  const double &fiber_radius = Params[0]; // effective radius, close to 1.5 pix
  const double &moffat_sigma = Params[1]; // sigma of Moffat
  double sigma_x_inv = 1./Params[2]; // gauss-hermit sigma
  double sigma_y_inv = 1./Params[3]; // gauss-hermit sigma
  double x = input_X*sigma_x_inv;
  double y = input_Y*sigma_y_inv;
  int first_hermite_param_index = 4; // first 4 params are 3 sigmas and one radius
  
  // compute hat moffat terms
  double pixel_radius = sqrt(input_X*input_X+input_Y*input_Y);
  
  double r=fiber_radius/moffat_sigma; 
  double u=pixel_radius/moffat_sigma; 
  double r2=r*r;
  double u2=u*u;
  double ru=r*u;
  
  double t1=(r2+u2-2*ru+1)*(r2+u2+2*ru+1);

  double num = atan(r+u)+atan(r-u)+2*r*(1+r2-u2)/t1;
  double denom = M_PI*moffat_sigma*moffat_sigma*( 0.5*M_PI*(r2-1)+2*r+(r2-1)*atan((r2-1)/(2*r))+2*atan(2*r/fabs(r2-1)) );

  double hat_moffat_scale = (1-Params[first_hermite_param_index]); // so that total integral = 1
  double hat_moffat_val   = num/denom;

  if(ParamDer) {
    
    double dnumdr = 4*(1+2*u2+2*r2+6*r2*u2+r2*r2+u2*u2)/(t1*t1);
    double ddendr = 0;
    if(r>1) {
      ddendr = M_PI*moffat_sigma*moffat_sigma*(M_PI*r+2+2*r*atan((r2-1)/(2*r))+2*(r2-3)/(r2+1));
    }else{
      ddendr = M_PI*moffat_sigma*moffat_sigma*(M_PI*r+2+2*r*atan((r2-1)/(2*r))+2);
    }
    
    double dnumdu = -16*ru*(r2+u2+1)/(t1*t1);
    double dvaldr = (dnumdr*denom-num*ddendr)/(denom*denom);
    double dvaldu = dnumdu/denom;
    
    (*ParamDer)(0) = hat_moffat_scale*dvaldr/moffat_sigma; // dvaldR
    (*ParamDer)(1) = -hat_moffat_scale*(dvaldr*r+dvaldu*u+2*hat_moffat_val)/moffat_sigma; // dvaldsigma
    
  }
  if(PosDer) {
    // do something here
    //SPECEX_ERROR("not implemented");
  }
  
  // add gauss-hermite psf
  int nx=(degree+1);
  int ny=(degree+1);
  int nc=nx*ny;
  harp::vector_double Monomials;
  harp::vector_double Monomials_dx;
  harp::vector_double Monomials_dy;
  double prefactor=0;
  
  double expfact=1./(2*M_PI)*sigma_x_inv*sigma_y_inv*exp(-0.5*(x*x+y*y));
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

    (*ParamDer)[2] = (x*x-1)*sigma_x_inv*gauss_hermite_val;
    (*ParamDer)[3] = (y*y-1)*sigma_y_inv*gauss_hermite_val;
    
    (*ParamDer)[2] -= (x*sigma_x_inv*expfact)*specex::dot(Monomials_dx,ublas::project(Params,ublas::range(first_hermite_param_index,first_hermite_param_index+nc)));
    (*ParamDer)[3] -= (y*sigma_y_inv*expfact)*specex::dot(Monomials_dy,ublas::project(Params,ublas::range(first_hermite_param_index,first_hermite_param_index+nc)));
    
    ublas::project(*ParamDer,ublas::range(first_hermite_param_index,first_hermite_param_index+nc)) = expfact*Monomials;
    
    // add change of hat_moffat_scale
    (*ParamDer)[first_hermite_param_index] -= hat_moffat_val;
  }

  return hat_moffat_scale*hat_moffat_val+gauss_hermite_val;
}

 
int specex::HatMoffatPSF::LocalNAllPar() const 
{
  int npar = 2; // hat-moffat radius sigma
  npar += 2; // gauss-hermite sigmax sigmay
  npar += (degree+1)*(degree+1); // gauss-hermite
  
#ifdef EXTERNAL_TAIL
  npar += 5;
#endif
  
  return npar;
}
 
harp::vector_double specex::HatMoffatPSF::DefaultParams() const 
{
  
  harp::vector_double Params(LocalNAllPar());
  Params.clear(); // all = zero at beginning 
  int index=0;
  Params(index++) = 1.5; // fiber radius in pixels
  Params(index++) = 0.5; // moffat sigma
  Params(index++) = 2; // this is sigma_x
  Params(index++) = 2; // this is sigma_y  
  
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

std::vector<std::string> specex::HatMoffatPSF::DefaultParamNames() const
{
  std::vector<std::string> paramNames;
  paramNames.push_back("RADIUS");
  paramNames.push_back("SIGMA");
  paramNames.push_back("GHSIGX2");
  paramNames.push_back("GHSIGY2");
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

void specex::HatMoffatPSF::Append(const specex::PSF_p other) {
  SPECEX_ERROR("specex::HatMoffatPSF::Append not implemented");
}
