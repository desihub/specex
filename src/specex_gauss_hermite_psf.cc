
#include <cmath>
#include <assert.h>
#include <ctime>

#include <harp.hpp>

#include "specex_hermite.h"
#include "specex_gauss_hermite_psf.h"
#include "specex_linalg.h"
#include "specex_message.h"
#include "specex_fits.h"
#include "specex_psf.h"

using namespace std;



specex::GaussHermitePSF::GaussHermitePSF(int ideg) : PSF() 
{
  name = "GaussHermitePSF";
  SetDegree(ideg);
  
}



void specex::GaussHermitePSF::SetDegree(const int ideg) {
  degree = ideg;
}



double specex::GaussHermitePSF::Profile(const double &input_X, const double &input_Y,
			  const harp::vector_double &Params,
			  harp::vector_double *PosDer,
			  harp::vector_double *ParamDer) const
{

  double sx = Params(0);
  double sy = Params(1);
  if(sx<0.1) {sx=0.1;} // to avoid failures in exploration of model params
  if(sy<0.1) {sy=0.1;} // to avoid failures in exploration of model params
  double sigma_x_inv = 1./sx;
  double sigma_y_inv = 1./sy;
  double x = input_X*sigma_x_inv;
  double y = input_Y*sigma_y_inv;

 
  int nx=(degree+1);
  int ny=(degree+1);
  int nc=nx*ny-1; // skip (0,0)
  
  // precompute to go faster
  harp::vector_double Monomials;
  harp::vector_double Monomials_dx;
  harp::vector_double Monomials_dy;
  double prefactor=1;
  int first_hermite_param_index = 2; // first 2 params are sigmas
  double expfact=1./(2*M_PI)*sigma_x_inv*sigma_y_inv*exp(-0.5*(x*x+y*y));


  if(PosDer==0 && ParamDer==0) {
    int param_index=first_hermite_param_index;
    for(int j=0;j<ny;j++) {
      double Hyj=HermitePol(j,y);
      int imin=0; if(j==0) imin=1; // skip (0,0)
      for(int i=imin;i<nx;i++,param_index++) {
	prefactor+=Params(param_index)*Hyj*HermitePol(i,x);
      }
    }
    return expfact*prefactor;
    
  } else if(ParamDer) {
    Monomials.resize(nc);
    Monomials_dx.resize(nc); // need this for derivative of sigma
    Monomials_dy.resize(nc);
    int index=0;
    for(int j=0;j<ny;j++) {
      double Hyj=HermitePol(j,y);
      double dHyj=HermitePolDerivative(j,y);
      int imin=0; if(j==0) imin=1; // skip (0,0)
      for(int i=imin;i<nx;i++,index++) {
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
      int imin=0; if(j==0) imin=1; // skip (0,0)
      for(int i=imin;i<nx;i++,index++,param_index++) {
	
	double Hxi=HermitePol(i,x);
	prefactor+=Params(param_index)*Hxi*Hyj;
	Monomials_dx[index]=Hyj*HermitePolDerivative(i,x);
	Monomials_dy[index]=dHyj*Hxi;
      }
    }
  }
  
  
  double psf_val = expfact*prefactor;
  
  if(ParamDer) {

    (*ParamDer)[0] = (x*x-1)*sigma_x_inv*psf_val; // exact ONLY if all gauss-hermite terms except zeroth order  = 0
    (*ParamDer)[1] = (y*y-1)*sigma_y_inv*psf_val; // exact ONLY if all gauss-hermite terms except zeroth order  = 0
    
    /* other terms 

       dPSD/sigmax += sum_ij Pij dH(x,i)/dsigmax H(y,j) * expfact
                   += dx/dsigmax * sum_ij Pij dH(x,i)dx *H(y,j) * expfact
                   -= x/sigmax * sum_i P_i Monomials_dx_i * expfact
     */
     
    (*ParamDer)[0] -= (x*sigma_x_inv*expfact)*specex::dot(Monomials_dx,ublas::project(Params,ublas::range(first_hermite_param_index,first_hermite_param_index+nc)));
    (*ParamDer)[1] -= (y*sigma_y_inv*expfact)*specex::dot(Monomials_dy,ublas::project(Params,ublas::range(first_hermite_param_index,first_hermite_param_index+nc)));
    
    ublas::project(*ParamDer,ublas::range(first_hermite_param_index,first_hermite_param_index+nc)) = expfact*Monomials;
  }
  
  if(PosDer) { // wrong sign on purpose (derivatives w.r.t -X)
    double& dvdx=(*PosDer)(0);
    double& dvdy=(*PosDer)(1);  
   

    dvdx=x*sigma_x_inv*psf_val;
    dvdy=y*sigma_y_inv*psf_val;
    
    double d_poly_dx = specex::dot(Monomials_dx,ublas::project(Params,ublas::range(first_hermite_param_index,first_hermite_param_index+nc)));
    double d_poly_dy = specex::dot(Monomials_dy,ublas::project(Params,ublas::range(first_hermite_param_index,first_hermite_param_index+nc)));
    dvdx -= d_poly_dx*expfact*sigma_x_inv; // minus sign cause derivative wrt -x
    dvdy -= d_poly_dy*expfact*sigma_y_inv;
}
  
  return psf_val;
}

double specex::GaussHermitePSF::PixValue(const double &Xc, const double &Yc,
					 const double &XPix, const double &YPix,
					 const harp::vector_double &Params,
					 harp::vector_double *PosDer,
					 harp::vector_double *ParamDer) const {

  
  // shut down analytic calculation
  // return specex::PSF::PixValue(Xc,Yc,XPix,YPix,Params,PosDer,ParamDer);
  
  

  
  /*
    P_n(x) = 1./sqrt(2*pi)*exp(-x**2/2)*H_n(x)
    P_n(x) = -P_{n-1}'(x) for n>0, 
    can be demonstrated with the recursive relation of P_n(P_(n-1),P_(n-2))
    
    int P_n(x) dx = P_{n-1}(x-hx)-P_{n-1}(x+hx)
    
    For n=0, 
    int P_0(x) dx = e(x) = 0.5*(erf((x+hx)/sq2)-erf((x-hx)/sq2))
    
    int dx dy P(x,y)
    = int dx dy sum_ij cij Pi(x) Pj(y) + A 
    = sum_ij cij [int dx Pi(x)][int dy Pj(y)] + A
    = sum_ij cij [P(i-1)(x-hx)-P(i-1)(x+hx)][P(j-1)(y-hy)-P(j-1)(y+hy)] + A
    = sum_ij c(i+1,j+1) [Pi(x-hx)-Pi(x+hx)][Pj(y-hy)-Pj(y+hy)] + A
    = prof(c+,x-hx,y-hy)-prof(c+,x+hx,y-hy)-prof(c+,x-hx,y+hy)+prof(c+,x+hx,y+hy) + A
    
    A = e(x) * sum_j c(0,j+1) [Pj(y-hy)-Pj(y+hy)] + e(y) * sum_i c(i+1,0) [Py(x-hx)-Pi(x+hx)] + c(0,0)*e(x)*e(y)
  */

  
  // sigmas of Gaussian
  double sx = Params(0);
  double sy = Params(1);
  if(sx<0.1) {sx=0.1;} // to avoid failures in exploration of model params
  if(sy<0.1) {sy=0.1;} // to avoid failures in exploration of model params
  
  const double isx = 1./sx;
  const double isy = 1./sy;
  
  // reduced coordinates of edges of pixel
  const double x1 = (floor(XPix+0.5) - Xc - 0.5)*isx;
  const double x2 = (floor(XPix+0.5) - Xc + 0.5)*isx;
  const double y1 = (floor(YPix+0.5) - Yc - 0.5)*isy;
  const double y2 = (floor(YPix+0.5) - Yc + 0.5)*isy;
  
  const double isq2 = 1./sqrt(2.);
  const double isq2pi = 1./sqrt(2.*M_PI);

  // gaussians values on edges of pixel
  const double gx1=isq2pi*isx*exp(-0.5*(x1*x1));
  const double gx2=isq2pi*isx*exp(-0.5*(x2*x2));
  const double gy1=isq2pi*isy*exp(-0.5*(y1*y1));
  const double gy2=isq2pi*isy*exp(-0.5*(y2*y2));

  // integral of gaussians in pixel  
  const double ex = 0.5*( erf(x2*isq2)-erf(x1*isq2) );
  const double ey = 0.5*( erf(y2*isq2)-erf(y1*isq2) );
  
  
  
  // indexation
  const int first_hermite_param_index = 2; // first 2 params are sigmas
  const int nx=(degree+1);
  const int ny=(degree+1);
  const int nc=nx*ny-1; // skip (0,0)

  // tmp variables
  double Pi=0;
  double Pj=0;
  
  double psfval=0;

  if(PosDer==0 && ParamDer==0) { // no computation of derivatives, just value
    
    int param_index=first_hermite_param_index;
    
    psfval = ex*ey; // 0th order
    
    for(int j=0;j<ny;j++) {
      
      
      if(j==0)
	Pj=ey;
      else 
	Pj=sy*( gy1*HermitePol(j-1,y1)-gy2*HermitePol(j-1,y2) );

      int imin=0; if(j==0) imin=1; // skip (0,0)
      for(int i=imin;i<nx;i++,param_index++) {
	
	if(i==0)
	  Pi=ex;
	else 
	  Pi=sx*( gx1*HermitePol(i-1,x1)-gx2*HermitePol(i-1,x2) );
	
	psfval+=Params(param_index)*Pj*Pi;
      }
    }
    
    return psfval;
    
  } 

  
  // derivative of integral of gaussians in pixel wrt sigmax and sigmay
  const double dexdsx = (x1 * gx1 - x2 * gx2); // VALIDATED
  const double deydsy = (y1 * gy1 - y2 * gy2); // VALIDATED
  
  // derivative of integral of gaussians in pixel wrt xc and yc ( - that of x or y)
  const double dexdx = (gx1 -gx2); // VALIDATED
  const double deydy = (gy1 -gy2); // VALIDATED
  

  // full derivative computation
  // PiPj = int_{x1,y1}^{x2,y2} dx dy exp(-(x**2+y**2)/2) * Hi(x) * Hj(y)
  harp::vector_double PiPj(nc);    
  harp::vector_double dPiPjdsx(nc);// derivative wrt sigma
  harp::vector_double dPiPjdsy(nc);
  harp::vector_double dPiPjdx(nc);// derivative wrt x
  harp::vector_double dPiPjdy(nc);

  /*
  PiPj.resize(nc);
  dPiPjdsx.resize(nc); // derivative wrt sigma
  dPiPjdsy.resize(nc);
  dPiPjdx.resize(nc); // derivative wrt x
  dPiPjdy.resize(nc);
  */
 
 
  
  
    double dPidsx=0;
    double dPjdsy=0;
    double dPidx=0;
    double dPjdy=0;
    double t1=0;
    double t2=0;
    int index=0;
    for(int j=0;j<ny;j++) {
      
      if(j==0) {
	Pj=ey;// VALIDATED
	dPjdsy=deydsy;// VALIDATED
	dPjdy=deydy;
	
      }
      else {
	t1=sy*gy1*HermitePol(j-1,y1);
	t2=sy*gy2*HermitePol(j-1,y2);
	Pj= t1 - t2;
	dPjdsy = ( -gy1*y1*HermitePolDerivative(j-1,y1)+gy2*y2*HermitePolDerivative(j-1,y2) )
	  + ( t1*y1*y1*isy - t2*y2*y2*isy ); // VALIDATED	
	dPjdy = -isy*( sy*gy1*HermitePolDerivative(j-1,y1)-sy*gy2*HermitePolDerivative(j-1,y2) - y1*t1 + y2*t2  ); // minus sign because derivatice wrt yc  , VALIDATED 
      }
      
      
      int imin=0; if(j==0) imin=1; // skip (0,0)
      for(int i=imin;i<nx;i++,index++) {
	
	if(i==0) {
	  Pi=ex;// VALIDATED
	  dPidsx=dexdsx;// VALIDATED
	  dPidx=dexdx;
	}
	else {
	  t1=sx*gx1*HermitePol(i-1,x1); 
	  t2=sx*gx2*HermitePol(i-1,x2);
	  Pi= t1 - t2 ; // VALIDATED
	  // t = sx*gx*H(x) = sx*isq2pi*isx*exp(-0.5*x**2)*H(x) 
	  // t = a*g(x)*H(x)
	  // dH/ds = dH/dx * dx/ds = dH/dx * -x/s
	  // dg/ds = dg/dx * dx/s  = -x*g * -x/s = x**2/s * g
	  // dt/ds = a*g*dH/ds + a*dg/ds*H = a*g * dH/dx * -x/s + a*g*H*x**2/s 
	  dPidsx = ( -gx1*x1*HermitePolDerivative(i-1,x1)+gx2*x2*HermitePolDerivative(i-1,x2) )
	    + ( t1*x1*x1*isx - t2*x2*x2*isx ); // VALIDATED

	  // dt/dx = a*g*dH/dx + a*dg/dx*H = a*g * dH/dx - a*g*H*x 	  
	  dPidx = -isx*( sx*gx1*HermitePolDerivative(i-1,x1)-sx*gx2*HermitePolDerivative(i-1,x2) - x1*t1 + x2*t2  ); // minus sign because derivatice wrt xc  , VALIDATED
	 
	}
	PiPj[index]=Pi*Pj;
	dPiPjdsx[index]=dPidsx*Pj;
	dPiPjdsy[index]=Pi*dPjdsy;
	dPiPjdx[index]=dPidx*Pj;
	dPiPjdy[index]=Pi*dPjdy;
	
      }
    }
    
    psfval = ex*ey + specex::dot(PiPj,ublas::project(Params,ublas::range(first_hermite_param_index,first_hermite_param_index+nc)));
    
    if(ParamDer) {
      // derivative wrt gauss-hermite coefficients , VALIDATED
      ublas::project(*ParamDer,ublas::range(first_hermite_param_index,first_hermite_param_index+nc)) = PiPj;
      
      // derivative wrt sigmax, VALIDATED
      (*ParamDer)[0] = dexdsx*ey + specex::dot(dPiPjdsx,ublas::project(Params,ublas::range(first_hermite_param_index,first_hermite_param_index+nc)));
      // derivative wrt sigmay, VALIDATED
      (*ParamDer)[1] += ex*deydsy + specex::dot(dPiPjdsy,ublas::project(Params,ublas::range(first_hermite_param_index,first_hermite_param_index+nc)));
    }
    if(PosDer) {
      // derivative wrt x
      (*PosDer)[0] = dexdx*ey + specex::dot(dPiPjdx,ublas::project(Params,ublas::range(first_hermite_param_index,first_hermite_param_index+nc)));
      // derivative wrt y
      (*PosDer)[1] += ex*deydy + specex::dot(dPiPjdy,ublas::project(Params,ublas::range(first_hermite_param_index,first_hermite_param_index+nc)));
    }
    
    
    return psfval;

}
  
int specex::GaussHermitePSF::LocalNAllPar() const {
    
  int npar = 2; // sigma_x and sigma_y
  npar += (degree+1)*(degree+1)-1;// skip (0,0) : -1 because normalized
  
#ifdef EXTERNAL_TAIL
  npar += 5;
#endif

  return npar;
}

harp::vector_double specex::GaussHermitePSF::DefaultParams() const 
{
  
  harp::vector_double Params(LocalNAllPar());
  Params.clear(); // all = zero at beginning = a pure gaussian
  int index=0;
  Params(index++) = 1.0; // this is sigma_x ; value of 1. tuned on CCDS1R (EM-spectro)
  Params(index++) = 1.0; // this is sigma_y ; value of 1. tuned on CCDS1R (EM-spectro)
  index += ((degree+1)*(degree+1)-1);
  
#ifdef EXTERNAL_TAIL
  
  Params(index++) = 0.; // tail amplitude
  Params(index++) = 1.; // tail core size
  Params(index++) = 1.; // tail x scale
  Params(index++) = 1.; // tail y scale
  Params(index++) = 2.; // tail power law index
#endif

  return Params;
}

std::vector<std::string> specex::GaussHermitePSF::DefaultParamNames() const
{
  std::vector<std::string> paramNames;
  paramNames.push_back("GHSIGX");
  paramNames.push_back("GHSIGY");
      
  char n[10];
  for(int j=0;j<degree+1;j++) {
    for(int i=0;i<degree+1;i++) {
      if(i==0 && j==0) continue;
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

void specex::GaussHermitePSF::Append(const specex::PSF_p other_p) {
  SPECEX_INFO("GaussHermitePSF::Append starting");

  if(Name() != other_p->Name()) SPECEX_ERROR("Cannot append two different kind of PSF " << Name() << " " << other_p->Name());
  // casting
  const specex::GaussHermitePSF& other = (const specex::GaussHermitePSF&)(*other_p);
  if(degree != other.degree) SPECEX_ERROR("Cannot append two GaussHermitePSF of different Hermite degree " << degree << " " << other.degree );
  if(LocalNAllPar() != other.LocalNAllPar()) SPECEX_ERROR("Cannot append two PSF with different number of parameters " << LocalNAllPar()  << " " << other.LocalNAllPar() );
 
  if(arc_exposure_id != other.arc_exposure_id) SPECEX_ERROR("Cannot append two PSF with different arc exposure id " <<  arc_exposure_id << " " << other.arc_exposure_id);
  
  // append fibertraces
  for(std::map<int,Trace>::const_iterator it=other.FiberTraces.begin(); it!=other.FiberTraces.end(); it++) {
    
    if(FiberTraces.find(it->first) != FiberTraces.end()) {
      SPECEX_WARNING("Ignore same fiber trace info in GaussHermitePSF::Append for fiber " << it->first);
      continue;
    }
    SPECEX_INFO("Appending trace of fiber " << it->first);
    FiberTraces[it->first] = it->second;
  } 

  // append parameters
  for(std::map<int,PSF_Params>::const_iterator it = other.ParamsOfBundles.begin(); it != other.ParamsOfBundles.end(); ++it) {
    if(ParamsOfBundles.find(it->first) != ParamsOfBundles.end()) {
      SPECEX_WARNING("Ignore parameters of same bundle in GaussHermitePSF::Append for bundle " << it->first);
      continue;
    }
    SPECEX_INFO("Appending parameters of bundle " << it->first);
    ParamsOfBundles[it->first] = it->second;
  }

  
  SPECEX_INFO("GaussHermitePSF::Append successful");
  
}
