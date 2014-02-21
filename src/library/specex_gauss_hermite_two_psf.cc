
#include <cmath>
#include <assert.h>
#include <ctime>



#include "harp.hpp"

#include "specex_hermite.h"
#include "specex_gauss_hermite_two_psf.h"
#include "specex_linalg.h"
#include "specex_message.h"
#include "specex_fits.h"
#include "specex_psf.h"

using namespace std;



specex::GaussHermite2PSF::GaussHermite2PSF(int i_core_degree, int i_second_degree)  : PSF() {
  name = "GaussHermite2PSF";
  SetDegree(i_core_degree,i_second_degree);
  
}

void specex::GaussHermite2PSF::SetDegree(const int i_core_degree, const int i_second_degree) {
  core_degree = i_core_degree;
  second_degree = i_second_degree;
}



double specex::GaussHermite2PSF::Profile(const double &input_X, const double &input_Y,
			  const harp::vector_double &Params,
			  harp::vector_double *PosDer,
			  harp::vector_double *ParamDer) const
{

  double sigma_x_1_inv = 1./Params(0);
  double sigma_y_1_inv = 1./Params(1);
  double x_1 = input_X*sigma_x_1_inv;
  double y_1 = input_Y*sigma_y_1_inv;

  double sigma_x_2_inv = 1./Params(2);
  double sigma_y_2_inv = 1./Params(3);
  double x_2 = input_X*sigma_x_2_inv;
  double y_2 = input_Y*sigma_y_2_inv;
  
  int first_hermite1_param_index = 4;
  int first_hermite2_param_index = 4+((core_degree+1)*(core_degree+1)-1);
  
  double scale2       = 1;
  double scale1       = 1-Params(first_hermite2_param_index);
  
  int nx1=(core_degree+1);
  int ny1=(core_degree+1);
  int nc1=nx1*ny1-1; // skip (0,0)
  
  int nx2=(second_degree+1);
  int ny2=(second_degree+1);
  int nc2=nx2*ny2; // don't skip (0,0)
  
  // precompute to go faster
  harp::vector_double Monomials1;
  harp::vector_double Monomials1_dx;
  harp::vector_double Monomials1_dy;
  harp::vector_double Monomials2;
  harp::vector_double Monomials2_dx;
  harp::vector_double Monomials2_dy;
  
  double prefactor1=1; // includes order 0
  double prefactor2=0;
  double expfact1=1./(2*M_PI)*sigma_x_1_inv*sigma_y_1_inv*exp(-0.5*(x_1*x_1+y_1*y_1));
  double expfact2=1./(2*M_PI)*sigma_x_2_inv*sigma_y_2_inv*exp(-0.5*(x_2*x_2+y_2*y_2));
  
  
  if(PosDer==0 && ParamDer==0) {
    int param_index;
    param_index=first_hermite1_param_index;
    for(int j=0;j<ny1;j++) {
      double Hyj=HermitePol(j,y_1);
      int imin=0; if(j==0) imin=1; // skip (0,0)
      for(int i=imin;i<nx1;i++,param_index++) {
	prefactor1+=Params(param_index)*Hyj*HermitePol(i,x_1);
      }
    }
    param_index=first_hermite2_param_index;
    for(int j=0;j<ny2;j++) {
      double Hyj=HermitePol(j,y_2);
      for(int i=0;i<nx2;i++,param_index++) {
	prefactor2+=Params(param_index)*Hyj*HermitePol(i,x_2);
      }
    }
    
    return scale1*expfact1*prefactor1+scale2*expfact2*prefactor2;
    
  } else if(ParamDer) {
    Monomials1.resize(nc1);
    Monomials1_dx.resize(nc1); // need this for derivative of sigma
    Monomials1_dy.resize(nc1);
    Monomials2.resize(nc2);
    Monomials2_dx.resize(nc2); // need this for derivative of sigma
    Monomials2_dy.resize(nc2);
    
    int index;
    
    index=0;
    for(int j=0;j<ny1;j++) {
      double Hyj=HermitePol(j,y_1);
      double dHyj=HermitePolDerivative(j,y_1);
      
      int imin=0; if(j==0) imin=1; // skip (0,0)
      for(int i=imin;i<nx1;i++,index++) {
	double Hxi=HermitePol(i,x_1);
	Monomials1[index]=Hyj*Hxi;
	Monomials1_dx[index]=Hyj*HermitePolDerivative(i,x_1);
	Monomials1_dy[index]=dHyj*Hxi;
      }
    }
    
    index=0;
    for(int j=0;j<ny2;j++) {
      double Hyj=HermitePol(j,y_2);
      double dHyj=HermitePolDerivative(j,y_2);
      
      for(int i=0;i<nx2;i++,index++) {
	double Hxi=HermitePol(i,x_2);
	Monomials2[index]=Hyj*Hxi;
	Monomials2_dx[index]=Hyj*HermitePolDerivative(i,x_2);
	Monomials2_dy[index]=dHyj*Hxi;
      }
    }
    prefactor1 += specex::dot(Monomials1,ublas::project(Params,ublas::range(first_hermite1_param_index,first_hermite1_param_index+nc1)));
    prefactor2 += specex::dot(Monomials2,ublas::project(Params,ublas::range(first_hermite2_param_index,first_hermite2_param_index+nc2)));
    
  } else if(PosDer && ParamDer==0) {
    
    Monomials1_dx.resize(nc1);
    Monomials1_dy.resize(nc1);
    Monomials2_dx.resize(nc2);
    Monomials2_dy.resize(nc2);
    

    int param_index,index;
    
    param_index=first_hermite1_param_index;
    index=0;
    for(int j=0;j<ny1;j++) {
      double Hyj=HermitePol(j,y_1);
      double dHyj=HermitePolDerivative(j,y_1);
      
      int imin=0; if(j==0) imin=1; // skip (0,0)
      for(int i=imin;i<nx1;i++,index++,param_index++) {
	
	double Hxi=HermitePol(i,x_1);
	prefactor1+=Params(param_index)*Hxi*Hyj;
	Monomials1_dx[index]=Hyj*HermitePolDerivative(i,x_1);
	Monomials1_dy[index]=dHyj*Hxi;
      }
    }

    param_index=first_hermite2_param_index;
    index=0;
    for(int j=0;j<ny2;j++) {
      double Hyj=HermitePol(j,y_2);
      double dHyj=HermitePolDerivative(j,y_2);
      
      for(int i=0;i<nx2;i++,index++,param_index++) {
	
	double Hxi=HermitePol(i,x_2);
	prefactor2+=Params(param_index)*Hxi*Hyj;
	Monomials2_dx[index]=Hyj*HermitePolDerivative(i,x_2);
	Monomials2_dy[index]=dHyj*Hxi;
      }
    }
  }

  double psf_val1 = scale1*expfact1*prefactor1;
  double psf_val2 = scale2*expfact2*prefactor2;

  double psf_val = psf_val1+psf_val2;
  
  if(ParamDer) {
    
    
    (*ParamDer)[0] = (x_1*x_1-1)*sigma_x_1_inv*psf_val1;
    (*ParamDer)[1] = (y_1*y_1-1)*sigma_y_1_inv*psf_val1;
    (*ParamDer)[2] = (x_2*x_2-1)*sigma_x_2_inv*psf_val2;
    (*ParamDer)[3] = (y_2*y_2-1)*sigma_y_2_inv*psf_val2;
    
    /* other terms 

       dPSD/sigmax += sum_ij Pij dH(x,i)/dsigmax H(y,j) * expfact
                   += dx/dsigmax * sum_ij Pij dH(x,i)dx *H(y,j) * expfact
                   -= x/sigmax * sum_i P_i Monomials_dx_i * expfact
     */
     
    (*ParamDer)[0] -= (x_1*sigma_x_1_inv*scale1*expfact1)*specex::dot(Monomials1_dx,ublas::project(Params,ublas::range(first_hermite1_param_index,first_hermite1_param_index+nc1)));
    (*ParamDer)[1] -= (y_1*sigma_y_1_inv*scale1*expfact1)*specex::dot(Monomials1_dy,ublas::project(Params,ublas::range(first_hermite1_param_index,first_hermite1_param_index+nc1)));
    (*ParamDer)[2] -= (x_2*sigma_x_2_inv*scale2*expfact2)*specex::dot(Monomials2_dx,ublas::project(Params,ublas::range(first_hermite2_param_index,first_hermite2_param_index+nc2)));
    (*ParamDer)[3] -= (y_2*sigma_y_2_inv*scale2*expfact2)*specex::dot(Monomials2_dy,ublas::project(Params,ublas::range(first_hermite2_param_index,first_hermite2_param_index+nc2)));
    
    ublas::project(*ParamDer,ublas::range(first_hermite1_param_index,first_hermite1_param_index+nc1)) = (scale1*expfact1)*Monomials1;
    ublas::project(*ParamDer,ublas::range(first_hermite2_param_index,first_hermite2_param_index+nc2)) = (scale2*expfact2)*Monomials2;

    // add one thing, change of scale1
    (*ParamDer)[first_hermite2_param_index] -= expfact1*prefactor1;
    
  }
  
  if(PosDer) { // wrong sign on purpose (derivatives w.r.t -X)
    double& dvdx=(*PosDer)(0);
    double& dvdy=(*PosDer)(1);  
   
    double d_poly_dx,d_poly_dy;

    dvdx=x_1*sigma_x_1_inv*psf_val1;
    dvdy=y_1*sigma_y_1_inv*psf_val1;
    
    d_poly_dx = specex::dot(Monomials1_dx,ublas::project(Params,ublas::range(first_hermite1_param_index,first_hermite1_param_index+nc1)));
    d_poly_dy = specex::dot(Monomials1_dy,ublas::project(Params,ublas::range(first_hermite1_param_index,first_hermite1_param_index+nc1)));
    dvdx -= d_poly_dx*scale1*expfact1*sigma_x_1_inv; // minus sign cause derivative wrt -x
    dvdy -= d_poly_dy*scale1*expfact1*sigma_y_1_inv;
    
    dvdx += x_2*sigma_x_2_inv*psf_val2;
    dvdy += y_2*sigma_y_2_inv*psf_val2;
    
    d_poly_dx = specex::dot(Monomials2_dx,ublas::project(Params,ublas::range(first_hermite2_param_index,first_hermite2_param_index+nc2)));
    d_poly_dy = specex::dot(Monomials2_dy,ublas::project(Params,ublas::range(first_hermite2_param_index,first_hermite2_param_index+nc2)));
    dvdx -= d_poly_dx*scale2*expfact2*sigma_x_2_inv; // minus sign cause derivative wrt -x
    dvdy -= d_poly_dy*scale2*expfact2*sigma_y_2_inv;
}
  
  return psf_val;
}

  
int specex::GaussHermite2PSF::LocalNAllPar() const {
    
  int npar = 4; // sigma_x, sigma_y, sigma_x_2, sigma_y_2
  
  npar += (core_degree+1)*(core_degree+1)-1;// skip (0,0) : -1 because normalized
  npar += (second_degree+1)*(second_degree+1);// don't skip (0,0)
  
#ifdef EXTERNAL_TAIL
  npar += 5;
#endif

  return npar;
}

harp::vector_double specex::GaussHermite2PSF::DefaultParams() const 
{
  
  harp::vector_double Params(LocalNAllPar());
  Params.clear(); // all = zero at beginning = a pure gaussian
  int index=0;
  Params(index++) = 1.1; // this is sigma_x
  Params(index++) = 1.1; // this is sigma_y  
  Params(index++) = 4.; // this is sigma_x_2
  Params(index++) = 4.; // this is sigma_y_2
  
  index += ((core_degree+1)*(core_degree+1)-1);
  index += ((second_degree+1)*(second_degree+1)); // here norm is a free param
  
#ifdef EXTERNAL_TAIL
  
  Params(index++) = 0.; // tail amplitude
  Params(index++) = 1.; // tail core size
  Params(index++) = 1.; // tail x scale
  Params(index++) = 1.; // tail y scale
  Params(index++) = 2.; // tail power law index
#endif

  return Params;
}

std::vector<std::string> specex::GaussHermite2PSF::DefaultParamNames() const
{
  std::vector<std::string> paramNames;
  paramNames.push_back("GHSIGX");
  paramNames.push_back("GHSIGY");
  paramNames.push_back("GHSIGX2");
  paramNames.push_back("GHSIGY2");

    
  char n[10];
  for(int j=0;j<core_degree+1;j++) {
    for(int i=0;i<core_degree+1;i++) {
      if(i==0 && j==0) continue;
      sprintf(n,"GH-%d-%d",i,j);
      paramNames.push_back(n);
    }
  }
  
  for(int j=0;j<second_degree+1;j++) {
    for(int i=0;i<second_degree+1;i++) {
      sprintf(n,"GH2-%d-%d",i,j);
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

void specex::GaussHermite2PSF::Append(const specex::PSF_p other_p) {
  SPECEX_INFO("GaussHermite2PSF::Append starting");

  if(Name() != other_p->Name()) SPECEX_ERROR("Cannot append two different kind of PSF " << Name() << " " << other_p->Name());
  // casting
  const specex::GaussHermite2PSF& other = (const specex::GaussHermite2PSF&)(*other_p);
  if(core_degree != other.core_degree) SPECEX_ERROR("Cannot append two GaussHermite2PSF of different Hermite core_degree " << core_degree << " " << other.core_degree );
  if(second_degree != other.second_degree) SPECEX_ERROR("Cannot append two GaussHermite2PSF of different Hermite second_degree " << second_degree << " " << other.second_degree );
  if(LocalNAllPar() != other.LocalNAllPar()) SPECEX_ERROR("Cannot append two PSF with different number of parameters " << LocalNAllPar()  << " " << other.LocalNAllPar() );
 
  if(arc_exposure_id != other.arc_exposure_id) SPECEX_ERROR("Cannot append two PSF with different arc expsure id " <<  arc_exposure_id << " " << other.arc_exposure_id);
  
  // append fibertraces
  for(std::map<int,Trace>::const_iterator it=other.FiberTraces.begin(); it!=other.FiberTraces.end(); it++) {
    
    if(FiberTraces.find(it->first) != FiberTraces.end()) {
      SPECEX_WARNING("Ignore same fiber trace info in GaussHermite2PSF::Append for fiber " << it->first);
      continue;
    }
    SPECEX_INFO("Appending trace of fiber " << it->first);
    FiberTraces[it->first] = it->second;
  } 

  // append parameters
  for(std::map<int,PSF_Params>::const_iterator it = other.ParamsOfBundles.begin(); it != other.ParamsOfBundles.end(); ++it) {
    if(ParamsOfBundles.find(it->first) != ParamsOfBundles.end()) {
      SPECEX_WARNING("Ignore parameters of same bundle in GaussHermite2PSF::Append for bundle " << it->first);
      continue;
    }
    SPECEX_INFO("Appending parameters of bundle " << it->first);
    ParamsOfBundles[it->first] = it->second;
  }

  
  SPECEX_INFO("GaussHermite2PSF::Append successful");
  
}
