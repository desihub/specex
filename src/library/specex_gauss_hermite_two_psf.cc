
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

// if we want to refit the sigmas and trace when hermite pols are not zero 
#ifdef FULL_DERIVATIVES
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
  
  double prefactor1=1-Params(first_hermite2_param_index); // includes order 0
  double prefactor2=0;
  double expfact1=1./(2*M_PI)*sigma_x_1_inv*sigma_y_1_inv*exp(-0.5*(x_1*x_1+y_1*y_1));
  double expfact2=1./(2*M_PI)*sigma_x_2_inv*sigma_y_2_inv*exp(-0.5*(x_2*x_2+y_2*y_2));
  
  if(ParamDer!=0 || PosDer!=0) {
    
    Monomials1.resize(nc1);
    Monomials1_dx.resize(nc1); // need this for derivative of sigma
    Monomials1_dy.resize(nc1);
    harp::vector_double H1x,H1y,dH1dx,dH1dy;
    HermitePolsAndDerivatives(H1x,dH1dx,core_degree,x_1);
    HermitePolsAndDerivatives(H1y,dH1dy,core_degree,y_1);
    int index=0;
    for(int j=0;j<ny1;j++) {
      int imin=0; if(j==0) imin=1; // skip (0,0)
      for(int i=imin;i<nx1;i++,index++) {
	Monomials1[index]=H1y[j]*H1x[i];
	Monomials1_dx[index]=H1y[j]*dH1dx[i];
	Monomials1_dy[index]=dH1dy[j]*H1x[i];
      }
    }
    prefactor1=1-Params(first_hermite2_param_index)+specex::dot(ublas::project(Params,ublas::range(first_hermite1_param_index,first_hermite1_param_index+nc1)),Monomials1);
    
    Monomials2.resize(nc2);
    Monomials2_dx.resize(nc2); // need this for derivative of sigma
    Monomials2_dy.resize(nc2);
    harp::vector_double H2x,H2y,dH2dx,dH2dy;
    HermitePolsAndDerivatives(H2x,dH2dx,second_degree,x_2);
    HermitePolsAndDerivatives(H2y,dH2dy,second_degree,y_2);
    index=0;
    for(int j=0;j<ny2;j++) {
      for(int i=0;i<nx2;i++,index++) {
	Monomials2[index]=H2y[j]*H2x[i];
	Monomials2_dx[index]=H2y[j]*dH2dx[i];
	Monomials2_dy[index]=dH2dy[j]*H2x[i];
      }
    }
    prefactor2=specex::dot(ublas::project(Params,ublas::range(first_hermite2_param_index,first_hermite2_param_index+nc2)),Monomials2);
    
  }else{
    //Monomials1.resize(nc1);
    harp::vector_double H1x,H1y;
    HermitePols(H1x,core_degree,x_1);
    HermitePols(H1y,core_degree,y_1);
    int index=first_hermite1_param_index;
    prefactor1=1-Params(first_hermite2_param_index);
    for(int j=0;j<ny1;j++) {
      int imin=0; if(j==0) imin=1; // skip (0,0)
      for(int i=imin;i<nx1;i++,index++) {
	prefactor1+=Params(index)*H1y[j]*H1x[i];
      }
    }
    
    Monomials2.resize(nc2);
    harp::vector_double H2x,H2y;
    HermitePols(H2x,second_degree,x_2);
    HermitePols(H2y,second_degree,y_2);
    index=first_hermite2_param_index;
    prefactor2=0;
    for(int j=0;j<ny2;j++) {
      for(int i=0;i<nx2;i++,index++) {
	prefactor2+=Params(index)*H2y[j]*H2x[i];
      }
    }
    return expfact1*prefactor1+expfact2*prefactor2;
  }
  
  double psf_val1 = expfact1*prefactor1;
  double psf_val2 = expfact2*prefactor2;

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
     
    (*ParamDer)[0] -= (x_1*sigma_x_1_inv*expfact1)*specex::dot(Monomials1_dx,ublas::project(Params,ublas::range(first_hermite1_param_index,first_hermite1_param_index+nc1)));
    (*ParamDer)[1] -= (y_1*sigma_y_1_inv*expfact1)*specex::dot(Monomials1_dy,ublas::project(Params,ublas::range(first_hermite1_param_index,first_hermite1_param_index+nc1)));
    (*ParamDer)[2] -= (x_2*sigma_x_2_inv*expfact2)*specex::dot(Monomials2_dx,ublas::project(Params,ublas::range(first_hermite2_param_index,first_hermite2_param_index+nc2)));
    (*ParamDer)[3] -= (y_2*sigma_y_2_inv*expfact2)*specex::dot(Monomials2_dy,ublas::project(Params,ublas::range(first_hermite2_param_index,first_hermite2_param_index+nc2)));
    
    ublas::project(*ParamDer,ublas::range(first_hermite1_param_index,first_hermite1_param_index+nc1)) = expfact1*Monomials1;
    ublas::project(*ParamDer,ublas::range(first_hermite2_param_index,first_hermite2_param_index+nc2)) = expfact2*Monomials2;

    // add one thing, change of scale1
    (*ParamDer)[first_hermite2_param_index] -= expfact1;
    
  }
  
  if(PosDer) { // wrong sign on purpose (derivatives w.r.t -X)
    double& dvdx=(*PosDer)(0);
    double& dvdy=(*PosDer)(1);  
   
    double d_poly_dx,d_poly_dy;

    dvdx=x_1*sigma_x_1_inv*psf_val1;
    dvdy=y_1*sigma_y_1_inv*psf_val1;
    
    d_poly_dx = specex::dot(Monomials1_dx,ublas::project(Params,ublas::range(first_hermite1_param_index,first_hermite1_param_index+nc1)));
    d_poly_dy = specex::dot(Monomials1_dy,ublas::project(Params,ublas::range(first_hermite1_param_index,first_hermite1_param_index+nc1)));
    dvdx -= d_poly_dx*expfact1*sigma_x_1_inv; // minus sign cause derivative wrt -x
    dvdy -= d_poly_dy*expfact1*sigma_y_1_inv;
    
    dvdx += x_2*sigma_x_2_inv*psf_val2;
    dvdy += y_2*sigma_y_2_inv*psf_val2;
    
    d_poly_dx = specex::dot(Monomials2_dx,ublas::project(Params,ublas::range(first_hermite2_param_index,first_hermite2_param_index+nc2)));
    d_poly_dy = specex::dot(Monomials2_dy,ublas::project(Params,ublas::range(first_hermite2_param_index,first_hermite2_param_index+nc2)));
    dvdx -= d_poly_dx*expfact2*sigma_x_2_inv; // minus sign cause derivative wrt -x
    dvdy -= d_poly_dy*expfact2*sigma_y_2_inv;
  }
  
  return psf_val;
}
#else
// faster version where we don't compute everything
double specex::GaussHermite2PSF::Profile(const double &input_X, const double &input_Y,
			  const harp::vector_double &Params,
			  harp::vector_double *PosDer,
			  harp::vector_double *ParamDer) const
{

  const double& inner_core_n_sig = Params(2);
  double sigma_x_1_inv = 1./Params(0);
  double sigma_y_1_inv = 1./Params(1);
  double x_1 = input_X*sigma_x_1_inv;
  double y_1 = input_Y*sigma_y_1_inv;
  bool in_inner_core = (x_1*x_1+y_1*y_1<inner_core_n_sig*inner_core_n_sig);
    
  int first_hermite1_param_index = 5;
  int first_hermite2_param_index = 5+((core_degree+1)*(core_degree+1)-1);
  
  double Hx0,Hx1,Hx;
  double Hy0,Hy1,Hy;
  
  double psf_val1=0;
  
  if(ParamDer) ParamDer->clear();
  if(PosDer) PosDer->clear();

  if(in_inner_core) {
    
    int nx1=(core_degree+1);
    int ny1=(core_degree+1);
      
    double expfact1=1./(2*M_PI)*sigma_x_1_inv*sigma_y_1_inv*exp(-0.5*(x_1*x_1+y_1*y_1));
    int index=first_hermite1_param_index;
    psf_val1=1-Params(first_hermite2_param_index);

    Hy0=1;
    Hy1=y_1;
    for(int j=0;j<ny1;j++) {
      if(j==0) { 
	Hy=Hy0;
      } else if(j==1) {
	Hy=Hy1;
      } else {
	Hy=y_1*Hy1-(j-1)*Hy0;
	Hy0=Hy1;
	Hy1=Hy;
      }
      
      Hx0=1;
      Hx1=x_1;
      
      int imin=0; if(j==0) imin=1; // skip (0,0)
      for(int i=imin;i<nx1;i++,index++) {
	
	if(i==0) { 
	  Hx=Hx0;
	} else if(i==1) {
	  Hx=Hx1;
	} else {
	  Hx=x_1*Hx1-(i-1)*Hx0;
	  Hx0=Hx1;
	  Hx1=Hx;
      }
	
	psf_val1+=Params[index]*Hy*Hx;
	if(ParamDer) (*ParamDer)[index]=expfact1*Hy*Hx;
      }
    }
  
    psf_val1 *= expfact1;
    if(ParamDer) { 
      // ok only if all hermite terms with deg>0 are = 0
      (*ParamDer)[0] = (x_1*x_1-1)*sigma_x_1_inv*psf_val1;
      (*ParamDer)[1] = (y_1*y_1-1)*sigma_y_1_inv*psf_val1;
      // add one thing, change of scale1
      (*ParamDer)[first_hermite2_param_index] -= expfact1;
    }
    if(PosDer) { // wrong sign on purpose (derivatives w.r.t -X)
      // works ONLY if for deg>0 hermite coeff=0
      (*PosDer)(0) = x_1*sigma_x_1_inv*psf_val1;
      (*PosDer)(1) = y_1*sigma_y_1_inv*psf_val1;
    }
  } // end of test of inner core
    
  
  // second gauss-hermite
  
  double sigma_x_2_inv = 1./Params(3);
  double sigma_y_2_inv = 1./Params(4);
  double x_2 = input_X*sigma_x_2_inv;
  double y_2 = input_Y*sigma_y_2_inv;
  
  double expfact2=1./(2*M_PI)*sigma_x_2_inv*sigma_y_2_inv*exp(-0.5*(x_2*x_2+y_2*y_2));
  
  int nx2=(second_degree+1);
  int ny2=(second_degree+1);

  int index=first_hermite2_param_index;
  double psf_val2=0;  
  Hy0=1;
  Hy1=y_2;
  
  for(int j=0;j<ny2;j++) {
    if(j==0) { 
      Hy=Hy0;
    } else if(j==1) {
      Hy=Hy1;
    } else {
      Hy=y_2*Hy1-(j-1)*Hy0;
      Hy0=Hy1;
      Hy1=Hy;
    }
    
    Hx0=1;
    Hx1=x_2;
    
    for(int i=0;i<nx2;i++,index++) {
      
      if(i==0) { 
	Hx=Hx0;
      } else if(i==1) {
	Hx=Hx1;
      } else {
	Hx=x_2*Hx1-(i-1)*Hx0;
	Hx0=Hx1;
	Hx1=Hx;
      }
      
      psf_val2+=Params[index]*Hy*Hx;
      if(ParamDer) (*ParamDer)[index]+=expfact2*Hy*Hx;
    }
  }
  
  
  psf_val2 *= expfact2;
    
  if(ParamDer) { 
    // ok only if all hermite terms with deg>0 are = 0
    (*ParamDer)[3] = (x_2*x_2-1)*sigma_x_2_inv*psf_val2;
    (*ParamDer)[4] = (y_2*y_2-1)*sigma_y_2_inv*psf_val2;
  }
  if(PosDer) { // wrong sign on purpose (derivatives w.r.t -X)
    // works ONLY if for deg>0 hermite coeff=0
    (*PosDer)(0) += x_2*sigma_x_2_inv*psf_val2;
    (*PosDer)(1) += y_2*sigma_y_2_inv*psf_val2;  
  }
  
  return psf_val1+psf_val2;
}
#endif

  
int specex::GaussHermite2PSF::LocalNAllPar() const {
    
  int npar = 5; // sigma_x, sigma_y, inner_core, sigma_x_2, sigma_y_2
  
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
  Params(index++) = 1000; // inner core radius
  Params(index++) = 3.; // this is sigma_x_2
  Params(index++) = 3.; // this is sigma_y_2
  
  index += ((core_degree+1)*(core_degree+1)-1);
  index += ((second_degree+1)*(second_degree+1)); // here norm is a free param
  
#ifdef EXTERNAL_TAIL
  
  Params(index++) = 0.01; // tail amplitude
  Params(index++) = 1.; // tail core size
  Params(index++) = 1.; // tail x scale
  Params(index++) = 1.; // tail y scale
  Params(index++) = 2.2; // tail power law index
#endif

  return Params;
}

std::vector<std::string> specex::GaussHermite2PSF::DefaultParamNames() const
{
  std::vector<std::string> paramNames;
  paramNames.push_back("GHSIGX");
  paramNames.push_back("GHSIGY");
  paramNames.push_back("GHNSIG");
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

//! Access to the current PSF, with user provided Params.
double specex::GaussHermite2PSF::PSFValueWithParamsXY(const double &Xc, const double &Yc, 
					 const int IPix, const int JPix,
					 const harp::vector_double &Params,
					 harp::vector_double *PosDer, harp::vector_double *ParamDer,
					 bool with_core, bool with_tail) const {
  
  if(PosDer) PosDer->clear();
  if(ParamDer) ParamDer->clear();
  
  //here is a hack
#define NASTY_HACK_FOR_SPECTER
#ifdef NASTY_HACK_FOR_SPECTER
#warning nasty hack to get rid off asap (i am sure it will last forever)
  double  orig_nsig=Params(2);
  double& mod_nsig=const_cast<harp::vector_double &>(Params)(2);
  double xPixCenter = floor(IPix+0.5);
  double yPixCenter = floor(JPix+0.5);
  if(square(Xc-xPixCenter)+square(Yc-yPixCenter)<square(orig_nsig))
    mod_nsig=1000; // in inner core everywhere in the pixel
  else
    mod_nsig=0; // out of inner core everywhere in the pixel
#endif
  
  double val = 0;
  if(with_core) val += PixValue(Xc,Yc,IPix, JPix, Params, PosDer, ParamDer); 

#ifdef NASTY_HACK_FOR_SPECTER 
  mod_nsig=orig_nsig;
#endif     

#ifdef EXTERNAL_TAIL
#ifdef INTEGRATING_TAIL
  if(with_tail && !with_core) {
    double prof = TailProfile(IPix-Xc,JPix-Yc, Params, with_core);
    if(ParamDer) (*ParamDer)(psf_tail_amplitude_index) = prof;
    val += Params(psf_tail_amplitude_index)*prof;
  }
#else
  if(with_tail) {
    double prof = TailProfile(IPix-Xc,JPix-Yc, Params, with_core);
    if(ParamDer) (*ParamDer)(psf_tail_amplitude_index) = prof;
    val += Params(psf_tail_amplitude_index)*prof;
  }
#endif
#endif
  return val;
}
