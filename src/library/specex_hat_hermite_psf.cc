
#include <cmath>
#include <assert.h>
#include <ctime>



#include "harp.hpp"

#include "specex_hermite.h"
#include "specex_hat_hermite_psf.h"
#include "specex_linalg.h"
#include "specex_message.h"
#include "specex_fits.h"

using namespace std;



specex::HatHermitePSF::HatHermitePSF(int ideg) {
  name = "HatHermitePSF";
  SetDegree(ideg);
  
}

//#define HAT_AND_GAUSSIAN

void specex::HatHermitePSF::SetDegree(const int ideg) {
  degree = ideg;
}

int specex::HatHermitePSF::LocalNAllPar() const {
    
  int npar = 1; // sigma
  
#ifdef HAT_AND_GAUSSIAN
  npar += 3;
#endif

  // npar += (degree+1)*(degree+1)-4;// skip (0,0)(1,0)(0,1)(1,1) : -1 because normalized, -3 because centered 
  npar += (degree+1)*(degree+1)-1;// skip (0,0) : -1 because normalized
  

  return npar;
}

double specex::HatHermitePSF::Profile(const double &x, const double &y,
			  const harp::vector_double &Params,
			  harp::vector_double *PosDer,
			  harp::vector_double *ParamDer) const
{

  
  const double &sr = Params[0];
  double x2 = x*x;
  double y2 = y*y;
  double r2 = x2+y2;
  double r  = sqrt(r2);
  double sr2 = sr*sr;
  
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
  double fr = 1.53; // fiber radius, hardcoded for test purposes, this number of for b-channel
  double fr2 = fr*fr;
  
  double s2si     = 1./(sqrt(2)*sr); 
  double erfrs    = erf(fr/sr);
  double exprs    = exp(-fr2/sr2);


  /*
    maple
    assume(fr>0);
    assume(sr>0);
    
    hat_norm_inv:=2*sqrt(Pi)*fr*sr*exp(-fr**2/sr**2)+Pi*(2*fr**2+sr**2)*erf(fr/rs);



   */
  double hat_norm = 1/(2*sqrt(M_PI)*fr*sr*exprs+M_PI*(2*fr2+sr2)*erfrs);
  /*
    approximate radial profile of disk of radius fiber_radius and symmetric gaussian of sigma sr.
    the true profile is not analytic.
  */
  double hat = hat_norm*(erf(s2si*(fr-r))+erf(s2si*(fr+r)));
  double prof = hat;
  
  /* for hermite terms we use an effective sigma=1 anyway*/
  

#ifdef HAT_AND_GAUSSIAN
  double sigma_x_2_inv = 1./Params[1];
  double sigma_y_2_inv = 1./Params[2];
  double scale_2       = Params[3];
  double scale_1       = 1-scale_2;
  double x_2 = x*sigma_x_2_inv;
  double y_2 = y*sigma_y_2_inv;
  
  double gaus=1/(2*M_PI)*sigma_x_2_inv*sigma_y_2_inv*exp(-0.5*(x_2*x_2+y_2*y_2));
  
  prof = scale_1*hat+scale_2*gaus;
  
  
#endif
  
  int nx=(degree+1);
  int ny=(degree+1);
  //int nc=nx*ny-4; // skip (0,0)(1,0)(0,1)(1,1)
  int nc=nx*ny-1; // skip (0,0)
  
  // precompute to go faster
  harp::vector_double Monomials;
  harp::vector_double Monomials_dx;
  harp::vector_double Monomials_dy;
  
#ifndef HAT_AND_GAUSSIAN
  int first_hermite_param_index = 1; // first 1 param is sigma
#else
  int first_hermite_param_index = 4; // 
#endif
  
  double prefactor=1;
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
    return prof*prefactor;
    
  } else if(ParamDer) {
    Monomials.resize(nc);
    Monomials_dx.resize(nc); // need this for derivative of sigma
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
  }
  
  
  double psf_val = prof*prefactor;
  
  double exp1 = 0;
  double exp2 = 0;
  
  if(ParamDer || PosDer) {
    
    exp1 = exp(-square(s2si*(fr-r)));
    exp2 = exp(-square(s2si*(fr+r)));
  }
  
  if(ParamDer) {
    
    
    
    
    // der wrt erf terms:
    // MEMO psf_val = prefactor*hat_norm*(erf(s2si*(fr-r))+erf(s2si*(fr+r)));
    
    double dvdsr = (prefactor*hat_norm*s2si*2/sqrt(M_PI))*(-1./sr)*((fr-r)*exp1+(fr+r)*exp2); 
    
    // der wrt hat_norm:
    // MEMO hat_norm = 1/(2*sqrt(M_PI)*fr*sr*exprs+M_PI*(2*fr2+sr2)*erfrs);
    // MEMO s2si     = 1./(sqrt(2)*sr); 
    // MEMO exprs    = exp(-fr2/sr2);
    // MEMO erfrs    = erf(fr/sr);
    
    /*
    dvdsr += prefactor*hat*(-hat_norm)*( 2*sqrt(M_PI)*fr*exprs*(1+2*fr2/sr2)
					 + M_PI*2*sr*erfrs
					 + M_PI*(2*fr2+sr2)*(-fr/sr2)*2/sqrt(M_PI)*exp(-fr2/sr2)
					 );
    */
    
    // simplification :
    
    dvdsr += prefactor*hat*(-hat_norm)*( M_PI*2*sr*erfrs );
        

#ifdef HAT_AND_GAUSSIAN
    dvdsr *= scale_1;
#endif 


    (*ParamDer)[0] = dvdsr;
    
#ifdef HAT_AND_GAUSSIAN
    // need to compute derivatives here if we care
    (*ParamDer)[1] = 0;
    (*ParamDer)[2] = 0;
    (*ParamDer)[3] = (-hat+gaus)*prefactor;
#endif 
    
    ublas::project(*ParamDer,ublas::range(first_hermite_param_index,first_hermite_param_index+nc)) = prof*Monomials;
  }
  
  if(PosDer) { 
    double& dvdx=(*PosDer)(0);
    double& dvdy=(*PosDer)(1);  
   

   
    // MEMO psf_val = prefactor*hat_norm*erf(s2si*(fr-r))+erf(s2si*(fr+r));
    double dvdr  = (prefactor*hat_norm*s2si*2/sqrt(M_PI))*(-exp1+exp2);

#ifdef HAT_AND_GAUSSIAN 
    dvdr *= scale_1;
#endif
    
    dvdx=-(x/r)*dvdr; // minus sign cause derivative wrt -x
    dvdy=-(y/r)*dvdr; // minus sign cause derivative wrt -y
    
    double d_poly_dx = specex::dot(Monomials_dx,ublas::project(Params,ublas::range(first_hermite_param_index,first_hermite_param_index+nc)));
    double d_poly_dy = specex::dot(Monomials_dy,ublas::project(Params,ublas::range(first_hermite_param_index,first_hermite_param_index+nc)));
    dvdx -= d_poly_dx*prof; // minus sign cause derivative wrt -x
    dvdy -= d_poly_dy*prof; // minus sign cause derivative wrt -y

#ifdef HAT_AND_GAUSSIAN
    // derivative due to gaussian terms
    
    dvdx += x_2*sigma_x_2_inv*scale_2*gaus*prefactor;
    dvdy += y_2*sigma_y_2_inv*scale_2*gaus*prefactor;

#endif

  }
  
  return psf_val;
}

  
harp::vector_double specex::HatHermitePSF::DefaultParams() const 
{
  
  harp::vector_double Params(LocalNAllPar());
  Params.clear(); // all = zero at beginning 
  Params(0) = 0.3; // this is sigma_x of the hat
#ifdef HAT_AND_GAUSSIAN
  Params(1) = 3.; // this is sigma_x_2
  Params(2) = 3.; // this is sigma_y_2
  Params(3) = 0.1; // this is the amplitude of the second gaussian
#endif
  return Params;
}

std::vector<std::string> specex::HatHermitePSF::DefaultParamNames() const
{
  std::vector<std::string> paramNames;
  paramNames.push_back("GHSIGX");
  
#ifdef HAT_AND_GAUSSIAN
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

void specex::HatHermitePSF::Append(const specex::PSF_p other) {
  SPECEX_ERROR("specex::HatHermitePSF::Append not implemented");
}
