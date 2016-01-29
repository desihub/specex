#include <iostream>
#include <fstream>
#include <ctime>

#include <harp.hpp>

#include "specex_message.h"
#include "specex_fits.h"
#include "specex_psf.h"
#include "specex_psf_io.h"
#include "specex_image_data.h"
#include "specex_serialization.h"

using namespace std;

int main() {
  
  
  specex_set_verbose(true);
  
  //specex::PSF_p psf(new specex::HatHermitePSF(3));
  //specex::PSF_p psf(new specex::GaussHermite2PSF(3,3));
  specex::PSF_p psf(new specex::GaussHermitePSF(3));
  //specex::PSF_p psf(new specex::HatMoffatPSF(3));
  psf->AllocateDefaultParams();
  
  vector<string> pnames = psf->DefaultParamNames();

  double x = 0;
  double y = 0;
  
  int i=int(x)+1;
  int j=int(y)+1;

  harp::vector_double P = psf->DefaultParams();

  // test for derivatives of sigma
  P(0)=1;
  P(1)=1.5;
  P(17)=0;
  
  for(size_t i=2;i<min(6,int(P.size()));i++) P(i)=4;
  //P(2)=3.;
  //P(20)=0.5;
  
  cout << P << endl;

  harp::vector_double ParamDer(P.size());
  harp::vector_double PosDer(2);
  
  psf->PSFValueWithParamsXY(x,y,i,j,P,&PosDer,&ParamDer,true,true);

  for(int k=0;k<int(P.size());k++) {
    double eps = fabs(P(k))/10000.;
    if(eps<1e-5) eps = 1.e-5;
    harp::vector_double Pp = P; Pp(k)+=eps/2.;
    double valp = psf->PSFValueWithParamsXY(x,y,i,j,Pp,0,0,true,true);
    harp::vector_double Pm = P; Pm(k)-=eps/2.;
    double valm = psf->PSFValueWithParamsXY(x,y,i,j,Pm,0,0,true,true);
    cout << eps << " " << valp << " " << valm << endl;
    double numDer = (valp-valm)/eps;

    cout << k << " " << pnames[k] << " P=" << P(k) << " aD=" << ParamDer(k) << " nD=" << numDer 
	 << " diff=" <<  ParamDer(k)-numDer << " ratio=" << (ParamDer(k)-numDer)/numDer << endl;
  }
  cout << "----------" << endl;

  //return 0;
  
  double eps    = 0.0001;
  {
    double valxp  = psf->PSFValueWithParamsXY(x+0.5*eps,y,i,j,P,0,0,true,true);
    double valxm  = psf->PSFValueWithParamsXY(x-0.5*eps,y,i,j,P,0,0,true,true);
    double numDer = (valxp-valxm)/eps;
    cout << "dPdx aD=" << PosDer(0) << " nD=" << numDer
	 << " diff=" <<  PosDer(0)-numDer << " ratio=" << PosDer(0)/numDer-1 << endl;
  }
  {
    double valyp  = psf->PSFValueWithParamsXY(x,y+0.5*eps,i,j,P,0,0,true,true);
    double valym  = psf->PSFValueWithParamsXY(x,y-0.5*eps,i,j,P,0,0,true,true);
    double numDer = (valyp-valym)/eps;
    cout << "dPdy aD=" << PosDer(1) << " nD=" << numDer
	 << " diff=" <<  PosDer(1)-numDer << " ratio=" << PosDer(1)/numDer-1 << endl;
  }
  return 0;
}
