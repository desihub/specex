#include <iostream>
#include <fstream>
#include <ctime>

#include "harp.hpp"

#include "specex_message.h"
#include "specex_fits.h"
#include "specex_psf.h"
#include "specex_psf_io.h"
#include "specex_gauss_hermite_psf.h"
#include "specex_image_data.h"
#include "specex_serialization.h"

using namespace std;

int main() {
  
  
  specex_set_verbose(true);
  
  specex::PSF_p psf;
  specex::read_psf_xml(psf,"psf.xml");

  double x=400;
  double y=1000;
  int bundle=1;
  
  int i=int(x)+1;
  int j=int(y)+1;

  harp::vector_double P = psf->FixedCoordParams(x,y,bundle);
  cout << P << endl;

  harp::vector_double ParamDer(P.size());
  harp::vector_double PosDer(2);
  
  double val0 = psf->PSFValueWithParams(x,y,i,j,P,&PosDer,&ParamDer);
  
  for(int k=0;k<int(P.size());k++) {
    double eps = fabs(P(k))/1000.;
    harp::vector_double Pp = P; Pp(k)+=eps/2.;
    double valp = psf->PSFValueWithParams(x,y,i,j,Pp,0,0);
    harp::vector_double Pm = P; Pm(k)-=eps/2.;
    double valm = psf->PSFValueWithParams(x,y,i,j,Pm,0,0);
    
    double numDer = (valp-valm)/eps;

    cout << k << " P=" << P(k) << " aD=" << ParamDer(k) << " nD=" << numDer 
	 << " diff=" <<  ParamDer(k)-numDer << " ratio=" << (ParamDer(k)-numDer)/numDer << endl;
  }
   cout << "----------" << endl;

  
  double eps    = 0.001;
  {
    double valxp  = psf->PSFValueWithParams(x+0.5*eps,y,i,j,P,0,0);
    double valxm  = psf->PSFValueWithParams(x-0.5*eps,y,i,j,P,0,0);
    double numDer = (valxp-valxm)/eps;
    cout << "dPdx aD=" << PosDer(0) << " nD=" << numDer
	 << " diff=" <<  PosDer(0)-numDer << " ratio=" << PosDer(0)/numDer-1 << endl;
  }
  {
    double valyp  = psf->PSFValueWithParams(x,y+0.5*eps,i,j,P,0,0);
    double valym  = psf->PSFValueWithParams(x,y-0.5*eps,i,j,P,0,0);
    double numDer = (valyp-valym)/eps;
    cout << "dPdy aD=" << PosDer(1) << " nD=" << numDer
	 << " diff=" <<  PosDer(1)-numDer << " ratio=" << PosDer(1)/numDer-1 << endl;
  }
  
  
  return 0;
}
