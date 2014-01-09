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

  int fiber=20;
  double wave=6000;
  
  double x = psf->Xccd(fiber,wave);
  double y = psf->Yccd(fiber,wave);
  

  int bundle=1;
  
  int i=int(x)+1;
  int j=int(y)+1;

  harp::vector_double P = psf->FixedCoordParamsFW(fiber,wave,bundle);
  cout << P << endl;

  harp::vector_double ParamDer(P.size());
  harp::vector_double PosDer(2);
  
  double val0 = psf->PSFValueWithParamsXY(x,y,i,j,P,&PosDer,&ParamDer);
  
  for(int k=0;k<int(P.size());k++) {
    double eps = fabs(P(k))/1000.;
    harp::vector_double Pp = P; Pp(k)+=eps/2.;
    double valp = psf->PSFValueWithParamsXY(x,y,i,j,Pp,0,0);
    harp::vector_double Pm = P; Pm(k)-=eps/2.;
    double valm = psf->PSFValueWithParamsXY(x,y,i,j,Pm,0,0);
    
    double numDer = (valp-valm)/eps;

    cout << k << " P=" << P(k) << " aD=" << ParamDer(k) << " nD=" << numDer 
	 << " diff=" <<  ParamDer(k)-numDer << " ratio=" << (ParamDer(k)-numDer)/numDer << endl;
  }
   cout << "----------" << endl;

  
  double eps    = 0.001;
  {
    double valxp  = psf->PSFValueWithParamsXY(x+0.5*eps,y,i,j,P,0,0);
    double valxm  = psf->PSFValueWithParamsXY(x-0.5*eps,y,i,j,P,0,0);
    double numDer = (valxp-valxm)/eps;
    cout << "dPdx aD=" << PosDer(0) << " nD=" << numDer
	 << " diff=" <<  PosDer(0)-numDer << " ratio=" << PosDer(0)/numDer-1 << endl;
  }
  {
    double valyp  = psf->PSFValueWithParamsXY(x,y+0.5*eps,i,j,P,0,0);
    double valym  = psf->PSFValueWithParamsXY(x,y-0.5*eps,i,j,P,0,0);
    double numDer = (valyp-valym)/eps;
    cout << "dPdy aD=" << PosDer(1) << " nD=" << numDer
	 << " diff=" <<  PosDer(1)-numDer << " ratio=" << PosDer(1)/numDer-1 << endl;
  }

  {
  cout  << "testing global monomials" << endl;
  specex::PSF_Params &psf_global_params = psf->ParamsOfBundles[bundle];
  harp::vector_double gP(psf->VaryingCoordNPar(bundle));
  harp::vector_double gM(psf->VaryingCoordNPar(bundle));
  int index=0;
  for(int p=0;p<psf->NPar();p++) {
    harp::vector_double pM = psf_global_params.Polynomials[p].Monomials(x,wave);
    int p_size = pM.size();
    ublas::project(gM,ublas::range(index,index+p_size))=pM;
    ublas::project(gP,ublas::range(index,index+p_size))=psf_global_params.Polynomials[p].coeff;
    index += p_size;
  }
  harp::vector_double spot_params_1 = psf->FixedCoordParamsXW(x,wave,bundle,gP);
  harp::vector_double spot_params_2 = psf->FixedCoordParamsXW(x,wave,bundle);
  cout << "spot_params_1=" << spot_params_1 << endl;
  cout << "spot_params_2=" << spot_params_2 << endl;
  cout << "diff=" << spot_params_1-spot_params_2 << endl;

  double val0 = psf->PSFValueWithParamsXY(x,y, i, j, spot_params_1, &PosDer, &ParamDer);

  // compute analytic der
  harp::vector_double H(gP.size());
  index=0;
  for(int p=0;p<psf->NPar();p++) {
    int p_size = psf_global_params.Polynomials[p].coeff.size();
    ublas::project(H,ublas::range(index,index+p_size))=ParamDer[p]*ublas::project(gM,ublas::range(index,index+p_size));
    index += p_size;
  }

  for(int k=0;k<int(gP.size());k++) {
    double eps = fabs(gP(k))/1000.;

    harp::vector_double gPp = gP; gPp(k)+=eps/2.;
    harp::vector_double lPp = psf->FixedCoordParamsXW(x,wave,bundle,gPp);
    double valp = psf->PSFValueWithParamsXY(x,y,i,j,lPp,0,0);

    harp::vector_double gPm = gP; gPm(k)-=eps/2.;
    harp::vector_double lPm = psf->FixedCoordParamsXW(x,wave,bundle,gPm);
    double valm = psf->PSFValueWithParamsXY(x,y,i,j,lPm,0,0);
    
    double numDer = (valp-valm)/eps;
    
    cout << k << " P=" << gP(k) << " aD=" << H(k) << " nD=" << numDer 
	 << " diff=" <<  H(k)-numDer << " ratio=" << (H(k)-numDer)/numDer << endl;
  }
  
  }
  
  return 0;
}
