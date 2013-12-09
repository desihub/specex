#include "harp.hpp"
#include "specex_gaussian_analytic_psf.h"
#include "specex_linalg.h"

using namespace std;

specex::GaussPSF::GaussPSF() {
  paramNames.push_back("wxx"); 
  paramNames.push_back("wyy");
  paramNames.push_back("wxy");
}

std::string specex::GaussPSF::Name() const {return "GAUSSIAN";}

size_t specex::GaussPSF::NPar() const { return 3;}

bool specex::GaussPSF::CheckParams(const harp::vector_double &Params) const
{ return (Params(0)*Params(1)-square(Params(2))) > 0;}



double specex::GaussPSF::Profile(const double &X, const double &Y,
			  const harp::vector_double &Params,
			  harp::vector_double *PosDer,
			  harp::vector_double *ParamDer) const
{
  /* parameters (Params) wxx, wyy, wxy */
  double det = Params(0)*Params(1)-Params(2)*Params(2);
  if(det<0) {
    cout << "WARNING in specex::GaussPSF::Profile det=" << det << " return " << PSF_NAN_VALUE << endl;
    //abort();
    return PSF_NAN_VALUE;
  }
  double norm = sqrt(det)/(2*M_PI);
  double exponent = -0.5*(X*X*Params(0)+Y*Y*Params(1)+2*X*Y*Params(2));
  /*
  if(exponent<-40) {
    cout << "WARNING in specex::GaussPSF::Profile exponent=" << exponent << " brought back to -40" << endl;
    exponent=-40;
  }
  */
  if(exponent>40) {
    cout << "WARNING in specex::GaussPSF::Profile exponent=" << exponent << "return "<< PSF_NAN_VALUE << endl;
    return PSF_NAN_VALUE;
  }
  double val = exp(exponent)*norm;
  if (PosDer) // wrong sign on purpose (derivatives w.r.t -X)
    {
      (*PosDer)(0) = (X*Params(0)+Y*Params(2))*val;
      (*PosDer)(1) = (Y*Params(1)+X*Params(2))*val;
    }
  if (ParamDer)
    {
      (*ParamDer)(0) = -0.5*val*(X*X - Params(1)/det);
      (*ParamDer)(1) = -0.5*val*(Y*Y - Params(0)/det);
      (*ParamDer)(2) = -val*(X*Y  + Params(2)/det);
    }
  return val;
}
  
  
void specex::GaussPSF::InitParams(const double &sigma, harp::vector_double &Params)
{
  if (Params.size() < NPar()) Params.resize(NPar());
  Params(0) = 1/(square(sigma));
  Params(1) = Params(0);
  Params(2) = 0;
}
