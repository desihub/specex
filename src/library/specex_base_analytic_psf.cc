#include <vector>
#include <string>
#include <cmath>

#include "harp.hpp"

#include "specex_base_analytic_psf.h"


using namespace std;




/*! 
  For regular users, there is no need to worry about what is going on 
  here: the ImagePSF class handles the needed interfaces.

  AnalyticPSF is a virtual class with actual derived classes that
  implement a given profile. You can choose your PSF type by name using
  ChooseAnalyticPSF. You can also add your own PSF to the list,
  without even altering this file. There is an example of test
  code (in order to make sure that your analytical derivatives
  are OK) at the end of this file. A trick is that the sign
  of the position derivatives has to be inverted, because the user
  worries about derivatives of profile(x_pix-x_star, ...) w.r.t
  x_star.

    The implemented derived classes have an integral of 1 over
  the whole plane. The constraint is integrated into the 
  derivatives w.r.t the profile parameters. This has to be enforced
  if you add new profiles.

  Analytic PSF's are integrated over pixels using Gauss quadrature,
 (as in DAOPHOT, where the relevant coefficients were taken),
  in the PixValue routine. The number of integration steps
  is choosen via a define (NPT)

*/


/* weights and abcissa for gauss integrations (borrowed from DAOPHOT),
   explained in Numerical recipes */

static double Dx[4][4] ={{0.00000000,  0.0,        0.0       , 0.0       },
			 {-0.28867513,  0.28867513, 0.0       , 0.0       },
			 {-0.38729833,  0.00000000, 0.38729833, 0.0       },
			 {-0.43056816, -0.16999052, 0.16999052, 0.43056816}};
static double Wt[4][4]= {{1.00000000,  0.0       , 0.0       , 0.0       },
			 {0.50000000,  0.50000000, 0.0       , 0.0       },
			 {0.27777778,  0.44444444, 0.27777778, 0.0       },
			 {0.17392742,  0.32607258, 0.32607258, 0.17392742}};
 




// number of points (per coordinate) to integrate over a pixel.
#define NPT 3


double specex::AnalyticPSF::PixValue(const double &Xc, const double &Yc,
				     const double &XPix, const double &YPix,
				     const harp::vector_double &Params,
				     harp::vector_double *PosDer,
				     harp::vector_double *ParamDer) const
{
  double xPixCenter = floor(XPix+0.5);
  double yPixCenter = floor(YPix+0.5);
  
  double integrated_dPdx=0;
  double integrated_dPdy=0;
  
  double *dPdx=0;
  double *dPdy=0;
  if(PosDer) {
    dPdx=&((*PosDer)[0]);
    dPdy=&((*PosDer)[1]);
  }
  
  //double *param_der=0;
  //double *integrated_param_der=0; 
  int npar=0;
  harp::vector_double& tmpParamDer = const_cast<AnalyticPSF*>(this)->TmpParamDer;

  if(ParamDer) {
    npar=ParamDer->size();
    tmpParamDer.resize(npar);
    tmpParamDer *= 0; // is there a zero function?
  }
  
  double val = 0;
  for (int ix=0; ix<NPT; ++ix)
    {
      double x = xPixCenter+Dx[NPT-1][ix];
      double wx = Wt[NPT-1][ix];
      for (int iy=0; iy<NPT; ++iy)
	{
	  double y = yPixCenter+Dx[NPT-1][iy];
	  double weight = wx*Wt[NPT-1][iy];
	  double prof = Profile(x-Xc,y-Yc, Params, PosDer, ParamDer);

	  if(prof ==  PSF_NAN_VALUE) return  PSF_NAN_VALUE;
	  
	  val += weight*prof;
	  
	  if (PosDer) {
	    integrated_dPdx += (*dPdx)*weight;
	    integrated_dPdy += (*dPdy)*weight;
	  }
	  if (ParamDer)
	    for (int ipar = 0; ipar<npar; ++ipar) 
	      tmpParamDer[ipar] += weight*(*ParamDer)[ipar];
	}
    }
  if (PosDer) {
    *dPdx=integrated_dPdx;
    *dPdy=integrated_dPdy;
  }
  if (ParamDer) {
    for (int ipar = 0; ipar<npar; ++ipar) 
      (*ParamDer)[ipar] = tmpParamDer[ipar];
  }
#ifdef DEBUG
  cout << "PixValue " 
       << Xc << ' ' << Yc  << ' '
       << XPix << ' ' << YPix << ' ' 
       << Params
       << val << ' ' 
       << endl;
#endif
  return val;
}


#ifdef DEAD

/***************** GaussPSF *************************/

//! gaussian PSF with 3 parameters :wxx, wyy , wxy.
class GaussPSF : public AnalyticPSF {
 public :

  GaussPSF() {paramNames.push_back("wxx"); 
    paramNames.push_back("wyy");paramNames.push_back("wxy");}

  string Name() const {return "GAUSSIAN";}
  unsigned NPar() const { return 3;}
  double Profile(const double &X, const double &Y,
			 const harp::vector_double &Params,
			 harp::vector_double *PosDer = 0,
			 harp::vector_double *ParamGradient = 0) const;
  void InitParamsFromSeeing(const double &Seeing, harp::vector_double &Params) const;


  bool CheckParams(const harp::vector_double &Params) const
  { return (Params(0)*Params(1)-sq(Params(2))) > 0;}
		   


  virtual ~GaussPSF(){}; // warning killer

};
  


double GaussPSF::Profile(const double &X, const double &Y,
			  const harp::vector_double &Params,
			  harp::vector_double *PosDer,
			  harp::vector_double *ParamDer) const
{
  /* parameters (Params) wxx, wyy, wxy */
  double det = Params(0)*Params(1)-Params(2)*Params(2);
  if(det<0) {
    cout << "WARNING in GaussPSF::Profile det=" << det << " return " << PSF_NAN_VALUE << endl;
    //abort();
    return PSF_NAN_VALUE;
  }
  double norm = sqrt(det)/(2*M_PI);
  double exponent = -0.5*(X*X*Params(0)+Y*Y*Params(1)+2*X*Y*Params(2));
  /*
  if(exponent<-40) {
    cout << "WARNING in GaussPSF::Profile exponent=" << exponent << " brought back to -40" << endl;
    exponent=-40;
  }
  */
  if(exponent>40) {
    cout << "WARNING in GaussPSF::Profile exponent=" << exponent << "return "<< PSF_NAN_VALUE << endl;
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
  
  
void GaussPSF::InitParamsFromSeeing(const double &Seeing, harp::vector_double &Params) const
{
  if (Params.Size() < NPar()) Params = harp::vector_double(NPar());
  Params(0) = 1/(sq(Seeing));
  Params(1) = Params(0);
  Params(2) = 0;
}


/************* MoffatPSF *******************/


class MoffatPSFFixedBeta : public AnalyticPSF {
protected :
  const double exponent;
 public :

  MoffatPSFFixedBeta(const double &Beta) : exponent(Beta) 
  {paramNames.push_back("wxx"); 
    paramNames.push_back("wyy");paramNames.push_back("wxy");}

  unsigned NPar() const { return 3;}
  double Profile(const double &X, const double &Y,
			 const harp::vector_double &Params,
			 harp::vector_double *PosDer = 0,
			 harp::vector_double *ParamGradient = 0) const;
  void InitParamsFromSeeing(const double &Seeing, harp::vector_double &Params) const;

  bool CheckParams(const harp::vector_double &Params) const
  { return (Params(0)*Params(1)-sq(Params(2))) > 0;}


  virtual ~MoffatPSFFixedBeta(){}; // warning killer

};


/**************** Moffat20PSF ****************/

class Moffat20PSF : public MoffatPSFFixedBeta {
 public :

  Moffat20PSF() :MoffatPSFFixedBeta(2.0) {};

  string Name() const {return "MOFFAT20";}

  virtual ~Moffat20PSF(){}; // warning killer

};

class Moffat25PSF : public MoffatPSFFixedBeta {
 public :

  Moffat25PSF() :MoffatPSFFixedBeta(2.5) {};

  string Name() const {return "MOFFAT25";}

  virtual ~Moffat25PSF(){}; // warning killer

};

class Moffat30PSF : public MoffatPSFFixedBeta {
 public :

  Moffat30PSF() :MoffatPSFFixedBeta(3.0) {};

  string Name() const {return "MOFFAT30";}

  virtual ~Moffat30PSF(){}; // warning killer

};
  


double MoffatPSFFixedBeta::Profile(const double &X, const double &Y,
			  const harp::vector_double &Params,
			  harp::vector_double *PosDer,
			  harp::vector_double *ParamDer) const
{
  /* parameters (Params) wxx, wyy, wxy */
  double det = Params(0)*Params(1)-Params(2)*Params(2);
  // in principle we should use sqrt(det). fabs saves some nan.
  double norm = sqrt(fabs(det))*(exponent-1)/M_PI;
  double fact =  1./(1.+(X*X*Params(0)+Y*Y*Params(1)+2*X*Y*Params(2)));
  double val = pow(fact,exponent)*norm;
  if (PosDer)
    {
      (*PosDer)(0) = 2*exponent*val*fact*(X*Params(0)+Y*Params(2));
      (*PosDer)(1) = 2*exponent*val*fact*(Y*Params(1)+X*Params(2));
    }
  if (ParamDer)
    {
      (*ParamDer)(0) = -exponent*val*fact*(X*X)+ 0.5*val*Params(1)/det;
      (*ParamDer)(1) = -exponent*val*fact*(Y*Y)+ 0.5*val*Params(0)/det;
      (*ParamDer)(2) = -2*exponent*val*fact*(X*Y) - val*Params(2)/det;
    }
  return val;
}

void MoffatPSFFixedBeta::InitParamsFromSeeing(const double &Seeing, harp::vector_double &Params) const
{
  if (Params.Size() < NPar()) Params = harp::vector_double(NPar());
  Params(0) = sq(0.6547/Seeing);
  Params(1) = Params(0);
  Params(2) = 0;
}

/***********************     MoffatPSF with variable exponent **************/

class MoffatPSF : public AnalyticPSF {

private :
  double startingBeta;
  
 public :

  MoffatPSF(const double StartingBeta = 2.5)
  {
    paramNames.push_back("wxx"); 
    paramNames.push_back("wyy");
    paramNames.push_back("wxy");
    paramNames.push_back("beta");
    startingBeta = StartingBeta;
  }

  string Name() const {return "MOFFAT_VAR_BETA";}

  unsigned NPar() const { return 4;}
  double Profile(const double &X, const double &Y,
			 const harp::vector_double &Params,
			 harp::vector_double *PosDer = 0,
			 harp::vector_double *ParamGradient = 0) const;
  void InitParamsFromSeeing(const double &Seeing, harp::vector_double &Params) const;

  bool CheckParams(const harp::vector_double &Params) const
  { return (Params(0)*Params(1)-sq(Params(2))) > 0;}


  virtual ~MoffatPSF(){}; // warning killer

};


void  MoffatPSF::InitParamsFromSeeing(const double &Seeing, harp::vector_double &Params) const
{
  if (Params.Size() < NPar()) Params = harp::vector_double(NPar());
  Params(0) = sq(0.6547/Seeing);
  Params(1) = Params(0);
  Params(2) = 0;
  Params(3) = startingBeta;
}

double MoffatPSF::Profile(const double &X, const double &Y,
			  const harp::vector_double &Params,
			  harp::vector_double *PosDer,
			  harp::vector_double *ParamDer) const
{
  /* parameters (Params) wxx, wyy, wxy */
  double exponent = Params(3);
  double det = Params(0)*Params(1)-Params(2)*Params(2);
  // in principle we should use sqrt(det). fabs saves some nan.
  double norm = sqrt(fabs(det))*(exponent-1)/M_PI;
  double fact =  1./(1.+(X*X*Params(0)+Y*Y*Params(1)+2*X*Y*Params(2)));
  double val = pow(fact,exponent)*norm;
  if (PosDer)
    {
      (*PosDer)(0) = 2*exponent*val*fact*(X*Params(0)+Y*Params(2));
      (*PosDer)(1) = 2*exponent*val*fact*(Y*Params(1)+X*Params(2));
    }
  if (ParamDer)
    {
      (*ParamDer)(0) = -exponent*val*fact*(X*X)+ 0.5*val*Params(1)/det;
      (*ParamDer)(1) = -exponent*val*fact*(Y*Y)+ 0.5*val*Params(0)/det;
      (*ParamDer)(2) = -2*exponent*val*fact*(X*Y) - val*Params(2)/det;
      (*ParamDer)(3) = (log(fact)+1./(exponent-1))*val;
      
    }
  return val;
}



#endif

