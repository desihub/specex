#include <iomanip>
#include <cmath>
#include <string>

#include "harp.hpp"

#include "specex_psf.h"
//#include "specex_base_analytic_psf.h"
#include "specex_message.h"

using namespace std;


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


double specex::PSF::PixValue(const double &Xc, const double &Yc,
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
  harp::vector_double& tmpParamDer = const_cast<specex::PSF*>(this)->TmpParamDer;

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
  
  return val;
}



specex::PSF::PSF() {
  name = "unknown";
  hSizeX = hSizeY = 12;
  verbose = false;
  
  /*
  if(psffile_name!="") {
    bool ok = read(psffile_name);
    if(!ok) abort();
  }
  */
}

int specex::PSF::FixedCoordNPar() const {
  return NPar();
}

int specex::PSF::VaryingCoordNPar() const {
  int n=0;
  for(size_t p=0;p<Params.size();p++)
    n += Params[p].coeff.size();
  return n;
}

int specex::PSF::TracesNPar() const {
  int n=0;
  for(std::map<int,Trace>::const_iterator it=FiberTraces.begin(); it!=FiberTraces.end(); ++it) {
    n += it->second.X_vs_lW.coeff.size();
    n += it->second.Y_vs_lW.coeff.size();
  }
  return n;
}
specex::PSF::~PSF() {
}


/*
void specex::PSF::write(ostream &os) const {
  os << setprecision(10);
  os << "BeginPSF" << endl;
  os << analyticPSF->Name() << endl;
  if(analyticPSF->Name()=="GAUSSHERMITE") {
    os << "GAUSSHERMITE_SPECIFIC_PARAMETERS " << ((GaussHermitePSF*)analyticPSF)->Degree() << " "<< ((GaussHermitePSF*)analyticPSF)->sigma << endl;
  }

  os << Params.size() << " " << hSizeX << " " << hSizeY << endl;
  for(size_t p=0; p<Params.size(); p++) {
    const Legendre2DPol& pol=Params[p];
    pol.write(os);
  }
  os << FiberTraces.size() << endl;
  for(map<int,Spec2DTrace>::const_iterator it=FiberTraces.begin();
      it !=FiberTraces.end(); ++it) {
    it->second.write(os);
  }
  os << "EndPSF" << endl;
}
#define ERROR_MESSAGE {cout << "ERROR specex::PSF::read failed" << endl; abort();}

bool specex::PSF::read(istream& is) {
   string label;
  if(! (is >> label)) ERROR_MESSAGE;
  if(! (label=="BeginPSF")) ERROR_MESSAGE;
  
  if(! (is >> label)) ERROR_MESSAGE;
  // choose analytic psf here
  analyticPSF = ChooseAnalyticPSF(label);

  if(label=="GAUSSHERMITE") {
    double sigma; int degree;
    if(! (is >> label >> degree >> sigma)) ERROR_MESSAGE;
    if(label != "GAUSSHERMITE_SPECIFIC_PARAMETERS") ERROR_MESSAGE;
    ((GaussHermitePSF*)analyticPSF)->SetDegree(degree);
    ((GaussHermitePSF*)analyticPSF)->sigma = sigma; 
  }


  int npar;
  if(! (is >> npar >> hSizeX >> hSizeY)) ERROR_MESSAGE;
  Params.clear();
  for(int p=0;p<npar;p++) {
    Legendre2DPol pol;
    if(! pol.read(is)) ERROR_MESSAGE;
    Params.push_back(pol);
  }
  int ntraces;
  FiberTraces.clear();
  if(! (is >> ntraces)) ERROR_MESSAGE;
  //cout << "INFO specex::PSF::read starting reading " << ntraces << " traces" << endl;
  for(int t=0;t<ntraces;t++) {
    Spec2DTrace trace;
    if(! trace.read(is)) ERROR_MESSAGE;
    FiberTraces[trace.fiber]=trace;
    //cout << "INFO specex::PSF::read read trace for fiber " << trace.fiber << endl;
  }
  if(! (is >> label)) ERROR_MESSAGE;
  if(! (label=="EndPSF")) ERROR_MESSAGE;
  cout << "INFO specex::PSF::read ok" << endl;
  return true;
}

#undef ERROR_MESSAGE


void specex::PSF::write(const std::string &FileName) const {
  ofstream os(FileName.c_str());
  write(os);
  os.close();
}

bool specex::PSF::read(const std::string &FileName) {
  ifstream is(FileName.c_str());
  bool ok = read(is);
  is.close();
  return ok;
}
 */
 
void specex::PSF::StampLimits(const double &X, const double &Y,
			    int &BeginI, int &EndI,
			    int &BeginJ, int &EndJ) const {
  
  int iPix = int(floor(X+0.5));
  int jPix = int(floor(Y+0.5));
  BeginI = iPix-hSizeX;
  BeginJ = jPix-hSizeY;
  EndI = iPix+hSizeX+1;
  EndJ = jPix+hSizeY+1; 
}

//! Access to the current PSF, with user provided Params.
double specex::PSF::PSFValueWithParams(const double &Xc, const double &Yc, 
				     const int IPix, const int JPix,
				     const harp::vector_double &Params,
				     harp::vector_double *PosDer, harp::vector_double *ParamDer) const {
  
  //if (PosDer)   PosDer->Zero();
  //if (ParamDer) ParamDer->Zero();  
  return PixValue(Xc,Yc,IPix, JPix, Params, PosDer, ParamDer);
}

//! Access to the current PSF pixels.
double specex::PSF::PSFValue(const double &Xc, const double &Yc, 
			   const int IPix, const int JPix,
			   harp::vector_double *PosDer, harp::vector_double *ParamDer) const {
  return PSFValueWithParams(Xc,Yc,IPix,JPix,FixedCoordParams(Xc,Yc), PosDer, ParamDer);
}


harp::vector_double specex::PSF::FixedCoordParams(const double &X, const double &Y) const {
  harp::vector_double params(Params.size());
  for (size_t k =0; k < Params.size(); ++k)
    params(k) = Params[k].Value(X,Y);
  return params;
}

harp::vector_double specex::PSF::FixedCoordParams(const double &X, const double &Y, const harp::vector_double& ForThesePSFParams) const {
  harp::vector_double params(Params.size());
  
  if(VaryingCoordNPar()>ForThesePSFParams.size()) abort();
  
  int index=0;
  for (size_t k =0; k < Params.size(); ++k) {
    int nc=Params[k].coeff.size();
    harp::vector_double coeff(nc);
    for(int c=0;c<nc;c++,index++)
      coeff(c)=ForThesePSFParams(index);
    params(k)=specex::dot(coeff,Params[k].Monomials(X,Y));
  }
  return params;
}
  
bool specex::PSF::IsLinear() const {
  if(Name() == "GAUSSHERMITE") return true;
  return false;
}

