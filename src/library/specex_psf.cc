#include <iomanip>
#include <cmath>
#include <string>

#include "harp.hpp"

#include "specex_psf.h"
#include "specex_base_analytic_psf.h"

using namespace std;


specex::PSF::PSF() {
  analyticPSF = NULL;
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
  return analyticPSF->NPar();
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
				     harp::vector_double *PosDer, harp::vector_double *ParamDer, double *AnalyticValue) const {
  
  //if (PosDer)   PosDer->Zero();
  //if (ParamDer) ParamDer->Zero();  
  return analyticPSF->PixValue(Xc,Yc,IPix, JPix, Params, PosDer, ParamDer);
}

//! Access to the current PSF pixels.
double specex::PSF::PSFValue(const double &Xc, const double &Yc, 
			   const int IPix, const int JPix,
			   harp::vector_double *PosDer, harp::vector_double *ParamDer, double *AnalyticValue) const {
  return PSFValueWithParams(Xc,Yc,IPix,JPix,FixedCoordParams(Xc,Yc), PosDer, ParamDer, AnalyticValue);
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
  if(analyticPSF->Name() == "GAUSSHERMITE") return true;
  return false;
}
