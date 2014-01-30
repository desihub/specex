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
  
  harp::vector_double tmpPosDer;
  if(PosDer) {
    tmpPosDer = boost::numeric::ublas::zero_vector<double>(2);
  }
  harp::vector_double tmpParamDer;
  int npar=0;
  if(ParamDer) {
    npar = ParamDer->size();
    tmpParamDer = boost::numeric::ublas::zero_vector<double>(npar);
  }
  
  double val = 0;
  for (int ix=0; ix<NPT; ++ix)
    {
      double x = xPixCenter+Dx[NPT-1][ix] - Xc;
      double wx = Wt[NPT-1][ix];
      for (int iy=0; iy<NPT; ++iy)
	{
	  double y = yPixCenter+Dx[NPT-1][iy] -Yc;
	  double weight = wx*Wt[NPT-1][iy];
	  double prof = Profile(x,y,Params,PosDer,ParamDer);

	  if(prof ==  PSF_NAN_VALUE) return  PSF_NAN_VALUE;
	  
	  val += weight*prof;
	  
	  if (PosDer) 
	    tmpPosDer += weight*(*PosDer);
	  
	  if (ParamDer)
	    tmpParamDer += weight*(*ParamDer);
	  
	}
    }
  
  if (PosDer)
    *PosDer = tmpPosDer;
  
  if (ParamDer)
    *ParamDer = tmpParamDer; 
  
  
  return val;
}



specex::PSF::PSF() {
  name = "unknown";
  hSizeX = hSizeY = 12;

#ifdef EXTERNAL_TAIL
  //r_tail_amplitude = 0;
  r_tail_core_size = 1;
  r_tail_x_scale   = 1;
  r_tail_y_scale   = 1;
#ifdef EXTERNAL_Y_TAIL
  y_tail_amplitude = 0;
  y_tail_core_size = 1;
  y_tail_power_law_index = 1;
  y_tail_sigma_x =1.1;
#endif
#endif
#ifdef CONTINUUM
  continuum_sigma_x = 1;
#endif
}

#ifdef EXTERNAL_TAIL
#define NX_TAIL_PROFILE 1000
#define NY_TAIL_PROFILE 8000
#define TAIL_OVERSAMPLING 2.

void specex::PSF::ComputeTailProfile() {
  
#pragma omp critical 
  {
  SPECEX_INFO("specex::PSF::ComputeTailProfile ...");
  
  r_tail_profile.resize(NX_TAIL_PROFILE,NY_TAIL_PROFILE); // hardcoded
  
  double rc2 = square(r_tail_core_size);
  
  for(int j=0;j<NY_TAIL_PROFILE;j++) {
    for(int i=0;i<NX_TAIL_PROFILE;i++) {
      double r2 = square(i/TAIL_OVERSAMPLING*r_tail_x_scale)+square(j/TAIL_OVERSAMPLING*r_tail_y_scale);
	
      //r_tail_profile(i,j) = r2/(rc2+r2)*pow(rc2+r2,-1.5/2.);
      r_tail_profile(i,j) = r2/(rc2+r2)*pow(rc2+r2,-2./2.);
      
    }
  }
  
#ifdef EXTERNAL_Y_TAIL
  y_tail_profile.resize(NX_TAIL_PROFILE,NY_TAIL_PROFILE); // hardcoded
  
  double yc2 = square(y_tail_core_size);
  
  for(int j=0;j<NY_TAIL_PROFILE;j++) {
    
    double y2 = square(j/TAIL_OVERSAMPLING);
    
    double yprof = y2/(yc2+y2)*pow(yc2+y2,-y_tail_power_law_index/2.);
    
    for(int i=0;i<NX_TAIL_PROFILE;i++) {
      
      y_tail_profile(i,j) = yprof*exp(-0.5*square(i/TAIL_OVERSAMPLING/y_tail_sigma_x));
    }
  }
  
#endif
  
  SPECEX_INFO("specex::PSF::ComputeTailProfile done");
   
  }
  
}


double specex::PSF::TailValueW(const double& wavelength, 
			      const double& dx, const double &dy, 
			      harp::vector_double* derivative_r_tail_amplitude,
			      double* derivative_y_tail_amplitude) const {
  
  if(r_tail_profile.n_rows()==0) 
    const_cast<specex::PSF*>(this)->ComputeTailProfile();

  
  int di = int(fabs(dx*TAIL_OVERSAMPLING)+0.5);
  int dj = int(fabs(dy*TAIL_OVERSAMPLING)+0.5);
  if(di>NX_TAIL_PROFILE || dj>NY_TAIL_PROFILE) return 0.;

  double r_prof = r_tail_profile(di,dj);
  
  harp::vector_double monomials = RTailAmplitudePol.Monomials(wavelength);
#ifndef EXPONENTIAL_TAIL_AMPLITUDE
  if(derivative_r_tail_amplitude) *derivative_r_tail_amplitude = r_prof*monomials;  
  double res =  specex::dot(monomials,RTailAmplitudePol.coeff)*r_prof;
#else
  double exponent = specex::dot(monomials,RTailAmplitudePol.coeff);
  if(exponent<-10)  exponent=-10;
  else if (exponent>10) exponent=10;  
  double res =  exp(exponent)*r_prof;
  if(derivative_r_tail_amplitude) *derivative_r_tail_amplitude = res*monomials; 
#endif


#ifdef EXTERNAL_Y_TAIL
  double y_prof = y_tail_profile(di,dj);
  if(derivative_y_tail_amplitude) *derivative_y_tail_amplitude = y_prof;
  res += y_tail_amplitude*y_prof;
#endif

  return res;

}


#ifndef EXTERNAL_Y_TAIL
double specex::PSF::TailValueA(const double& r_tail_amplitude, 
			       const double& dx, const double &dy) const {
  
  if(r_tail_profile.n_rows()==0) 
    const_cast<specex::PSF*>(this)->ComputeTailProfile();

  int di = int(fabs(dx*TAIL_OVERSAMPLING)+0.5);
  int dj = int(fabs(dy*TAIL_OVERSAMPLING)+0.5);
  if(di>NX_TAIL_PROFILE || dj>NY_TAIL_PROFILE) return 0.;
#ifndef EXPONENTIAL_TAIL_AMPLITUDE  
  return r_tail_amplitude*r_tail_profile(di,dj);
#else
  double exponent = r_tail_amplitude;
  if(exponent<-10)  exponent=-10;
  else if (exponent>10) exponent=10; 
  return exp(exponent)*r_tail_profile(di,dj);
#endif
}
#endif

#endif

int specex::PSF::BundleNFitPar(int bundle_id) const {
  std::map<int,PSF_Params>::const_iterator it = ParamsOfBundles.find(bundle_id);
  if(it==ParamsOfBundles.end()) SPECEX_ERROR("no such bundle #" << bundle_id);
  const std::vector<Legendre2DPol_p>& P=it->second.FitParPolXW;
  int n=0;
  for(size_t p=0;p<P.size();p++)
    n += P[p]->coeff.size();
  return n;
}
int specex::PSF::BundleNAllPar(int bundle_id) const {
  std::map<int,PSF_Params>::const_iterator it = ParamsOfBundles.find(bundle_id);
  if(it==ParamsOfBundles.end()) SPECEX_ERROR("no such bundle #" << bundle_id);
  const std::vector<Legendre2DPol_p>& P=it->second.AllParPolXW;
  int n=0;
  for(size_t p=0;p<P.size();p++)
    n += P[p]->coeff.size();
  return n;
}

int specex::PSF::TracesNPar() const {
  int n=0;
  for(std::map<int,Trace>::const_iterator it=FiberTraces.begin(); it!=FiberTraces.end(); ++it) {
    n += it->second.X_vs_W.coeff.size();
    n += it->second.Y_vs_W.coeff.size();
  }
  return n;
}
specex::PSF::~PSF() {
}

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

const specex::Trace& specex::PSF::GetTrace(int fiber) const {
  std::map<int,specex::Trace>::const_iterator it = FiberTraces.find(fiber);
  if(it == FiberTraces.end()) SPECEX_ERROR("No trace for fiber " << fiber);
  return it->second;
}

specex::Trace& specex::PSF::GetTrace(int fiber)  {
  std::map<int,specex::Trace>::iterator it = FiberTraces.find(fiber);
  if(it == FiberTraces.end()) SPECEX_ERROR("No trace for fiber " << fiber);
  return it->second;
}

void specex::PSF::AddTrace(int fiber) {
  std::map<int,specex::Trace>::iterator it = FiberTraces.find(fiber);
  if(it != FiberTraces.end()) {
    SPECEX_WARNING("Trace for fiber " << fiber << " already exists");
  }else{
    FiberTraces[fiber]=specex::Trace();
    LoadXYPol();
  }
}

void specex::PSF::LoadXYPol() {
#pragma omp critical
  {
  XPol.clear();
  YPol.clear();
  for(std::map<int,specex::Trace>::iterator it = FiberTraces.begin(); it != FiberTraces.end(); ++it) {
    XPol[it->first] = &(it->second.X_vs_W);
    YPol[it->first] = &(it->second.Y_vs_W);
  }
  }
}

double specex::PSF::Xccd(int fiber, const double& wave) const {
  std::map<int,specex::Legendre1DPol*>::const_iterator it = XPol.find(fiber);
  if(it==XPol.end()) {
    const_cast<specex::PSF*>(this)->LoadXYPol();
    it = XPol.find(fiber);
    if(it == XPol.end()) SPECEX_ERROR("No trace for fiber " << fiber);
  }
  return it->second->Value(wave);
}

double specex::PSF::Yccd(int fiber, const double& wave) const {
  std::map<int,specex::Legendre1DPol*>::const_iterator it = YPol.find(fiber);
  if(it==YPol.end()) {
    const_cast<specex::PSF*>(this)->LoadXYPol();
    it = YPol.find(fiber);
    if(it == YPol.end()) SPECEX_ERROR("No trace for fiber " << fiber);
  }
  return it->second->Value(wave);
}  


//! Access to the current PSF, with user provided Params.
double specex::PSF::PSFValueWithParamsXY(const double &Xc, const double &Yc, 
				     const int IPix, const int JPix,
				     const harp::vector_double &Params,
				     harp::vector_double *PosDer, harp::vector_double *ParamDer) const {
  
  return PixValue(Xc,Yc,IPix, JPix, Params, PosDer, ParamDer); 
}

//! Access to the current PSF, with user provided Params.
double specex::PSF::PSFValueWithParamsFW(const int fiber, const double &wave,
				     const int IPix, const int JPix,
				     const harp::vector_double &Params,
				     harp::vector_double *PosDer, harp::vector_double *ParamDer) const {
  
  double X=Xccd(fiber,wave);
  double Y=Yccd(fiber,wave);
  return PixValue(X,Y,IPix, JPix, Params, PosDer, ParamDer);
}

//! Access to the current PSF 
double specex::PSF::PSFValueFW(const int fiber, const double &wave,
			       const int IPix, const int JPix, int bundle_id,
			       harp::vector_double *PosDer, harp::vector_double *ParamDer) const {
  double X=Xccd(fiber,wave);
  double Y=Yccd(fiber,wave);
  return PSFValueWithParamsXY(X,Y,IPix,JPix,AllLocalParamsXW(X,wave,bundle_id), PosDer, ParamDer);
}




harp::vector_double specex::PSF::AllLocalParamsXW(const double &X, const double &wave, int bundle_id) const {
  
  std::map<int,PSF_Params>::const_iterator it = ParamsOfBundles.find(bundle_id);
  if(it==ParamsOfBundles.end()) SPECEX_ERROR("no such bundle #" << bundle_id);
  const std::vector<Legendre2DPol_p>& P=it->second.AllParPolXW;
  
  harp::vector_double params(P.size());
  for (size_t k =0; k < P.size(); ++k)
    params(k) = P[k]->Value(X,wave);
  
  return params;
}

harp::vector_double specex::PSF::AllLocalParamsFW(const int fiber, const double &wave, int bundle_id) const {
  double X=Xccd(fiber,wave); 
  return AllLocalParamsXW(X,wave,bundle_id);
}

harp::vector_double specex::PSF::FitLocalParamsXW(const double &X, const double &wave, int bundle_id) const {
  
  std::map<int,PSF_Params>::const_iterator it = ParamsOfBundles.find(bundle_id);
  if(it==ParamsOfBundles.end()) SPECEX_ERROR("no such bundle #" << bundle_id);
  const std::vector<Legendre2DPol_p>& P=it->second.FitParPolXW;
  
  harp::vector_double params(P.size());
  for (size_t k =0; k < P.size(); ++k)
    params(k) = P[k]->Value(X,wave);
  return params;
}

harp::vector_double specex::PSF::FitLocalParamsFW(const int fiber, const double &wave, int bundle_id) const {
  double X=Xccd(fiber,wave); 
  return FitLocalParamsXW(X,wave,bundle_id);
}


harp::vector_double specex::PSF::AllLocalParamsXW_with_AllBundleParams(const double &X, const double &wave, int bundle_id, const harp::vector_double& ForThesePSFParams) const {
  
  std::map<int,PSF_Params>::const_iterator it = ParamsOfBundles.find(bundle_id);
  if(it==ParamsOfBundles.end()) SPECEX_ERROR("no such bundle #" << bundle_id);
  const std::vector<Legendre2DPol_p>& P=it->second.AllParPolXW;
  
  harp::vector_double params(LocalNAllPar());
  
  if(BundleNAllPar(bundle_id)>ForThesePSFParams.size()) SPECEX_ERROR("VaryingCoordNPar(bundle_id)<=ForThesePSFParams.size()");
  
  int index=0;
  for (size_t k =0; k < P.size(); ++k) {
    size_t c_size = P[k]->coeff.size();
    params(k)=specex::dot(ublas::project(ForThesePSFParams,ublas::range(index,index+c_size)),P[k]->Monomials(X,wave));
    index += c_size;
  }
  return params;
}

  
harp::vector_double specex::PSF::AllLocalParamsFW_with_AllBundleParams(const int fiber, const double &wave, int bundle_id, const harp::vector_double& ForThesePSFParams) const {
  double X=Xccd(fiber,wave); 
  return AllLocalParamsXW_with_AllBundleParams(X,wave,bundle_id,ForThesePSFParams);
}


harp::vector_double specex::PSF::AllLocalParamsXW_with_FitBundleParams(const double &X, const double &wave, int bundle_id, const harp::vector_double& ForThesePSFParams) const {
  
  
  harp::vector_double params(LocalNAllPar());
  
  std::map<int,PSF_Params>::const_iterator it = ParamsOfBundles.find(bundle_id);
  if(it==ParamsOfBundles.end()) SPECEX_ERROR("no such bundle #" << bundle_id);
  const std::vector<Legendre2DPol_p>& AP=it->second.AllParPolXW;
  const std::vector<Legendre2DPol_p>& FP=it->second.FitParPolXW;
  
  // whe need to find which param is fixed and which is not
  size_t fk=0;
  int index=0;
  for (size_t ak =0; ak < AP.size(); ++ak) { // loop on all params
    const Legendre2DPol_p APk = AP[ak];
    size_t c_size = APk->coeff.size();
    if(fk<FP.size() && APk==FP[fk]) { // this is a fit param because the addresses are the same
      params(ak)=specex::dot(ublas::project(ForThesePSFParams,ublas::range(index,index+c_size)),APk->Monomials(X,wave));
      
      //SPECEX_INFO("DEBUG all param " << ak << " and fit = " << fk << " are the same, param val =  " << params(ak));

      index += c_size;
      // change free param index for next iteration
      fk++;
    }else{ // this not a free param
      params(ak)=APk->Value(X,wave);
      //SPECEX_INFO("DEBUG all param " << ak << " is not it fit, param val = " << params(ak));
    }
  }
  return params;
}
harp::vector_double specex::PSF::AllLocalParamsFW_with_FitBundleParams(const int fiber, const double &wave, int bundle_id, const harp::vector_double& ForThesePSFParams) const {
  double X=Xccd(fiber,wave); 
  return AllLocalParamsXW_with_FitBundleParams(X,wave,bundle_id,ForThesePSFParams);
}

harp::vector_double specex::PSF::FitLocalParamsXW_with_FitBundleParams(const double &X, const double &wave, int bundle_id, const harp::vector_double& ForThesePSFParams) const {
  
  std::map<int,PSF_Params>::const_iterator it = ParamsOfBundles.find(bundle_id);
  if(it==ParamsOfBundles.end()) SPECEX_ERROR("no such bundle #" << bundle_id);
  const std::vector<Legendre2DPol_p>& P=it->second.FitParPolXW;
  
  harp::vector_double params(P.size());
  
  if(BundleNFitPar(bundle_id)>ForThesePSFParams.size()) SPECEX_ERROR("VaryingCoordNPar(bundle_id)<=ForThesePSFParams.size()");
  
  int index=0;
  for (size_t k =0; k < P.size(); ++k) {
    size_t c_size = P[k]->coeff.size();
    params(k)=specex::dot(ublas::project(ForThesePSFParams,ublas::range(index,index+c_size)),P[k]->Monomials(X,wave));
    index += c_size;
  }
  return params;
}

  
harp::vector_double specex::PSF::FitLocalParamsFW_with_FitBundleParams(const int fiber, const double &wave, int bundle_id, const harp::vector_double& ForThesePSFParams) const {
  double X=Xccd(fiber,wave); 
  return FitLocalParamsXW_with_FitBundleParams(X,wave,bundle_id,ForThesePSFParams);
}


bool specex::PSF::IsLinear() const {
  if(Name() == "GAUSSHERMITE") return true;
  return false;
}

void specex::PSF::WriteFits(const std::string& filename, int first_hdu) const {  
  fitsfile * fp;  
  harp::fits::create ( fp, filename );
  WriteFits(fp,first_hdu);
  harp::fits::close ( fp );
  
  SPECEX_INFO("wrote psf in " << filename);
}
    
void specex::PSF::ReadFits(const std::string& filename, int first_hdu)  {  
  fitsfile * fp;  
  harp::fits::open_read ( fp, filename );
  ReadFits(fp,first_hdu);
  harp::fits::close ( fp );
  
  SPECEX_INFO("read psf in " << filename);
}
    
