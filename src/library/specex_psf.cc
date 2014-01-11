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
    dPdx=&((*PosDer)(0));
    dPdy=&((*PosDer)(1));
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
	      tmpParamDer(ipar) += weight*(*ParamDer)(ipar);
	}
    }
  if (PosDer) {
    *dPdx=integrated_dPdx;
    *dPdy=integrated_dPdy;
  }
  if (ParamDer) {
    for (int ipar = 0; ipar<npar; ++ipar) 
      (*ParamDer)(ipar) = tmpParamDer(ipar);
  }
  
  return val;
}



specex::PSF::PSF() {
  name = "unknown";
  hSizeX = hSizeY = 12;

#ifdef EXTERNAL_TAIL
  r_tail_amplitude = 0;
  r_tail_core_size = 1;
  r_tail_x_scale   = 1;
  r_tail_y_scale   = 1;
  y_tail_amplitude = 0;
  y_tail_core_size = 1;
  y_tail_power_law_index = 1;
  y_tail_sigma_x =1.1;
#endif
}

#ifdef EXTERNAL_TAIL
#define NX_TAIL_PROFILE 1000
#define NY_TAIL_PROFILE 8000
#define TAIL_OVERSAMPLING 2.

double specex::PSF::TailValue(const double& dx, const double &dy,double* derivative_r_tail_amplitude,double* derivative_y_tail_amplitude) const {
  if(r_tail_amplitude==0 && y_tail_amplitude==0 && derivative_r_tail_amplitude==0) return 0;
  
  
  if(r_tail_profile.n_rows()==0) {
    SPECEX_INFO("specex::PSF::TailValue computing profile images ...");
    
    specex::PSF* psf = const_cast<specex::PSF*>(this);
    
    psf->r_tail_profile.resize(NX_TAIL_PROFILE,NY_TAIL_PROFILE); // hardcoded
    psf->y_tail_profile.resize(NX_TAIL_PROFILE,NY_TAIL_PROFILE); // hardcoded

    

    for(int j=0;j<NY_TAIL_PROFILE;j++) {
      for(int i=0;i<NX_TAIL_PROFILE;i++) {
	psf->r_tail_profile(i,j) = 1./(square(r_tail_core_size)+square(i/TAIL_OVERSAMPLING*r_tail_x_scale)+square(j/TAIL_OVERSAMPLING*r_tail_y_scale));
      }
    }
    for(int j=0;j<NY_TAIL_PROFILE;j++) {
      for(int i=0;i<NX_TAIL_PROFILE;i++) {
	psf->y_tail_profile(i,j) = pow(square(y_tail_core_size)+square(j/TAIL_OVERSAMPLING),-y_tail_power_law_index/2.)*exp(-0.5*square(i/TAIL_OVERSAMPLING/y_tail_sigma_x));
      }
    }
    SPECEX_INFO("specex::PSF::TailValue computing profile images done");
  }

  int di = int(fabs(dx*TAIL_OVERSAMPLING)+0.5);
  int dj = int(fabs(dy*TAIL_OVERSAMPLING)+0.5);
  if(di>NX_TAIL_PROFILE || dj>NY_TAIL_PROFILE) return 0.;
  
  double r_prof = r_tail_profile(di,dj);
  if(derivative_r_tail_amplitude) *derivative_r_tail_amplitude = r_prof;

  double y_prof = y_tail_profile(di,dj);
  if(derivative_y_tail_amplitude) *derivative_y_tail_amplitude = y_prof;
 
  return r_tail_amplitude*r_prof + y_tail_amplitude*y_prof;
}

#endif

int specex::PSF::BundleNPar(int bundle_id) const {
  std::map<int,PSF_Params>::const_iterator it = ParamsOfBundles.find(bundle_id);
  if(it==ParamsOfBundles.end()) SPECEX_ERROR("no such bundle #" << bundle_id);
  const std::vector<Legendre2DPol>& P=it->second.Polynomials;
  int n=0;
  for(size_t p=0;p<P.size();p++)
    n += P[p].coeff.size();
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
  XPol.clear();
  YPol.clear();
  for(std::map<int,specex::Trace>::iterator it = FiberTraces.begin(); it != FiberTraces.end(); ++it) {
    XPol[it->first] = &(it->second.X_vs_W);
    YPol[it->first] = &(it->second.Y_vs_W);
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
  return PSFValueWithParamsXY(X,Y,IPix,JPix,LocalParamsXW(X,wave,bundle_id), PosDer, ParamDer);
}




harp::vector_double specex::PSF::LocalParamsXW(const double &X, const double &wave, int bundle_id) const {
  
  std::map<int,PSF_Params>::const_iterator it = ParamsOfBundles.find(bundle_id);
  if(it==ParamsOfBundles.end()) SPECEX_ERROR("no such bundle #" << bundle_id);
  const std::vector<Legendre2DPol>& P=it->second.Polynomials;
  
  harp::vector_double params(P.size());
  for (size_t k =0; k < P.size(); ++k)
    params(k) = P[k].Value(X,wave);
  return params;
}

harp::vector_double specex::PSF::LocalParamsFW(const int fiber, const double &wave, int bundle_id) const {
  double X=Xccd(fiber,wave); 
  return LocalParamsXW(X,wave,bundle_id);
}
harp::vector_double specex::PSF::LocalParamsXW(const double &X, const double &wave, int bundle_id, const harp::vector_double& ForThesePSFParams) const {
  
  std::map<int,PSF_Params>::const_iterator it = ParamsOfBundles.find(bundle_id);
  if(it==ParamsOfBundles.end()) SPECEX_ERROR("no such bundle #" << bundle_id);
  const std::vector<Legendre2DPol>& P=it->second.Polynomials;
  
  harp::vector_double params(LocalNPar());
  
  if(BundleNPar(bundle_id)>ForThesePSFParams.size()) SPECEX_ERROR("VaryingCoordNPar(bundle_id)<=ForThesePSFParams.size()");
  
  int index=0;
  for (size_t k =0; k < P.size(); ++k) {
    size_t c_size = P[k].coeff.size();
    params(k)=specex::dot(ublas::project(ForThesePSFParams,ublas::range(index,index+c_size)),P[k].Monomials(X,wave));
    index += c_size;
  }
  return params;
}
  
harp::vector_double specex::PSF::LocalParamsFW(const int fiber, const double &wave, int bundle_id, const harp::vector_double& ForThesePSFParams) const {
  double X=Xccd(fiber,wave); 
  return LocalParamsXW(fiber,wave,bundle_id,ForThesePSFParams);
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
    
