#ifndef SPECEX_PSF__H
#define SPECEX_PSF__H

#include <vector>
#include <string>
#include <map>

#include "specex_legendre.h"
#include "specex_trace.h"
#include "specex_linalg.h"
#include "specex_image_data.h"

#define PSF_NAN_VALUE 9999999

#define EXTERNAL_TAIL

namespace specex {

  
  
  class Prior { 
  public :
    virtual double Chi2(const double& x) const = 0;
    virtual double hdChi2dx(const double& x) const = 0; // - 1/2 d(Chi2)/dx
    virtual double hd2Chi2dx2(const double& x) const = 0; // 1/2 d2(Chi2)/dx2
  };

  
  class GaussianPrior : public Prior {
  public :
    double mean;
    double sigma;
    GaussianPrior(const double& imean=0, const double& isigma=0) :
      mean(imean),
      sigma(isigma)
	{
	}
      
      double Chi2(const double& x)  const { return square((mean-x)/sigma);}
      double hdChi2dx(const double& x)  const {return (mean-x)/square(sigma);} 
      double hd2Chi2dx2(const double& x)  const {return 1./square(sigma);}
  };
  
  class PSF_Params  {

    friend class boost::serialization::access;

  public :

    

    // those are polynomials of x_cdd and lambda to allow continuous variation in CCD and at the same time
    // and easy projection per fiber for subsequent use
    
    std::vector<Pol_p> AllParPolXW; // all the parameters of the PSF , varying with x_ccd and lambda
    std::vector<Pol_p> FitParPolXW; // subsamples of parameters that are fit
    

    int bundle_id;
    int fiber_min; // first fiber this set of params applies to
    int fiber_max; // last fiber this set of params applies to

    double chi2; // fit results
    double chi2_in_core; // fit results
    int ndata;
    int ndata_in_core;
    int nparams;
    int fit_status; // -1=no fit performed 0=ok 1=cholesky error 2=no convergence 3=nan in fit
    int nspots_in_fit;

#define CONTINUUM

#ifdef CONTINUUM
    Legendre1DPol ContinuumPol;
    double continuum_sigma_x;
#endif

  PSF_Params() : 
    bundle_id(0), fiber_min(0), fiber_max(0), 
      chi2(0),ndata(0),fit_status(-1), nspots_in_fit(0)
#ifdef CONTINUUM
, continuum_sigma_x(1)
#endif
      {};
  
  private :

    template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
      ar & BOOST_SERIALIZATION_NVP(AllParPolXW);
      ar & BOOST_SERIALIZATION_NVP(bundle_id);
      ar & BOOST_SERIALIZATION_NVP(fiber_min);
      ar & BOOST_SERIALIZATION_NVP(fiber_max);
      ar & BOOST_SERIALIZATION_NVP(chi2);
      ar & BOOST_SERIALIZATION_NVP(ndata);
      ar & BOOST_SERIALIZATION_NVP(chi2_in_core);
      ar & BOOST_SERIALIZATION_NVP(ndata_in_core);
      ar & BOOST_SERIALIZATION_NVP(nparams);
      ar & BOOST_SERIALIZATION_NVP(fit_status);
      ar & BOOST_SERIALIZATION_NVP(nspots_in_fit);

#ifdef CONTINUUM
      ar & BOOST_SERIALIZATION_NVP(ContinuumPol);
      ar & BOOST_SERIALIZATION_NVP(continuum_sigma_x);
#endif 

    }
  };

  class PSF {

    friend class boost::serialization::access;

    // AnalyticPSF* analyticPSF;
    // BEGIN WAS BEFORE IN ANALYTIC PSF
    // ==================================================

  
 




#ifdef EXTERNAL_TAIL
    
  public :
    double TailProfile(const double& dx, const double &dy, const harp::vector_double &Params, bool full_calculation=false) const;
    
  protected :

    specex::image_data r_tail_profile; // to go much faster    
    double r2_tail_core_size;
    double r_tail_x_scale;
    double r_tail_y_scale;
    double r_tail_power_law_index;
    
    bool r_tail_profile_must_be_computed;
    void ComputeTailProfile(const harp::vector_double &Params);
    double TailProfileValue(const double& dx, const double &dy) const;
    int psf_tail_amplitude_index;    
    
    /*
    double TailAmplitudeXW(const double &x, const double& wavelength, int fiber_bundle) const;
    double TailAmplitudeFW(const int fiber, const double& wavelength, int fiber_bundle) const;
    
    double TailValueWithParamsXY(const double &Xc, const double &Yc, 
				 const int IPix, const int JPix,
				 const harp::vector_double &Params,
				 harp::vector_double* derivative_r_tail_amplitude) const;
   
    double TailValueFW(const int fiber, const double& wavelength, 
		       const int IPix, const int JPix, int bundle_id, 
		       harp::vector_double* derivative_r_tail_amplitude) const;    
    */
       
#endif



public :



    const std::string& ParamName(int p) const;
    int ParamIndex(const std::string& name) const;
    bool HasParam(const std::string& name) const;
    
  protected :
    
    
    
    //harp::vector_double TmpParamDer;
    std::string name;

    std::map<int,Legendre1DPol*> XPol; // Legendre1DPol X_vs_W of fibers, data are in FiberTraces
    std::map<int,Legendre1DPol*> YPol; // Legendre1DPol X_vs_W of fibers, data are in FiberTraces
    
    
  public :

    TraceSet FiberTraces; // fiber traces, one independent trace per fiber

    const Trace& GetTrace(int fiber) const;
    Trace& GetTrace(int fiber);
    void AddTrace(int fiber);
    void LoadXYPol();
    double Xccd(int fiber, const double& wave) const;
    double Yccd(int fiber, const double& wave) const;
    
    
    long long int arc_exposure_id;
    long long int mjd;
    long long int plate_id;
    size_t ccd_image_n_cols;
    size_t ccd_image_n_rows;
    std::string camera_id;
    
    virtual std::string Name() const {return name;};
    
    
  protected :
    //! integrates PSF and requested derivatives over the pixel that contains XPix and YPix (pixel limits are at integer values + 1/2)
    //! called by public functions PSFValue... that can recast parameters
    
    virtual double PixValue(const double &Xc, const double &Yc,
		    const double &XPix, const double &YPix,
		    const harp::vector_double &Params,
		    harp::vector_double *PosDer = 0,
		    harp::vector_double *ParamDer = 0) const;

  public :

    virtual double Profile(const double &X, const double &Y,
			   const harp::vector_double &Params,
			   harp::vector_double *PosDer = 0,
			   harp::vector_double *ParamGradient = 0) const {return 0;};
    
    virtual harp::vector_double DefaultParams() const = 0;
    virtual std::vector<std::string> DefaultParamNames() const = 0;
    virtual void AllocateDefaultParams();
    
    // check if parameter values are within bounds
    virtual bool CheckParams(const harp::vector_double &Params) const {return false;};
    
    
    // ==================================================
    // END WAS BEFORE IN ANALYTIC PSF
    
    //! half size of the PSF in pixels in ccd
    int hSizeX, hSizeY; 
    
    double gain; // for Poisson noise, with gain in e/ADU
    double readout_noise; // rms value
    double psf_error; // fraction, like 0.01, accounts for flat and psf relative error
    
    std::map<int,PSF_Params> ParamsOfBundles;
    

    std::map<int,Prior*> Priors;
    void SetPrior(int k,Prior* p) {
      if(Priors.find(k) != Priors.end()) {
	delete Priors.find(k)->second;
      }
      Priors[k]=p;
    }
    void SetPrior(const string& pname,Prior* p) {
    int k=ParamIndex(pname);
    if(k<0) {
      cout << "WARNING no param " << pname << " in PSF" << endl;
      exit(12);
    }
    if(Priors.find(k) != Priors.end()) {
      delete Priors.find(k)->second;
    }
    cout << "DEBUG SetPrior '" <<  pname << "' index=" << k << endl;
    Priors[k]=p;
    }
    void ClearPriors() {
      for(std::map<int,Prior*>::const_iterator it = Priors.begin(); it != Priors.end() ; ++it) {
	delete it->second;
      }
      Priors.clear();
    }
    
    //#warning NEED TO IMPLEMENT JUMPS IN WAVELENGTH SOLUTION, AND NEED TO SAVE IN FITS
  
    
  
    
    virtual int LocalNAllPar() const = 0; // number of all the parameters parameters needed to describe psf at fixed ccd position
    int BundleNFitPar(int bundle_id) const; // number of parameters to be fitted  needed to describe psf varying with xy ccd coordinates
    int BundleNAllPar(int bundle_id) const; // number of all the parameters needed to describe psf varying with xy ccd coordinates
    int TracesNPar() const;
    
    bool IsLinear() const; // true if PSF linear wrt PSF params


    //! Constructor. RI should have been allocated via "new".
    PSF();
  
    //! limits of the stamp around Where StartI is in the stamp, EndI is beyond the last column.
    void StampLimits(const double &X, const double &Y,
		     int &BeginI, int &EndI,
		     int &BeginJ, int &EndJ) const;
    
    //! Access to the current PSF, with user provided Params.
    
    // this is the fastest (no conversion fiber,wave -> x,y)
    virtual double PSFValueWithParamsXY(const double& X, const double &Y, 
				const int IPix, const int JPix,
				const harp::vector_double &Params,
				harp::vector_double *PosDer, harp::vector_double *ParamDer,
				bool with_core=true, bool with_tail=true) const;
    
    double PSFValueWithParamsFW(const int fiber, const double &wave, 
				const int IPix, const int JPix,
				const harp::vector_double &Params,
				harp::vector_double *PosDer, harp::vector_double *ParamDer,
				bool with_core=true, bool with_tail=true) const;
    
    

    //! Access to the current PSF pixels.
    double PSFValueFW(const int fiber, const double &wave, 
		    const int IPix, const int JPix, int bundle_id,
		      harp::vector_double *PosDer = 0, harp::vector_double *ParamDer = 0,
		      bool with_core=true, bool with_tail=true) const;
    
    

    
    //! Access to current analytical PSF params (which may depend on position in the frame).
    int GetBundleOfFiber(int fiber) const;
    harp::vector_double AllLocalParamsFW(const int fiber, const double &wave, int bundle_id=-1) const;
    harp::vector_double AllLocalParamsXW(const double& x, const double &wave, int bundle_id) const;
    harp::vector_double AllLocalParamsFW_with_AllBundleParams(const int fiber, const double &wave, int bundle_id, const harp::vector_double& ForThesePSFParams) const;
    harp::vector_double AllLocalParamsXW_with_AllBundleParams(const double& x, const double &wave, int bundle_id, const harp::vector_double& ForThesePSFParams) const;
    harp::vector_double AllLocalParamsFW_with_FitBundleParams(const int fiber, const double &wave, int bundle_id, const harp::vector_double& ForThesePSFParams) const;
    harp::vector_double AllLocalParamsXW_with_FitBundleParams(const double& x, const double &wave, int bundle_id, const harp::vector_double& ForThesePSFParams) const;

    harp::vector_double FitLocalParamsFW(const int fiber, const double &wave, int bundle_id) const;
    harp::vector_double FitLocalParamsXW(const double& x, const double &wave, int bundle_id) const;
    harp::vector_double FitLocalParamsFW_with_FitBundleParams(const int fiber, const double &wave, int bundle_id, const harp::vector_double& ForThesePSFParams) const;
    harp::vector_double FitLocalParamsXW_with_FitBundleParams(const double& x, const double &wave, int bundle_id, const harp::vector_double& ForThesePSFParams) const;
    
    
    //! I/O for interface with specter
    /*
    virtual void WriteFits(fitsfile* fp, int first_hdu=1) const {};
    virtual void ReadFits(fitsfile* fp, int first_hdu=1) {};  
    virtual void WriteFits(const std::string& filename, int first_hdu=1) const;
    virtual void ReadFits(const std::string& filename, int first_hdu=1);
    */

    virtual void Append(const boost::shared_ptr < specex::PSF > other) = 0;
    

    ~PSF();
    
  private :

    template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
      ar & BOOST_SERIALIZATION_NVP(name);
      ar & BOOST_SERIALIZATION_NVP(hSizeX);
      ar & BOOST_SERIALIZATION_NVP(hSizeY);
      ar & BOOST_SERIALIZATION_NVP(ParamsOfBundles);
      ar & BOOST_SERIALIZATION_NVP(FiberTraces);
      ar & BOOST_SERIALIZATION_NVP(arc_exposure_id);
      ar & BOOST_SERIALIZATION_NVP(mjd);
      ar & BOOST_SERIALIZATION_NVP(plate_id);
      ar & BOOST_SERIALIZATION_NVP(camera_id);
      ar & BOOST_SERIALIZATION_NVP(ccd_image_n_cols);
      ar & BOOST_SERIALIZATION_NVP(ccd_image_n_rows);

      ar & BOOST_SERIALIZATION_NVP(gain);
      ar & BOOST_SERIALIZATION_NVP(readout_noise);
      ar & BOOST_SERIALIZATION_NVP(psf_error);

 
      return;
    }
    
    
    
  };
  
  BOOST_SERIALIZATION_ASSUME_ABSTRACT(PSF)
  
  BOOST_SERIALIZATION_SHARED_PTR(PSF)  
    
  typedef boost::shared_ptr < specex::PSF > PSF_p;
  typedef boost::weak_ptr < specex::PSF > PSF_wp;

};

#endif
