#ifndef SPECEX_PSF__H
#define SPECEX_PSF__H

#include <vector>
#include <string>
#include <map>

#include "specex_legendre.h"
#include "specex_trace.h"
#include "specex_linalg.h"

#define PSF_NAN_VALUE 9999999

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
    
    std::vector<Legendre2DPol> Polynomials;
    int bundle_id;
    int fiber_min; // first fiber this set of params applies to
    int fiber_max; // last fiber this set of params applies to
    PSF_Params() : bundle_id(0), fiber_min(0), fiber_max(0) {};
  
  private :

    template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
      ar & BOOST_SERIALIZATION_NVP(Polynomials);
      ar & BOOST_SERIALIZATION_NVP(bundle_id);
      ar & BOOST_SERIALIZATION_NVP(fiber_min);
      ar & BOOST_SERIALIZATION_NVP(fiber_max);
    }
  };

  class PSF {

    friend class boost::serialization::access;

    // AnalyticPSF* analyticPSF;
    // BEGIN WAS BEFORE IN ANALYTIC PSF
    // ==================================================
 
  protected :
    
    std::vector<std::string> paramNames;
    harp::vector_double TmpParamDer;
    std::string name;

  public :

    long long int arc_exposure_id;
    long long int mjd;
    long long int plate_id;
    size_t ccd_image_n_cols;
    size_t ccd_image_n_rows;
    std::string camera_id;
    
    virtual std::string Name() const {return name;};
    virtual std::string ParamName(const int Rank) const {return paramNames.at(unsigned(Rank));};
    virtual int ParamIndex(const std::string& name) const {
      for(size_t i=0;i<paramNames.size();++i) {
	if(paramNames[i] == name) return int(i);
      }
      return -1;
    };
    virtual bool HasParam(const std::string& name) const { return (ParamIndex(name)>=0);}
    
    virtual size_t NPar() const {return 0;};
    
    //! integrates PSF and requested derivatives over the pixel that contains XPix and YPix (pixel limits are at integer values + 1/2)
    double PixValue(const double &Xc, const double &Yc,
		    const double &XPix, const double &YPix,
		    const harp::vector_double &Params,
		    harp::vector_double *PosDer = 0,
		    harp::vector_double *ParamDer = 0) const;
    
    virtual double Profile(const double &X, const double &Y,
			   const harp::vector_double &Params,
			   harp::vector_double *PosDer = 0,
			   harp::vector_double *ParamGradient = 0) const {return 0;};
    
    virtual void InitParams(const double &sigma, harp::vector_double &Params) {};
    
    // check if parameter values are within bounds
    virtual bool CheckParams(const harp::vector_double &Params) const {return false;};
    
    
    // ==================================================
    // END WAS BEFORE IN ANALYTIC PSF
    
    //! half size of the PSF in pixels in ccd
    int hSizeX, hSizeY; 
    
    //! parameters of psf shape varying with xy coordinates
    //std::vector<Legendre2DPol> Params;
    
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
  
    //! fiber traces, one independent trace per fiber
    std::map<int,Trace> FiberTraces;
    
    int FixedCoordNPar() const; // set of parameters needed to describe psf at fixed ccd position
    int VaryingCoordNPar(int bundle_id) const; // set of parameters needed to describe psf varying with xy ccd coordinates
    int TracesNPar() const;
  
    bool IsLinear() const; // true if PSF linear wrt PSF params


    //! Constructor. RI should have been allocated via "new".
    PSF();
  
    //! limits of the stamp around Where StartI is in the stamp, EndI is beyond the last column.
    void StampLimits(const double &X, const double &Y,
		     int &BeginI, int &EndI,
		     int &BeginJ, int &EndJ) const;
    
    //! Access to the current PSF, with user provided Params.
    double PSFValueWithParams(const double &Xc, const double &Yc, 
			      const int IPix, const int JPix,
			      const harp::vector_double &Params,
			      harp::vector_double *PosDer, harp::vector_double *ParamDer) const;

    //! Access to the current PSF pixels.
    double PSFValue(const double &Xc, const double &Yc, 
		    const int IPix, const int JPix, int bundle_id,
		    harp::vector_double *PosDer = 0, harp::vector_double *ParamDer = 0) const;
    
    //! Access to current analytical PSF params (which may depend on position in the frame).
    harp::vector_double FixedCoordParams(const double &X, const double &Y, int bundle_id) const;
    harp::vector_double FixedCoordParams(const double &X, const double &Y, int bundle_id, const harp::vector_double& ForThesePSFParams) const;
    
    //! I/O for interface with specter
    virtual void WriteFits(fitsfile* fp, int first_hdu=1) const {};
    virtual void ReadFits(fitsfile* fp, int first_hdu=1) {};  
    virtual void WriteFits(const std::string& filename, int first_hdu=1) const;
    
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
      
      return;
    }
    
    
    
  };
  
  BOOST_SERIALIZATION_ASSUME_ABSTRACT(PSF)
  
  BOOST_SERIALIZATION_SHARED_PTR(PSF)  
    
  typedef boost::shared_ptr < specex::PSF > PSF_p;
  typedef boost::weak_ptr < specex::PSF > PSF_wp;

};

#endif
