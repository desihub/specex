#ifndef SPECEX_PSF__H
#define SPECEX_PSF__H

#include <vector>
#include <string>
#include <map>

#include "specex_legendre.h"
#include "specex_trace.h"
#include "specex_linalg.h"

#include "specex_base_analytic_psf.h"

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
  
  
  class PSF {
    
  protected :
    
    
  public :
    
    //! analytic function of x and y
    AnalyticPSF* analyticPSF;
    
    //! half size of the PSF in pixels in ccd
    int hSizeX, hSizeY; 
    
    //! parameters of psf shape varying with xy coordinates
    std::vector<Legendre2DPol> Params;
    
    std::map<int,Prior*> Priors;
    void SetPrior(int k,Prior* p) {
      if(Priors.find(k) != Priors.end()) {
	delete Priors.find(k)->second;
      }
      Priors[k]=p;
    }
    void SetPrior(const string& pname,Prior* p) {
    int k=analyticPSF->ParamIndex(pname);
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
  
  
    //! fiber traces, one independent trace per fiber
    std::map<int,Trace> FiberTraces;
  
    int FixedCoordNPar() const; // set of parameters needed to describe psf at fixed ccd position
    int VaryingCoordNPar() const; // set of parameters needed to describe psf varying with xy ccd coordinates
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
			      harp::vector_double *PosDer, harp::vector_double *ParamDer, double *AnalyticValue = 0) const;

    //! Access to the current PSF pixels.
    double PSFValue(const double &Xc, const double &Yc, 
		    const int IPix, const int JPix,
		    harp::vector_double *PosDer = 0, harp::vector_double *ParamDer = 0, double *AnalyticValue = NULL) const;
    
    //! Access to current analytical PSF params (which may depend on position in the frame).
    harp::vector_double FixedCoordParams(const double &X, const double &Y) const;
    harp::vector_double FixedCoordParams(const double &X, const double &Y, const harp::vector_double& ForThesePSFParams) const;
    
    ~PSF();

    bool verbose;

  /*
    void write(std::ostream &os) const;
    bool read(std::istream &is);
    void write(const std::string &FileName) const;
    bool read(const std::string &FileName) ;
  */

  };

};




#endif
