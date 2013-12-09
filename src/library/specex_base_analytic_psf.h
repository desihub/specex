#ifndef SPECEX_BASE_ANALYTIC_PSF__H
#define SPECEX_BASE_ANALYTIC_PSF__H


#include <string>
#include <vector>



#define PSF_NAN_VALUE 9999999


namespace specex {

  //! The virtual (interface) class. It holds the PixValue routine that integrates over pixels.
  class AnalyticPSF {
  
  protected :
    //  ProfileFunc *profile;
    // const int npar;
    std::vector<std::string> paramNames;
    
    
    harp::vector_double TmpParamDer;
    
  public :
    
    //     AnalyticPSF(const string &Name, ProfileFunc *Profile, 
    //      const int NPar , const vector<string> ParamNames);
    
    virtual std::string Name() const = 0;
    virtual std::string ParamName(const int Rank) const {return paramNames.at(unsigned(Rank));};
    virtual int ParamIndex(const std::string& name) const {
      for(size_t i=0;i<paramNames.size();++i) {
	if(paramNames[i] == name) return int(i);
      }
      return -1;
    };
    virtual bool HasParam(const std::string& name) const { return (ParamIndex(name)>=0);}
    
    virtual size_t NPar() const  = 0;
    
    //! integrates PSF and requested derivatives over the pixel that contains XPix and YPix (pixel limits are at integer values + 1/2)
    double PixValue(const double &Xc, const double &Yc,
		    const double &XPix, const double &YPix,
		    const harp::vector_double &Params,
		    harp::vector_double *PosDer = 0,
		    harp::vector_double *ParamDer = 0) const;
    
    virtual double Profile(const double &X, const double &Y,
			   const harp::vector_double &Params,
			   harp::vector_double *PosDer = 0,
			   harp::vector_double *ParamGradient = 0) const = 0;
    
    virtual void InitParams(const double &sigma, harp::vector_double &Params) = 0;
    
    // check if parameter values are within bounds
    virtual bool CheckParams(const harp::vector_double &Params) const = 0;
    
    virtual ~AnalyticPSF() {};
  };
}



#endif /* ANALYTICPSF__H */
