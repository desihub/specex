#ifndef SPECEX_PSF_FITTER__H
#define SPECEX_PSF_FITTER__H

#include <vector>
#include <string>
#include <map>

#include "harp.hpp"

#include "specex_psf.h"
#include "specex_spot.h"
#include "specex_stamp.h"
#include "specex_mask.h"
#include "specex_image_data.h"

namespace specex {

class PSF_Fitter {

 private :
  
  // internal to fitseveralspots
  int nspots;
  unsigned npar_fixed_coord;
  unsigned npar_varying_coord;
  std::map<int,int> XvsW_index_of_fiber;
  std::map<int,int> YvsW_index_of_fiber;
  harp::vector_double gradPar,gradPos;
  harp::vector_double spots_flux,spots_x,spots_y;
  size_t nparTot;
  int npar_traces;
  vector<harp::vector_double> TraceXvsW_Monomials_of_spots;
  vector<harp::vector_double> TraceYvsW_Monomials_of_spots;
  vector<Stamp> spot_stamps;
  
 public :
  // internal set of parameters and matrices
  harp::vector_double Params; // parameters that are fit (PSF, fluxes, XY CCD positions)
  harp::matrix_double A; // for Gauss-Newton solving
  harp::vector_double B; // for Gauss-Newton solving
  harp::matrix_double fitWeight; // saved weight matrix of fitter parameters
  
 public :
  
  PSF_p psf;
  const image_data& image;
  const image_data& weight;
  image_data footprint_weight; // weight x psf footprint for global fit
  Stamp stamp; // rectangle in image where the fit occurs
  bool fit_psf;
  bool fit_trace;
  bool fit_flux;
  bool fit_position;
  double chi2_precision;
  bool include_signal_in_weight;
  bool verbose;
  bool fatal;
  double gain; // for Poisson noise, with gain in e/ADU
  double readout_noise; // rms value
  double flatfield_error; // fraction, like 0.01, accounts for flat and psf relative error
  // model parameters 
  int polynomial_degree_along_x;
  int polynomial_degree_along_y;
  
  
  Mask mask;
  


  PSF_Fitter(PSF_p i_psf, const image_data& i_image, const image_data& i_weight) :
    
    psf(i_psf),
    image(i_image),
    weight(i_weight),
    stamp(i_image), 
    fit_psf(true),
    fit_trace(false),
    fit_flux(true),
    fit_position(true), 
    chi2_precision(0.1),
    include_signal_in_weight(false),
    fatal(true),
    verbose(true),
    gain(1),
    readout_noise(0),
    flatfield_error(0),
    polynomial_degree_along_x(4),
    polynomial_degree_along_y(4)
      {
      };
    
    void SetStampLimitsFromPSF(Stamp& stamp, const PSF_p psf, const double &X, const double &Y);
    void SetStampLimitsFromPSF(Stamp& stamp, const PSF_p psf, const double &xc_min, const double &xc_max, const double &yc_min, const double &yc_max);

    int NPar(int nspots) const;
   //int Index_PSF() const;
   //int Index_Flux(int spotid, int nspots) const;
   

   double ComputeChi2AB(std::vector<Spot*>& spots, bool compute_ab) ;

  

   void SetPSFParams(const harp::vector_double &Params); 

  bool FitOneSpot(Spot& spot, double *chi2_val=0, int *n_iterations=0);
  bool FitSeveralSpots(std::vector<Spot*>& spots, double *chi2_val=0, int *n_pixels=0, int *n_iterations=0);
  
  bool InterpolateSpotPSFs(std::vector<Spot*>& spots, double *chi2_val=0, int *n_iterations=0);
  bool FitTraces(std::vector<Spot*>& spots, int *nok=0);

  
  bool FitIndividualSpotFluxes(std::vector<Spot*>& spots);
  bool FitIndividualSpotPositions(std::vector<Spot*>& spots);
  bool FitEverything(std::vector<Spot*>& spots, bool init_psf=false);
  
};


}

#endif
