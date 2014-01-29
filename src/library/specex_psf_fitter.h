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

  class SpotTmpData  {
  
  public :

    
    double x,y,wavelength,flux;
    double frozen_x,frozen_y,frozen_flux; // for tails, keep fixed during minimization to avoid fitting flux with tails
    
    harp::vector_double trace_x_monomials;
    harp::vector_double trace_y_monomials;
    harp::vector_double psf_monomials;
    harp::vector_double psf_all_params;
    
    int flux_parameter_index;
    int x_parameter_index;
    int y_parameter_index;
    
    int trace_x_parameter_index;
    int trace_y_parameter_index;
    
    Stamp stamp;
    
  };

class PSF_Fitter {

 private :
  
  // internal to fitseveralspots
  unsigned npar_fixed_coord;
  unsigned npar_varying_coord;
  size_t nparTot;
  
  std::vector<SpotTmpData> spot_tmp_data;
#ifdef EXTERNAL_TAIL
  size_t psf_r_tail_amplitude_index;
#ifdef EXTERNAL_Y_TAIL
  size_t psf_y_tail_amplitude_index;
#endif
#endif
#ifdef CONTINUUM
  size_t continuum_index;
#endif

 public :
  // internal set of parameters and matrices
  harp::vector_double Params; // parameters that are fit (PSF, fluxes, XY CCD positions)
  std::vector<harp::matrix_double> A_of_chunk; // for Gauss-Newton solving
  std::vector<harp::vector_double> B_of_chunk; // for Gauss-Newton solving
  std::vector<double> chi2_of_chunk; // for Gauss-Newton solving
  harp::matrix_double fitWeight; // saved weight matrix of fitter parameters
  
 public :
  
  PSF_p psf;

 private :

  PSF_Params* psf_params; // this is a pointer set by the fitter to the current psf bundle being fit
  
 public :
  void SelectFiberBundle(int bundle); // this sets bundle_id and psf_global_params


  int number_of_image_chuncks; // for parallel processing (automatically set = to the variable OMP_NUM_THREADS of openmp)

  const image_data& image;
  const image_data& weight;
  image_data footprint_weight; // weight x psf footprint for global fit
  Stamp stamp; // rectangle in image where the fit occurs
  bool fit_psf;
  bool fit_trace;
  bool fit_flux;
  bool fit_position;
  bool scheduled_fit_of_traces;
#ifdef EXTERNAL_TAIL
  bool fit_psf_tail;
  bool scheduled_fit_of_psf_tail;
#endif
#ifdef CONTINUUM
  bool fit_continuum;
  bool scheduled_fit_of_continuum;
#endif


  double chi2_precision;
  bool include_signal_in_weight;
  bool verbose;
  bool fatal;
  double gain; // for Poisson noise, with gain in e/ADU
  double readout_noise; // rms value
  double psf_error; // fraction, like 0.01, accounts for flat and psf relative error
  
  double polynomial_degree_along_x;
  double polynomial_degree_along_wave;
  
  Mask mask;
  


 PSF_Fitter(PSF_p i_psf, const image_data& i_image, const image_data& i_weight) :
    
  psf(i_psf),
    image(i_image),
    weight(i_weight),
    stamp(i_image), 
    fit_psf(false),
    fit_trace(false),
    fit_flux(false),
    fit_position(false), 
    chi2_precision(0.1),
    include_signal_in_weight(false),
    scheduled_fit_of_traces(true),
#ifdef EXTERNAL_TAIL
    fit_psf_tail(false),
    scheduled_fit_of_psf_tail(false),
#endif
#ifdef CONTINUUM
    fit_continuum(false),
    scheduled_fit_of_continuum(false),
#endif
    fatal(true),
    verbose(true),
    gain(1),
    readout_noise(0),
      psf_error(0.01),
      polynomial_degree_along_x(1),
      polynomial_degree_along_wave(4)
      {
	


      };
    
    void SetStampLimitsFromPSF(Stamp& stamp, const PSF_p psf, const double &X, const double &Y);
    void SetStampLimitsFromPSF(Stamp& stamp, const PSF_p psf, const double &xc_min, const double &xc_max, const double &yc_min, const double &yc_max);

    int NPar(int nspots) const;
   //int Index_PSF() const;
   //int Index_Flux(int spotid, int nspots) const;
    
    void UpdateTmpData(bool compute_ab);
    double ParallelizedComputeChi2AB(bool compute_ab);
    double ComputeChi2AB(bool compute_ab, int begin_j=0, int end_j=0, harp::matrix_double* Ap=0, harp::vector_double* Bp=0, bool update_tmp_data=true) const;

   void ComputeWeigthImage(std::vector<specex::Spot_p>& spots, int* npix);

   // void SetAllPSFParams(const harp::vector_double &Params); 

  bool FitOneSpot(Spot_p& spot, double *chi2_val=0, int *n_iterations=0);
  bool FitSeveralSpots(std::vector<Spot_p>& spots, double *chi2_val=0, int *n_pixels=0, int *n_iterations=0);
  
  //bool InterpolateSpotPSFs(std::vector<Spot_p>& spots, double *chi2_val=0, int *n_iterations=0);
  bool FitTraces(std::vector<Spot_p>& spots, int *nok=0);

  
  bool FitIndividualSpotFluxes(std::vector<Spot_p>& spots);
  bool FitIndividualSpotPositions(std::vector<Spot_p>& spots);
  bool FitEverything(std::vector<Spot_p>& spots, bool init_psf=false);
  
};


}

#endif
