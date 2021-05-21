#ifndef SPECEX_PSF_FITTER__H
#define SPECEX_PSF_FITTER__H

#include <vector>
#include <string>
#include <map>

#include <specex_unbls.h>

#include "specex_psf.h"
#include "specex_spot.h"
#include "specex_stamp.h"
#include "specex_mask.h"
#include "specex_image_data.h"




namespace specex {

  class SpotTmpData  {
  
  public :

    
    double x,y,wavelength,fiber,fiber_bundle,flux;
    double frozen_x,frozen_y,frozen_flux; // for tails, keep fixed during minimization to avoid fitting flux with tails
    
    unbls::vector_double trace_x_monomials;
    unbls::vector_double trace_y_monomials;
    unbls::vector_double psf_monomials;
    unbls::vector_double psf_all_params;
    
    int flux_parameter_index;
    int x_parameter_index;
    int y_parameter_index;
    
    int trace_x_parameter_index;
    int trace_y_parameter_index;
    
#ifdef EXTERNAL_TAIL
    double tail_amplitude;
    unbls::vector_double tail_monomials;
#endif


    Stamp stamp;

    bool can_measure_flux;
    bool ignore;
  };

class PSF_Fitter {

 private :
  
  // internal to fitseveralspots
  unsigned npar_fixed_coord;
  unsigned npar_varying_coord;
  unsigned npar_trace;
  size_t nparTot;
  int index_of_spots_parameters;
  
  std::vector<SpotTmpData> spot_tmp_data;

#ifdef CONTINUUM
  size_t continuum_index;
#endif

 public :
  // internal set of parameters and matrices
  unbls::vector_double Params; // parameters that are fit (PSF, fluxes, XY CCD positions)
  std::vector<unbls::matrix_double> A_of_band; // for Gauss-Newton solving
  std::vector<unbls::vector_double> B_of_band; // for Gauss-Newton solving
  unbls::matrix_double fitWeight; // saved weight matrix of fitter parameters
  
 public :
  
  PSF_p psf;

 private :

  PSF_Params* psf_params; // this is a pointer set by the fitter to the current psf bundle being fit
  
 public :
  void SelectFiberBundle(int bundle); // this sets bundle_id and psf_global_params


  int number_of_image_bands; // for parallel processing (automatically set = to the variable OMP_NUM_THREADS of openmp)

  const image_data& image;
  const image_data& weight;
  const image_data& readnoise;
  image_data footprint_weight; // weight x psf footprint for global fit
  image_data corefootprint;  
  Stamp stamp; // rectangle in image where the fit occurs
  
  double corefootprint_weight_boost;
  bool fit_psf;
  bool fit_trace;
  bool fit_flux;
  bool fit_position;
  bool scheduled_fit_of_traces;
  bool scheduled_fit_of_sigmas;
  bool scheduled_fit_of_psf;
#ifdef EXTERNAL_TAIL
  bool fit_psf_tail;
  bool scheduled_fit_of_psf_tail;
#endif
#ifdef CONTINUUM
  bool fit_continuum;
  bool scheduled_fit_of_continuum;
#endif
  bool scheduled_fit_with_weight_model;
  bool sparse_pol;
  bool direct_simultaneous_fit;
  bool write_tmp_results;
  int trace_prior_deg;
  
  double chi2_precision;
  bool include_signal_in_weight;
  bool recompute_weight_in_fit;
  bool force_positive_flux;
  bool increase_weight_of_side_bands;
  bool fatal;
  bool parallelized;
  double polynomial_degree_along_x;
  double polynomial_degree_along_wave;
  
  int max_number_of_lines;
  
  Mask mask;
  
  std::map<std::string,Prior*> priors;
  
  std::map<int,int> tmp_trace_x_parameter;
  std::map<int,int> tmp_trace_y_parameter;
  
 PSF_Fitter(PSF_p i_psf, const image_data& i_image, const image_data& i_weight, const image_data& i_readnoise) :
    
  psf(i_psf),
    image(i_image),
    weight(i_weight),
    readnoise(i_readnoise),
    stamp(i_image),
    corefootprint_weight_boost(0), 
    fit_psf(false),
    fit_trace(false),
    fit_flux(false),
    fit_position(false), 
    chi2_precision(0.1),
    include_signal_in_weight(false),
    recompute_weight_in_fit(false),
   force_positive_flux(false),
   increase_weight_of_side_bands(false),
    scheduled_fit_of_traces(true),
    scheduled_fit_of_sigmas(true),
    scheduled_fit_of_psf(true),
#ifdef EXTERNAL_TAIL
    fit_psf_tail(false),
    scheduled_fit_of_psf_tail(false),
#endif
#ifdef CONTINUUM
    fit_continuum(false),
    scheduled_fit_of_continuum(false),
#endif
    scheduled_fit_with_weight_model(false),
    sparse_pol(true),
    direct_simultaneous_fit(false),
    write_tmp_results(false),
    trace_prior_deg(0),
    fatal(true),
    parallelized(true),
        
    polynomial_degree_along_x(1),
    polynomial_degree_along_wave(4),
    max_number_of_lines(0)


      {
      };
    
    void SetStampLimitsFromPSF(Stamp& stamp, const PSF_p psf, const double &X, const double &Y);
    void SetStampLimitsFromPSF(Stamp& stamp, const PSF_p psf, const double &xc_min, const double &xc_max, const double &yc_min, const double &yc_max);

    int NPar(int nspots) const;
   //int Index_PSF() const;
   //int Index_Flux(int spotid, int nspots) const;
    
    void InitTmpData(const std::vector<Spot_p>& spots);
    void UpdateTmpData(bool compute_ab);
    double ParallelizedComputeChi2AB(bool compute_ab);
    double ComputeChi2AB(bool compute_ab, int begin_j=0, int end_j=0, unbls::matrix_double* Ap=0, unbls::vector_double* Bp=0, bool update_tmp_data=true) const;

  void ComputeWeigthImage(std::vector<specex::Spot_p>& spots, int* npix);

  bool FitOneSpot(Spot_p& spot, double *chi2_val=0, int *n_iterations=0);
  bool FitSeveralSpots(std::vector<Spot_p>& spots, double *chi2_val=0, int *n_pixels=0, int *n_iterations=0);
  
  bool FitIndividualSpotFluxes(std::vector<Spot_p>& spots);
  bool FitIndividualSpotPositions(std::vector<Spot_p>& spots);
  bool FitEverything(std::vector<Spot_p>& spots, bool init_psf=false);
  
  void compare_spots_chi2_and_mask(std::vector<specex::Spot_p>& spots, const double& nsig=4.);
  std::vector<specex::Spot_p> select_spots(std::vector<specex::Spot_p>& input_spots, double minimum_signal_to_noise, double min_wave_dist=0, double chi2_nsig=4);
  

};


}

#endif
