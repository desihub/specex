#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <map>

#include <boost/program_options.hpp>

#include <harp.hpp>

#include <specex_message.h>

#include <specex_serialization.h>

#include <specex_trace.h>
#include <specex_spot.h>
#include <specex_spot_array.h>
#include <specex_lamp_lines_utils.h>
#include <specex_fits.h>
#include <specex_desi_io.h>
#include <specex_psf_fitter.h>

#include <specex_psf_io.h>

#include <specex_tokens.h>

#include <specex_pyoptions.h>
#include <specex_pyio.h>
#include <specex_pyprior.h>
#include <specex_pyimage.h>
#include <specex_pyfitting.h>

using namespace std;
using namespace specex;

#define _GNU_SOURCE 1
#ifndef __USE_GNU
#define __USE_GNU
#endif
#include <portable_fenv.h>

/*
  input format of fits file
  HDU 0 : blank
  HDU 1 ELECTRONS : 2D image of noisy electrons
  - RDNOISE is in header
  HDU 2 IVAR : Inverse variance [1/electrons^2]
  - ivar = 1/(pix.clip(0) + rdnoise**2)
  HDU 3 MASK : 0=good.  Currently all 0.
  HDU 4 XCOEFF : Legendre coefficients for mapping wavelength -> x
  - WAVEMIN, WAVEMAX : domain for mapping to [-1,1] for Legendre polynomials
  - image is coefficients for each fiber
  HDU 5 YCOEFF : Legendre coefficients for mapping wavelength -> y
  HDU 5 TRUE_ELECTRONS : original noiseless image in electrons
  ELECTRONS = poisson(TRUE_ELECTRONS) + gaussian(rdnoise)
 */

int specex::PyFitting::fit_psf(
			      specex::PyOptions opts,
			      specex::PyIO      pyio,
			      specex::PyPrior   pypr,
			      specex::PyImage   pymg,
			      specex::PyPSF&    pyps
			      ){
  
  // to crash when NaN

  feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);

  // copy to local psf pointer 
  specex::PSF_p psf = pyps.psf;

  // set logging info

  specex_set_debug(opts.vm.count("debug")>0);
  specex_set_verbose(opts.vm.count("quiet")==0);
  specex_set_dump_core(opts.vm.count("core")>0);

  bool fit_traces = (opts.vm.count("no-trace-fit")==0);
  bool fit_sigmas = (opts.vm.count("no-sigma-fit")==0);
  bool fit_thepsf = (opts.vm.count("no-psf-fit")==0);
  if (opts.vm.count("only-trace-fit")>0) {
    if (opts.vm.count("no-trace-fit")>0) {
      SPECEX_ERROR("cannot have both options --only-trace-fit and --no-trace-fit");
      return EXIT_FAILURE;
    }
    fit_traces = true;
    fit_sigmas = false;
    fit_thepsf = false;      
  }
  
  bool single_bundle = (opts.vm.count("single-bundle")>0);
  
  bool write_tmp_results = (opts.vm.count("tmp-results")>0);
  
  if(opts.trace_prior_deg>0) SPECEX_INFO("Will apply a prior on the traces high order terms in a bundle");
    
#ifdef EXTERNAL_TAIL
  bool fit_psf_tails = (opts.vm.count("fit-psf-tails")>0);
#endif
#ifdef CONTINUUM
  bool fit_continuum = (opts.vm.count("fit-continuum")>0);
#endif
  bool use_variance_model = (opts.vm.count("variance-model")>0);
  bool fit_individual_spots_position = opts.vm.count("positions");
  
  SPECEX_INFO("using lamp lines file " << opts.lamp_lines_filename); 
  
  // set broken fibers
  
  vector<string> tokens = split(opts.broken_fibers_string,',');
  vector<int> broken_fibers;
  for(size_t i=0;i<tokens.size();i++) {
    broken_fibers.push_back(atoi(tokens[i].c_str())%500);
  }
  SPECEX_DEBUG("number of input broken fibers: " << broken_fibers.size());
  for(size_t i=0;i<broken_fibers.size();i++) {
    SPECEX_DEBUG("input broken fiber #" << broken_fibers[i]);
  }
  
  for(size_t i=0;i<broken_fibers.size();i++) {
    psf->FiberTraces[broken_fibers[i]].mask = 3;
    SPECEX_DEBUG("Fiber trace " << broken_fibers[i] << " OFF=" << psf->FiberTraces[broken_fibers[i]].Off());
  }
        
  psf->ccd_image_n_cols = pymg.image.n_cols();
  psf->ccd_image_n_rows = pymg.image.n_rows();
  
  // bundle sizes
  int number_of_fibers_per_bundle=0;
  if(single_bundle) {
    number_of_fibers_per_bundle = psf->FiberTraces.size();
  }else {
    number_of_fibers_per_bundle = eval_bundle_size(psf->FiberTraces);
  }
  int number_of_fiber_bundles_per_ccd=psf->FiberTraces.size()/number_of_fibers_per_bundle;          
  if(opts.first_fiber_bundle<0 || opts.first_fiber_bundle>= number_of_fiber_bundles_per_ccd) {
    SPECEX_ERROR("invalid first fiber bundle");
  }
  if(opts.last_fiber_bundle<opts.first_fiber_bundle || opts.last_fiber_bundle>= number_of_fiber_bundles_per_ccd) {
    SPECEX_ERROR("invalid last fiber bundle");
  }   
    
  // main body
  
  try {
    // init PSF fitter
    // -------------------------------------------- 
    PSF_Fitter fitter(psf,pymg.image,pymg.weight,pymg.rdnoise);
    
    fitter.polynomial_degree_along_x    = opts.legendre_deg_x;
    fitter.polynomial_degree_along_wave = opts.legendre_deg_wave;
    fitter.psf->psf_error               = opts.psf_error;
    fitter.corefootprint_weight_boost   = opts.psf_core_wscale;
    fitter.write_tmp_results            = write_tmp_results;
    fitter.trace_prior_deg              = opts.trace_prior_deg;
    
#ifdef EXTERNAL_TAIL
    fitter.scheduled_fit_of_psf_tail    = fit_psf_tails;
#endif

#ifdef CONTINUUM
    fitter.scheduled_fit_of_continuum   = fit_continuum;
#endif

    fitter.scheduled_fit_with_weight_model  = use_variance_model;
    
    fitter.scheduled_fit_of_traces      = fit_traces;
    fitter.scheduled_fit_of_sigmas      = fit_sigmas;
    fitter.scheduled_fit_of_psf         = fit_thepsf;
    fitter.direct_simultaneous_fit      = true; // use_input_specex_psf;
    fitter.max_number_of_lines          = opts.max_number_of_lines;
    
    fitter.psf->gain = 1; // images are already in electrons
    fitter.psf->readout_noise = 0; // readnoise is a property of image, not PSF
        
    if(pymg.header.find("CAMERA") != pymg.header.end()) {
      fitter.psf->camera_id = pymg.header["CAMERA"];
      SPECEX_INFO("CAMERA = " << fitter.psf->camera_id );
    }else{
      SPECEX_WARNING("CAMERA Id not found in header");
    }
    
    fitter.priors = pypr.priors;
        
    fitter.mask.Clear();

    SPECEX_INFO(
		"PSF '" << opts.psf_model << "' stamp size = " 
		<< psf->hSizeX << "x" << psf->hSizeY );    
    
    int first_fitted_fiber=-1;
    int last_fitted_fiber=-1;
    
    // loop on fiber bundles 
    // -------------------------------------------- 
    for(int bundle = opts.first_fiber_bundle; bundle <= opts.last_fiber_bundle ; bundle ++) {
    
      // allocate bundle in PSF if necessary
      if(psf->ParamsOfBundles.find(bundle)==psf->ParamsOfBundles.end()) {
	psf->ParamsOfBundles[bundle] = specex::PSF_Params();
	psf->ParamsOfBundles[bundle].bundle_id = bundle;
	psf->ParamsOfBundles[bundle].fiber_min = number_of_fibers_per_bundle*bundle;
	psf->ParamsOfBundles[bundle].fiber_max = psf->ParamsOfBundles[bundle].fiber_min+number_of_fibers_per_bundle-1; // included
      }
      
      // now check mask ?
      if(psf->ParamsOfBundles[bundle].fiber_min < opts.first_fiber) {
	psf->ParamsOfBundles[bundle].fiber_min = opts.first_fiber;
	SPECEX_INFO("restricting fiber range first fiber = " << opts.first_fiber);
      }
      if(psf->ParamsOfBundles[bundle].fiber_max > opts.last_fiber) {
	psf->ParamsOfBundles[bundle].fiber_max = opts.last_fiber;
	SPECEX_INFO("restricting fiber range last fiber = " << opts.last_fiber);
      }
      
      fitter.SelectFiberBundle(bundle);
      
      // loading arc lamp spots belonging to this bundle
      // --------------------------------------------
      
      int ymin = 0; // range of usable CCD coordinates, hard coded for now
      int ymax = pymg.image.n_rows(); // range of usable CCD coordinates, hard coded for now
      
      /*
      SPECEX_WARNING("RESTRICTING Y RANGE !!!!!");
      if(psf->camera_id=="b1") {ymin=696; ymax = 3516;};
      if(psf->camera_id=="b2") {ymin=696; ymax = 3516;};
      if(psf->camera_id=="r1") {ymin=200; ymax = 3668;}; 
      if(psf->camera_id=="r2") {ymin=200; ymax = 3668;};
      */

      int margin = -psf->hSizeY+1; // we need to include spots that contribute to the image signal
      ymin+=margin;
      ymax-=margin;

      SPECEX_INFO("valid y(=rows) range = " << ymin << " " << ymax);

      vector<Spot_p> spots;
      
      double min_wavelength = 0;
      double max_wavelength = 1e6;  
      allocate_spots_of_bundle(spots,opts.lamp_lines_filename,psf->FiberTraces,
			       bundle,psf->ParamsOfBundles[bundle].fiber_min,
			       psf->ParamsOfBundles[bundle].fiber_max,ymin,ymax,
			       min_wavelength,max_wavelength);
      SPECEX_INFO("number of spots = " << spots.size());
      
      //if(write_tmp_results)
      //  write_spots_xml(spots,"spots-init.xml");
      
      // starting fit
      // --------------------------------------------
      bool init_psf = (!pyio.use_input_specex_psf);
      fitter.FitEverything(spots,init_psf);

      int ndf = psf->ParamsOfBundles[bundle].ndata - psf->ParamsOfBundles[bundle].nparams;
      SPECEX_INFO("Bundle " << bundle << " PSF fit status   = " <<
		  psf->ParamsOfBundles[bundle].fit_status);
      SPECEX_INFO("Bundle " << bundle << " PSF fit chi2/ndf = " <<
		  psf->ParamsOfBundles[bundle].chi2 << "/" << ndf <<
		  " = " << psf->ParamsOfBundles[bundle].chi2/ndf);
      SPECEX_INFO("Bundle " << bundle << " PSF fit ndata    = "<<
		  psf->ParamsOfBundles[bundle].ndata);
      SPECEX_INFO("Bundle " << bundle << " PSF fit nspots   = "<<
		  psf->ParamsOfBundles[bundle].nspots_in_fit);
      SPECEX_INFO("Bundle " << bundle << " PSF fit chi2/ndata (core) = "<<
		  psf->ParamsOfBundles[bundle].chi2_in_core << "/" <<
		  psf->ParamsOfBundles[bundle].ndata_in_core << " = " <<
		  psf->ParamsOfBundles[bundle].chi2_in_core /
		  psf->ParamsOfBundles[bundle].ndata_in_core);
      
      if(bundle == opts.first_fiber_bundle) {
	first_fitted_fiber=psf->ParamsOfBundles[bundle].fiber_min;
	last_fitted_fiber=psf->ParamsOfBundles[bundle].fiber_max;
      }
      
      first_fitted_fiber=min(first_fitted_fiber,psf->ParamsOfBundles[bundle].fiber_min);
      last_fitted_fiber=max(last_fitted_fiber,psf->ParamsOfBundles[bundle].fiber_max);
          
      for(size_t s=0;s<spots.size();s++) {
	pyps.fitted_spots.push_back(spots[s]);
      }

      if(fit_individual_spots_position) // for debugging
	fitter.FitIndividualSpotPositions(spots);

    } // end of loop on bundles

    pyps.psf = fitter.psf;

  // ending
  // --------------------------------------------   
   
  } catch(harp::exception e) {
    cerr << "FATAL ERROR (harp) " << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::exception& e) {
    cerr << "FATAL ERROR (other std) " << e.what() <<endl;
    return EXIT_FAILURE;

  }catch (...) {
    cerr << "FATAL ERROR (unknown)" << endl;
    return EXIT_FAILURE;
 }
  
  // may prevent crashing on non-floating point exceptions outside this function
  fedisableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
  return EXIT_SUCCESS;
}


#include <specex_serialization_implement.h>


