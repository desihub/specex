#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <map>

#include <boost/program_options.hpp>

#include <boost/archive/xml_oarchive.hpp>

#include <harp.hpp>

#include <specex_message.h>
#include <specex_psf.h>
#include <specex_trace.h>
#include <specex_spot.h>
#include <specex_spot_array.h>
#include <specex_spectrograph.h>
#include <specex_lamp_lines_utils.h>
#include <specex_fits.h>
#include <specex_desi_io.h>
#include <specex_psf_fitter.h>

#include <specex_psf_io.h>
#include <specex_serialization.h>

using namespace std;
using namespace specex;

namespace popts = boost::program_options;

#define _GNU_SOURCE 1
#ifndef __USE_GNU
#define __USE_GNU
#endif
#include <fenv.h>

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



int main ( int argc, char *argv[] ) {
  
  // to crash when NaN
  feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
  

  // default arguments
  // --------------------------------------------
  string psf_model = "GAUSSHERMITE";
  string spectrograph_name = "DESI";  
  int    first_fiber_bundle=1;
  int    last_fiber_bundle=1;
  int    first_fiber=0;
  int    last_fiber=100000;
  int    half_size_x=4;
  int    half_size_y=4;
  
  string arc_image_filename="";
  //string xy_trace_fits_name="";
  //string wy_trace_fits_name="";
  string lamp_lines_filename="";  
  double min_wavelength = 0;
  double max_wavelength = 1e6;
  int gauss_hermite_deg  = 3;
  int gauss_hermite_deg2 = 2;
  //double gauss_hermite_sigma = 1.1;
  int legendre_deg_wave = 4;
  int legendre_deg_x = 1;
  
  double psf_error = 0;
  double psf_core_wscale = 0;
  
  string output_xml_filename="";
  string output_fits_filename="";
  string output_spots_filename="";
  
  bool write_tmp_results = false;
  int flux_hdu=1;
  int ivar_hdu=2;
  
  if(getenv("SPECEXDATA"))
    lamp_lines_filename = string(getenv("SPECEXDATA"))+"/lamplines-specex.par";
   
  // reading arguments
  // --------------------------------------------
  popts::options_description desc ( "Allowed Options" );
  desc.add_options()
    ( "help,h", "display usage information" )
    ( "arc,a", popts::value<string>( &arc_image_filename ), "arc pre-reduced fits image file name (mandatory), ex:  sdProc-b1-00108382.fits" )
    ( "flux-hdu", popts::value<int>( &flux_hdu ), " flux hdu in input arc fits")
    ( "ivar-hdu", popts::value<int>( &ivar_hdu ), " ivar hdu in input arc fits")
    ( "first_bundle", popts::value<int>( &first_fiber_bundle ), "first fiber bundle to fit")
    ( "last_bundle", popts::value<int>( &last_fiber_bundle ), "last fiber bundle to fit")
    ( "first_fiber", popts::value<int>( &first_fiber ), "first fiber (must be in bundle)")
    ( "last_fiber", popts::value<int>( &last_fiber ), "last fiber (must be in bundle)")
    ( "half_size_x", popts::value<int>( &half_size_x ), "half size of PSF stamp (full size is 2*half_size+1)")
    ( "half_size_y", popts::value<int>( &half_size_y ), "half size of PSF stamp (full size is 2*half_size+1)")
    ( "psfmodel", popts::value<string>( &psf_model ), "PSF model, default is GAUSSHERMITE")
    ( "positions", "fit positions of each spot individually after global fit for debugging")
    ( "verbose,v", "turn on verbose mode" )
    ( "lamplines", popts::value<string>( &lamp_lines_filename ), "lamp lines ASCII file name (def. is $SPECEXDATA/opfiles/lamplines.par)" )
    ( "core", "dump core files when harp exception is thrown" )
    ( "gauss_hermite_deg",  popts::value<int>( &gauss_hermite_deg ), "degree of Hermite polynomials (same for x and y, only if GAUSSHERMITE psf)")
    ("gauss_hermite_deg2",  popts::value<int>( &gauss_hermite_deg2 ), "degree of Hermite polynomials (same for x and y, only if GAUSSHERMITE2 psf)")
    //( "gauss_hermite_sigma",  popts::value<double>( &gauss_hermite_sigma ), "sigma of Gauss-Hermite PSF (same for x and y, only if GAUSSHERMITE psf)")
    ( "legendre_deg_wave",  popts::value<int>( &legendre_deg_wave ), "degree of Legendre polynomials along wavelength (can be reduced if missing data)")
    ( "legendre_deg_x",  popts::value<int>( &legendre_deg_x ), "degree of Legendre polynomials along x_ccd (can be reduced if missing data)")
    ( "psf_error",  popts::value<double>( &psf_error ), "psf fractional uncertainty (default is 0.01, for weights in the fit)")
    ( "psf_core_wscale",  popts::value<double>( &psf_core_wscale ), "scale up the weight of pixels in 5x5 PSF core")
#ifdef EXTERNAL_TAIL
    ( "fit_psf_tails", "unable fit of psf tails")
#endif
#ifdef CONTINUUM
    ( "fit_continuum", "unable fit of continuum")
#endif
    ( "no_trace_fit", "do not fit traces")
    ( "out_xml", popts::value<string>( &output_xml_filename ), " output psf xml file name")
    ( "out_fits", popts::value<string>( &output_fits_filename ), " output psf fits file name")  
    ( "out_spots", popts::value<string>( &output_spots_filename ), " output spots file name")  
    //( "out", popts::value<string>( &outfile ), "output image file" )
    ;

  popts::variables_map vm;

  try {
    
    popts::store(popts::command_line_parser( argc, argv ).options(desc).run(), vm);
    popts::notify(vm);
    
    if ( ( argc < 2 ) || vm.count( "help" ) || ( ! vm.count( "arc" ) ) ) {
      cerr << endl;
      cerr << desc << endl;
      cerr << "example:" << endl;
      cerr << argv[0] << " -v" << endl;
      return EXIT_FAILURE;
    }

    if(lamp_lines_filename == "") {
      cerr << endl;
      cerr << "missing lamp_lines_filename either define env. variable SPECEXDATA or use option --lamplines" << endl;
      return EXIT_FAILURE;
    }


  }catch(std::exception e) {
    cerr << "error in arguments" << endl;
    cerr << endl;
    cerr << desc << endl;
    return EXIT_FAILURE;
  }
  
  try {
    specex_set_verbose(vm.count("verbose")>0);
    specex_set_dump_core(vm.count("core")>0);
    bool fit_traces = (vm.count("no_trace_fit")==0);
#ifdef EXTERNAL_TAIL
    bool fit_psf_tails = (vm.count("fit_psf_tails")>0);
#endif
#ifdef CONTINUUM
    bool fit_continuum = (vm.count("fit_continuum")>0);
#endif
    bool fit_individual_spots_position = vm.count("positions");
    
    SPECEX_INFO("using lamp lines file " << lamp_lines_filename); 
    
  

    // define spectrograph (this should be improved)
    // --------------------------------------------
    TraceSet traceset;
    map<string,string> image_infos;
    Spectrograph *spectro = 0;
    if(spectrograph_name == "BOSS") {
      SPECEX_ERROR("do not deal with BOSS here");
    }else if(spectrograph_name == "DESI"){
      spectro = new DESI_Spectrograph();
      /* read traces in arc file */
      read_DESI_traceset_in_fits(traceset,arc_image_filename,5,6);
      spectro->AutoConfigure(traceset);
      read_DESI_keywords(arc_image_filename,image_infos);
    }else{
      SPECEX_ERROR("unknown spectrograph");
    }
    
    

    
    if(first_fiber_bundle<0 || first_fiber_bundle>= spectro->number_of_fiber_bundles_per_ccd) {
      SPECEX_ERROR("invalid first fiber bundle");
    }
    if(last_fiber_bundle<first_fiber_bundle || last_fiber_bundle>= spectro->number_of_fiber_bundles_per_ccd) {
      SPECEX_ERROR("invalid last fiber bundle");
    }
    
  
    
    // open image
    // --------------------------------------------
    image_data image,weight;
    //read_fits_images(arc_image_filename,image,weight);
    read_fits_image(arc_image_filename,flux_hdu,image);
    read_fits_image(arc_image_filename,ivar_hdu,weight);
    
    SPECEX_INFO("image  size = " << image.n_cols() << "x" << image.n_rows());
    SPECEX_INFO("weight size = " << weight.n_cols() << "x" << weight.n_rows());
    

    // weight is 0 or 1
    /*
      for(int j=0;j<weight.Ny();j++) {
      for(int i=0;i<weight.Nx();i++) {
	if(weight(i,j)>0)
	  weight(i,j)=1;
	else
	  weight(i,j)=0;
      }
    }
    */
    
    // init PSF
    // --------------------------------------------
    specex::PSF_p psf;

    if(psf_model=="GAUSSHERMITE")
      psf = PSF_p(new specex::GaussHermitePSF(gauss_hermite_deg));
    else if(psf_model=="GAUSSHERMITE2")
      psf = PSF_p(new specex::GaussHermite2PSF(gauss_hermite_deg,gauss_hermite_deg2));
    else if(psf_model=="HATHERMITE")
      psf = PSF_p(new specex::HatHermitePSF(gauss_hermite_deg));
    else if(psf_model=="HATMOFFAT")
      psf = PSF_p(new specex::HatMoffatPSF(gauss_hermite_deg));
    else
      SPECEX_ERROR("don't know this psf model");
    
    psf->ccd_image_n_cols = image.n_cols();
    psf->ccd_image_n_rows = image.n_rows();
    psf->hSizeX = half_size_x;
    psf->hSizeY = half_size_y;
    
        
    
    psf->FiberTraces.clear();
    
    
    // init PSF fitter
    // -------------------------------------------- 
    PSF_Fitter fitter(psf,image,weight);
    
    fitter.polynomial_degree_along_x    = legendre_deg_x;
    fitter.polynomial_degree_along_wave = legendre_deg_wave;
    fitter.psf->psf_error               = psf_error;
    fitter.corefootprint_weight_boost   = psf_core_wscale;
    
#ifdef EXTERNAL_TAIL
    fitter.scheduled_fit_of_psf_tail    = fit_psf_tails;
#endif

#ifdef CONTINUUM
    fitter.scheduled_fit_of_continuum   = fit_continuum;
#endif
    
    fitter.scheduled_fit_of_traces      = fit_traces;

    fitter.psf->gain = 1; // images are already in electrons
    // fitter.readout_noise = 2; // b1, evaluated using spatial variance in sdProc-b1-00108382.fits[1:4000,1:600]

    // compute mean readout noise
    fitter.psf->readout_noise = 2;
    if(image_infos.find("RDNOISE") != image_infos.end()) {
      fitter.psf->readout_noise = atof(image_infos["RDNOISE"].c_str());
      if(fitter.psf->readout_noise<=0) {
	SPECEX_ERROR("Absurd RDNOISE read in header = " << fitter.psf->readout_noise );
      }
      SPECEX_INFO("Using RDNOISE from image header = " << fitter.psf->readout_noise );
    }
    
    
    for(int j=0;j<weight.Ny();j++) {
      for(int i=0;i<weight.Nx();i++) {
	if(weight(i,j)>0)
	  weight(i,j)=1/square(psf->readout_noise);
	else
	  weight(i,j)=0;
      }
    }

    fitter.mask.Clear();

    
    SPECEX_INFO(
		"PSF '" << psf_model << "' stamp size = " 
		<< psf->hSizeX << "x" << psf->hSizeY );
    
    
    
    int first_fitted_fiber=-1;
    int last_fitted_fiber=-1;
    vector<Spot_p> fitted_spots;
    
    // loop on fiber bundles 
    // -------------------------------------------- 
    for(int bundle = first_fiber_bundle; bundle <= last_fiber_bundle ; bundle ++) {
    
      // allocate bundle in PSF if necessary
      if(psf->ParamsOfBundles.find(bundle)==psf->ParamsOfBundles.end()) {
	psf->ParamsOfBundles[bundle] = specex::PSF_Params();
	psf->ParamsOfBundles[bundle].bundle_id = bundle;
	psf->ParamsOfBundles[bundle].fiber_min = spectro->number_of_fibers_per_bundle*bundle;
	psf->ParamsOfBundles[bundle].fiber_max = psf->ParamsOfBundles[bundle].fiber_min+spectro->number_of_fibers_per_bundle-1; // included
	
	// now check mask ?


	if(psf->ParamsOfBundles[bundle].fiber_min < first_fiber) {
	  psf->ParamsOfBundles[bundle].fiber_min = first_fiber;
	  SPECEX_INFO("restricting fiber range first fiber = " << first_fiber);
	}
	if(psf->ParamsOfBundles[bundle].fiber_max > last_fiber) {
	  psf->ParamsOfBundles[bundle].fiber_max = last_fiber;
	  SPECEX_INFO("restricting fiber range last fiber = " << last_fiber);
	}
	
      }
      
      // load traces in PSF 
      // --------------------------------------------  
      for(int fiber=psf->ParamsOfBundles[bundle].fiber_min; fiber<=psf->ParamsOfBundles[bundle].fiber_max; fiber++) {
	psf->FiberTraces[fiber]=traceset[fiber];
      } 
      
      fitter.SelectFiberBundle(bundle);
      
      
      // loading arc lamp spots belonging to this bundle
      // --------------------------------------------
      
      int ymin = 0; // range of usable CCD coordinates, hard coded for now
      int ymax = image.n_rows(); // range of usable CCD coordinates, hard coded for now
      
      if(psf->camera_id=="b1") {ymin=696; ymax = 3516;};
      if(psf->camera_id=="b2") {ymin=696; ymax = 3516;};
      if(psf->camera_id=="r1") {ymin=200; ymax = 3668;}; 
      if(psf->camera_id=="r2") {ymin=200; ymax = 3668;};

      
      int margin = 1; // we want to keep as much as possible to minimize extrapolation
      ymin+=margin;
      ymax-=margin;
      

      SPECEX_INFO("valid y(=rows) range = " << ymin << " " << ymax);

      vector<Spot_p> spots;
      allocate_spots_of_bundle(spots,*spectro,lamp_lines_filename,traceset,bundle,psf->ParamsOfBundles[bundle].fiber_min,psf->ParamsOfBundles[bundle].fiber_max,ymin,ymax,min_wavelength,max_wavelength);
      SPECEX_INFO("number of spots = " << spots.size());
      
      if(write_tmp_results)
	write_spots_xml(spots,"spots-init.xml");
      
      
      //exit(12);

      // starting fit
      // --------------------------------------------
      fitter.FitEverything(spots,true);

      int ndf = psf->ParamsOfBundles[bundle].ndata - psf->ParamsOfBundles[bundle].nparams;
      SPECEX_INFO("Bundle " << bundle << " PSF fit status   = "<< psf->ParamsOfBundles[bundle].fit_status);
      SPECEX_INFO("Bundle " << bundle << " PSF fit chi2/ndf = "<< psf->ParamsOfBundles[bundle].chi2 << "/" << ndf << " = " << psf->ParamsOfBundles[bundle].chi2/ndf);
      SPECEX_INFO("Bundle " << bundle << " PSF fit ndata    = "<< psf->ParamsOfBundles[bundle].ndata);
      SPECEX_INFO("Bundle " << bundle << " PSF fit nspots   = "<< psf->ParamsOfBundles[bundle].nspots_in_fit);
      SPECEX_INFO("Bundle " << bundle << " PSF fit chi2/ndata (core) = "<< psf->ParamsOfBundles[bundle].chi2_in_core << "/" << psf->ParamsOfBundles[bundle].ndata_in_core << " = " << psf->ParamsOfBundles[bundle].chi2_in_core/psf->ParamsOfBundles[bundle].ndata_in_core);
      

      if(bundle == first_fiber_bundle) {
	first_fitted_fiber=psf->ParamsOfBundles[bundle].fiber_min;
	last_fitted_fiber=psf->ParamsOfBundles[bundle].fiber_max;
      }
      
      first_fitted_fiber=min(first_fitted_fiber,psf->ParamsOfBundles[bundle].fiber_min);
      last_fitted_fiber=max(last_fitted_fiber,psf->ParamsOfBundles[bundle].fiber_max);
      
      
      for(size_t s=0;s<spots.size();s++) {
	fitted_spots.push_back(spots[s]);
      }


      

      if(fit_individual_spots_position) // for debugging
	fitter.FitIndividualSpotPositions(spots);

      

    } // end of loop on bundles
    
    
    if(output_xml_filename != "")
      write_psf_xml(fitter.psf,output_xml_filename);
    if(output_fits_filename != "")
      write_psf_fits(fitter.psf,output_fits_filename);
    if(output_spots_filename != "")
      write_spots_xml(fitted_spots,output_spots_filename);
    


    //write_psf_fits_image(fitter.psf,"psf.fits",500,2000,1,4);


    
    
    
    
    
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
  
  return EXIT_SUCCESS;
}



