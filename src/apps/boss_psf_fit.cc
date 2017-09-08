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
#include <specex_boss_io.h>
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
#include <portable_fenv.h>



int main ( int argc, char *argv[] ) {
  
  // to crash when NaN
  feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
  

  // default arguments
  // --------------------------------------------
  string psf_model = "GAUSSHERMITE";
  string spectrograph_name = "BOSS";
  int    first_fiber_bundle=1;
  int    last_fiber_bundle=1;
  int    first_fiber=0;
  int    last_fiber=100000;
  int    half_size_x=4;
  int    half_size_y=4;
  
  string arc_image_filename="";
  string xy_trace_fits_name="";
  string wy_trace_fits_name="";
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
  

  if(getenv("IDLSPEC2D_DIR")) 
    lamp_lines_filename = string(getenv("IDLSPEC2D_DIR"))+"/opfiles/lamplines.par";
   
  // reading arguments
  // --------------------------------------------
  popts::options_description desc ( "Allowed Options" );
  desc.add_options()
    ( "help,h", "display usage information" )
    ( "arc,a", popts::value<string>( &arc_image_filename ), "arc pre-reduced fits image file name (mandatory), ex:  sdProc-b1-00108382.fits" )
    ( "xy", popts::value<string>( &xy_trace_fits_name ), " fits file name where is stored x vs y trace (mandatory), ex: spFlat-b1-00108381.fits.gz" )
    ( "wy", popts::value<string>( &wy_trace_fits_name ), " fits file name where is stored w vs y trace (mandatory), ex: spArc-b1-00108382.fits.gz" )
    ( "first_bundle", popts::value<int>( &first_fiber_bundle ), "first fiber bundle to fit")
    ( "last_bundle", popts::value<int>( &last_fiber_bundle ), "last fiber bundle to fit")
    ( "first_fiber", popts::value<int>( &first_fiber ), "first fiber (must be in bundle)")
    ( "last_fiber", popts::value<int>( &last_fiber ), "last fiber (must be in bundle)")
    ( "half_size_x", popts::value<int>( &half_size_x )->default_value(4), "half size of PSF stamp (full size is 2*half_size+1)")
    ( "half_size_y", popts::value<int>( &half_size_y )->default_value(4), "half size of PSF stamp (full size is 2*half_size+1)")
    ( "psfmodel", popts::value<string>( &psf_model ), "PSF model, default is GAUSSHERMITE")
    ( "positions", "fit positions of each spot individually after global fit for debugging")
    ( "verbose,v", "turn on verbose mode" )
    ( "lamplines", popts::value<string>( &lamp_lines_filename ), "lamp lines ASCII file name (def. is $IDLSPEC2D_DIR/opfiles/lamplines.par)" )
    ( "core", "dump core files when harp exception is thrown" )
    ( "gauss_hermite_deg",  popts::value<int>( &gauss_hermite_deg )->default_value(3), "degree of Hermite polynomials (same for x and y, only if GAUSSHERMITE psf)")
    ("gauss_hermite_deg2",  popts::value<int>( &gauss_hermite_deg2 )->default_value(2), "degree of Hermite polynomials (same for x and y, only if GAUSSHERMITE2 psf)")
    //( "gauss_hermite_sigma",  popts::value<double>( &gauss_hermite_sigma ), "sigma of Gauss-Hermite PSF (same for x and y, only if GAUSSHERMITE psf)")
    ( "legendre_deg_wave",  popts::value<int>( &legendre_deg_wave )->default_value(4), "degree of Legendre polynomials along wavelength (can be reduced if missing data)")
    ( "legendre_deg_x",  popts::value<int>( &legendre_deg_x )->default_value(1), "degree of Legendre polynomials along x_ccd (can be reduced if missing data)")
    ( "psf_error",  popts::value<double>( &psf_error )->default_value(0), "psf fractional uncertainty (default is 0.01, for weights in the fit)")
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
    //( "out", popts::value<string>( &outfile ), "output image file" )
    ;

  popts::variables_map vm;

  try {
    
    popts::store(popts::command_line_parser( argc, argv ).options(desc).run(), vm);
    popts::notify(vm);
    
    if ( ( argc < 2 ) || vm.count( "help" ) || ( ! vm.count( "arc" ) )  || ( ! vm.count( "xy" ) ) || ( ! vm.count( "wy" ) ) ) {
      cerr << endl;
      cerr << desc << endl;
      cerr << "example:" << endl;
      cerr << argv[0] << " --arc sdProc-b1-00108382.fits --xy redux/spFlat-b1-00108381.fits.gz --wy redux/spArc-b1-00108382.fits.gz -v" << endl;
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
    Spectrograph *spectro = 0;
    if(spectrograph_name == "BOSS") {
      spectro = new BOSS_Spectrograph();
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
    read_fits_images(arc_image_filename,image,weight);

    SPECEX_WARNING("Hardcoded readnoise for BOSS, need to change this to get reliable errors");
    image_data readnoise(image.n_cols(),image.n_rows());
    for(int j=0;j<readnoise.Ny();j++) {
      for(int i=0;i<readnoise.Nx();i++) {
	readnoise(i,j)=2;
      }
    }
    
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
    
        
    

    // read trace set derived in BOSS pipeline (this should be improved)
    // --------------------------------------------    
    TraceSet traceset;
    
    map<string,string> image_infos;
    
    if(spectrograph_name == "BOSS") {
      int xy_trace_hdu = 1; // in the spFlat file , warning : this HDU numbering starts at 0 (to check)
      int wy_trace_hdu = 2; // int the spArc file , warning : this HDU numbering starts at 0 (to check)
      read_BOSS_traceset_in_fits(traceset,wy_trace_fits_name,wy_trace_hdu,xy_trace_fits_name,xy_trace_hdu);
      
      // SPECEX_INFO("EXIT FOR DEBUG"); return 0;
      read_BOSS_keywords(psf,arc_image_filename,image_infos);
    }
    
    
    
    psf->FiberTraces.clear();
    
    
    // init PSF fitter
    // --------------------------------------------

    
    
    PSF_Fitter fitter(psf,image,weight,readnoise);
    
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
    fitter.psf->readout_noise = 0; // property of image, not PSF
        
    for(int j=0;j<weight.Ny();j++) {
      for(int i=0;i<weight.Nx();i++) {
	if(weight(i,j)>0)
	  weight(i,j)=1/square(readnoise(i,j));
	else
	  weight(i,j)=0;
      }
    }

    fitter.mask.Clear();

    /*
    SPECEX_INFO(
		"PSF '" << psf_model << "' stamp size = " 
		<< psf->hSizeX << "x" << psf->hSizeY << " npar(loc) = " 
		<< psf->FixedCoordNPar() << " npar(glob) = " << psf->VaryingCoordNPar(bundle_id)
		);
    */
    
    
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
      int ymax = 0; // range of usable CCD coordinates, hard coded for now
      
      if(psf->camera_id=="b1") {ymin=696; ymax = 3516;};
      if(psf->camera_id=="b2") {ymin=696; ymax = 3516;};
      if(psf->camera_id=="r1") {ymin=200; ymax = 3668;}; 
      if(psf->camera_id=="r2") {ymin=200; ymax = 3668;};

#warning "need to know where to get xy ccd range info"

      int margin = 1; // we want to keep as much as possible to minimize extrapolation
      ymin+=margin;
      ymax-=margin;
      

      vector<Spot_p> spots;
      allocate_spots_of_bundle(spots,lamp_lines_filename,traceset,bundle,psf->ParamsOfBundles[bundle].fiber_min,psf->ParamsOfBundles[bundle].fiber_max,ymin,ymax,min_wavelength,max_wavelength);
      SPECEX_INFO("number of spots = " << spots.size());
      
      
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

      {
	// writing spots as xml
	
	char filename[100];
	sprintf(filename,"spots-%s-%08d-%03d-%03d.xml",psf->camera_id.c_str(),(int)psf->arc_exposure_id,first_fitted_fiber,last_fitted_fiber);
	write_spots_xml(fitted_spots,filename);
      }
      {
	// writing psf as xml
	char filename[100];
	
	sprintf(filename,"psf-%s-%08d-%03d-%03d.xml",psf->camera_id.c_str(),(int)psf->arc_exposure_id,first_fitted_fiber,last_fitted_fiber);
	write_psf_xml(fitter.psf,filename);
	
	//if(psf_model == "GAUSSHERMITE" || psf_model == "GAUSSHERMITE2") {
	if(psf_model == "GAUSSHERMITE2") {
	  sprintf(filename,"psf-%s-%08d-%03d-%03d.fits",psf->camera_id.c_str(),(int)psf->arc_exposure_id,first_fitted_fiber,last_fitted_fiber);
	  write_psf_fits(fitter.psf,filename);
	}
	

      }

      if(fit_individual_spots_position) // for debugging
	fitter.FitIndividualSpotPositions(spots);

      

    } // end of loop on bundles
    
    
    if(output_xml_filename != "")
      write_psf_xml(fitter.psf,output_xml_filename);
    if(output_fits_filename != "")
      write_psf_fits(fitter.psf,output_fits_filename);
    


    //write_psf_fits_image(fitter.psf,"psf.fits",500,2000,1,4);


    
    
    
    
    
  // ending
  // --------------------------------------------
  } catch(harp::exception e) {
    cerr << "FATAL ERROR (harp) " << e.what() << endl;
    return EXIT_FAILURE;
  }
  /*
  catch (std::exception& e) {
    cerr << "FATAL ERROR (other std) " << e.what() <<endl;
    return EXIT_FAILURE;
  }catch (...) {
    cerr << "FATAL ERROR (unknown)" << endl;
    return EXIT_FAILURE;
  }
  */
  return EXIT_SUCCESS;
}



