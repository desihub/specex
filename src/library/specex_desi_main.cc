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
#include <specex_spectrograph.h>
#include <specex_lamp_lines_utils.h>
#include <specex_fits.h>
#include <specex_desi_io.h>
#include <specex_psf_fitter.h>

#include <specex_psf_io.h>


using namespace std;
using namespace specex;

namespace popts = boost::program_options;

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



int specex_desi_psf_fit_main ( int argc, char *argv[] ) {
  
  // to crash when NaN
  feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
  
  // default arguments
  // --------------------------------------------
  string psf_model = "GAUSSHERMITE";
  string spectrograph_name = "DESI";  
  int    first_fiber_bundle=0;
  int    last_fiber_bundle=0;
  int    first_fiber=0;
#define DEFAULT_LAST_FIBER 100000
  int    last_fiber=DEFAULT_LAST_FIBER;
  int    half_size_x=4;
  int    half_size_y=4;
  
  string arc_image_filename="";
  string trace_filename="";
  string xtrace_filename="";
  string ytrace_filename="";
  int xtrace_hdu=0;
  int ytrace_hdu=0;
  //string wy_trace_fits_name="";
  string lamp_lines_filename="";  
  double min_wavelength = 0;
  double max_wavelength = 1e6;
  int gauss_hermite_deg  = 3;
  int gauss_hermite_deg2 = 2;
  //double gauss_hermite_sigma = 1.1;
  int legendre_deg_wave = 4;
  int legendre_deg_x = 1;
  int trace_deg_wave = 0;
  int trace_deg_x = 0;
  
  double psf_error = 0;
  double psf_core_wscale = 0;
  
  string output_xml_filename="";
  string output_fits_filename="";
  string output_spots_filename="";
  string input_psf_filename="";
  
  int flux_hdu=0;
  int ivar_hdu=0;
  int mask_hdu=0;
  int header_hdu=0;
  
  if(getenv("SPECEXDATA"))
    lamp_lines_filename = string(getenv("SPECEXDATA"))+"/specex_linelist_desi.txt";
   
  // reading arguments
  // --------------------------------------------
  popts::options_description desc ( "Allowed Options" );

  vector<string> argurment_priors;

  desc.add_options()
    ( "help,h", "display usage information" )
    ( "arc,a", popts::value<string>( &arc_image_filename ), "arc pre-reduced fits image file name (mandatory), ex:  sdProc-b1-00108382.fits" )
    ( "flux-hdu", popts::value<int>( &flux_hdu )->default_value(1), " flux hdu in input arc fits")
    ( "ivar-hdu", popts::value<int>( &ivar_hdu )->default_value(2), " ivar hdu in input arc fits")
    ( "mask-hdu", popts::value<int>( &mask_hdu )->default_value(3), " mask hdu in input arc fits")
    ( "header-hdu", popts::value<int>( &header_hdu )->default_value(1), " header hdu in input arc fits")

    ( "trace", popts::value<string>( &trace_filename ), "fits image file name with image extensions XTRACE and YTRACE for legendre polynomial of wavelength (mandatory)" )
    
    ( "xcoord-file", popts::value<string>( &xtrace_filename ), "fits image file name with x trace legendre polynomial of wavelength (ignored if --trace given)" )
    ( "xcoord-hdu", popts::value<int>( &xtrace_hdu )->default_value(-1), "hdu of x trace legendre polynomial of wavelength (default is looking for extension XTRACE)" )
    ( "ycoord-file", popts::value<string>( &ytrace_filename ), "fits image file name with y trace legendre polynomial of wavelength (ignored if --trace given)" )
    ( "ycoord-hdu", popts::value<int>( &ytrace_hdu )->default_value(-1), "hdu of y trace legendre polynomial of wavelength (default is looking for extension YTRACE)" )
    ( "first_bundle", popts::value<int>( &first_fiber_bundle ), "first fiber bundle to fit")
    ( "last_bundle", popts::value<int>( &last_fiber_bundle ), "last fiber bundle to fit")
    ( "first_fiber", popts::value<int>( &first_fiber ), "first fiber (must be in bundle)")
    ( "last_fiber", popts::value<int>( &last_fiber ), "last fiber (must be in bundle)")
    ( "half_size_x", popts::value<int>( &half_size_x )->default_value(8), "half size of PSF stamp (full size is 2*half_size+1)")
    ( "half_size_y", popts::value<int>( &half_size_y )->default_value(5), "half size of PSF stamp (full size is 2*half_size+1)")
    ( "psfmodel", popts::value<string>( &psf_model )->default_value("GAUSSHERMITE"), "PSF model, default is GAUSSHERMITE")
    ( "positions", "fit positions of each spot individually after global fit for debugging")
    ( "verbose,v", "turn on verbose mode (deprecated, true by default)" )
    ( "quiet", "no info message, only warning" )
    ( "debug", "turn on debug mode" )
    ( "lamplines", popts::value<string>( &lamp_lines_filename ), "lamp lines ASCII file name (def. is $SPECEXDATA/specex_linelist_desi.txt)" )
    ( "core", "dump core files when harp exception is thrown" )
    ( "gauss_hermite_deg",  popts::value<int>( &gauss_hermite_deg )->default_value(6), "degree of Hermite polynomials (same for x and y, only if GAUSSHERMITE psf)")
    ("gauss_hermite_deg2",  popts::value<int>( &gauss_hermite_deg2 )->default_value(2), "degree of Hermite polynomials (same for x and y, only if GAUSSHERMITE2 psf)")
    //( "gauss_hermite_sigma",  popts::value<double>( &gauss_hermite_sigma ), "sigma of Gauss-Hermite PSF (same for x and y, only if GAUSSHERMITE psf)")
    ( "legendre_deg_wave",  popts::value<int>( &legendre_deg_wave )->default_value(4), "degree of Legendre polynomials along wavelength (can be reduced if missing data)")
    ( "legendre_deg_x",  popts::value<int>( &legendre_deg_x )->default_value(1), "degree of Legendre polynomials along x_ccd (can be reduced if missing data)")
    ( "trace_deg_wave",  popts::value<int>( &trace_deg_wave )->default_value(6), "degree of Legendre polynomials along wavelength for fit of traces")
    ( "trace_deg_x",  popts::value<int>( &trace_deg_x )->default_value(6), "degree of Legendre polynomials along x_ccd for fit of traces")
    ( "psf_error",  popts::value<double>( &psf_error )->default_value(0), "psf fractional uncertainty (default is 0.01, for weights in the fit)")
    ( "psf_core_wscale",  popts::value<double>( &psf_core_wscale ), "scale up the weight of pixels in 5x5 PSF core")
    ( "per_fiber", "as many X parameters as fibers")
#ifdef EXTERNAL_TAIL
    ( "fit_psf_tails", "unable fit of psf tails")
#endif
#ifdef CONTINUUM
    ( "fit_continuum", "unable fit of continuum")
#endif
    ( "variance_model", "refit at the end with a model of the variance to avoid Poisson noise bias")
    ( "no_trace_fit", "do not fit traces")
    ( "no_sigma_fit", "do not fit the gaussian sigma")
    ( "no_psf_fit", "do not fit the psf")
    ( "out_xml", popts::value<string>( &output_xml_filename ), " output psf xml file name")
    ( "out_fits", popts::value<string>( &output_fits_filename ), " output psf fits file name")  
    ( "input-psf", popts::value<string>( &input_psf_filename ), " input psf file name (fits or xml)")  
    ( "out_spots", popts::value<string>( &output_spots_filename ), " output spots file name")  
    ( "prior", popts::value< vector<string> >( &argurment_priors )->multitoken(), " gaussian prior on a param : 'name' value error")  
    ( "tmp_results", " write tmp results")  
    //( "out", popts::value<string>( &outfile ), "output image file" )
    ;

  popts::variables_map vm;

  try {
    
    popts::store(popts::command_line_parser( argc, argv ).options(desc).run(), vm);
    popts::notify(vm);
    
    if ( ( argc < 2 ) || vm.count( "help" ) || ( ! vm.count( "arc" ) ) || ( (! vm.count("trace")) && (( ! vm.count( "xcoord-file" )  || ( ! vm.count( "ycoord-file" ) ) ) ) && ( ! vm.count( "input-psf" ) ))) {
      cerr << endl;
      cerr << desc << endl;
      cerr << "example:" << endl;
      cerr << argv[0] << " -v" << endl;
      return EXIT_FAILURE;
    }
    if(vm.count("trace") && vm.count( "input-psf" )) {
      cout << "will use input-psf and ignore --trace option" << endl;
    }
    if(vm.count( "input-psf" )) {
      if(first_fiber>0 || last_fiber<DEFAULT_LAST_FIBER) {
	cerr << "cannot specify first and last fiber with input psf option, because this would mean to change number of parameters" << endl;
	exit(12);
      }
    }
    
    if(lamp_lines_filename == "") {
      cerr << endl;
      cerr << "missing lamp_lines_filename either define env. variable SPECEXDATA or use option --lamplines" << endl;
      return EXIT_FAILURE;
    }


  }catch(std::exception& e) {
    cerr << "error in arguments : " << e.what() << endl;
    cerr << endl;
    cerr << desc << endl;
    return EXIT_FAILURE;
  }
  
  // dealing with priors
  map<string,Prior*> priors;
  {
    int np=int(argurment_priors.size());
    if( ! (np%3==0) ) {
      cerr << "error in parsing, priors must be of the form 'name value error' (np=" << np << ")" << endl;
      cerr << desc << endl;
      return EXIT_FAILURE;
    }
    try {
      for(int i=0;i<np/3;i++) {
	string pname=argurment_priors[3*i];
	double val=atof(argurment_priors[3*i+1].c_str());
	double err=atof(argurment_priors[3*i+2].c_str());
	cout << "priors[" << i << "]= " << pname << " " << val << " " << err << endl;
	priors[pname]= new GaussianPrior(val,err);
      }
    }catch(std::exception e) {
      cerr << "error in parsing arguments of priors" << endl;
      cerr << "priors must be of the form 'name value error'" << endl;
      cerr << desc << endl;
      return EXIT_FAILURE;
    } 
  }
  


  try {
    specex_set_debug(vm.count("debug")>0);
    //specex_set_verbose(vm.count("verbose")>0);
    specex_set_verbose(vm.count("quiet")==0);
    
    specex_set_dump_core(vm.count("core")>0);
    bool fit_traces = (vm.count("no_trace_fit")==0);
    bool fit_sigmas = (vm.count("no_sigma_fit")==0);
    bool fit_psf    = (vm.count("no_psf_fit")==0);
    bool write_tmp_results = (vm.count("tmp_results")>0);
    
#ifdef EXTERNAL_TAIL
    bool fit_psf_tails = (vm.count("fit_psf_tails")>0);
#endif
#ifdef CONTINUUM
    bool fit_continuum = (vm.count("fit_continuum")>0);
#endif
    bool use_variance_model = (vm.count("variance_model")>0);
    bool fit_individual_spots_position = vm.count("positions");
    
    SPECEX_INFO("using lamp lines file " << lamp_lines_filename); 
    
    bool per_fiber=vm.count("per_fiber")>0;
    
      

    // read data
    // --------------------------------------------
    map<string,string> header;

    image_data image,weight,mask,rdnoise;
    read_DESI_preprocessed_image(arc_image_filename,image,weight,mask,rdnoise,header);
    
    Spectrograph *spectro = new DESI_Spectrograph();
    
    // read trace
    // --------------------------------------------
    /* read traces in arc file */
    if(trace_filename !="") {
      xtrace_filename=trace_filename;
      ytrace_filename=trace_filename;
    }

    specex::PSF_p psf;
    TraceSet traceset; // only filled if no input PSF
    if( input_psf_filename == "") {

      read_DESI_traceset_in_fits(traceset,xtrace_filename,xtrace_hdu,ytrace_filename,ytrace_hdu,trace_deg_x,trace_deg_wave);
      if(per_fiber) {
	SPECEX_INFO("Use as many parameters along X as there are fibers (solve all at once, ie a single bundle, and overwrite legendre_deg_x)");
	int nfibers = traceset.size();
	spectro->number_of_fiber_bundles_per_ccd=nfibers;
	spectro->number_of_fibers_per_bundle=nfibers;
	SPECEX_DEBUG("number of fibers  = " << nfibers);
	legendre_deg_x = nfibers - 1;
      }else{
	spectro->AutoConfigure(traceset);
      }
      
      
      if(first_fiber_bundle<0 || first_fiber_bundle>= spectro->number_of_fiber_bundles_per_ccd) {
	SPECEX_ERROR("invalid first fiber bundle");
      }
      if(last_fiber_bundle<first_fiber_bundle || last_fiber_bundle>= spectro->number_of_fiber_bundles_per_ccd) {
	SPECEX_ERROR("invalid last fiber bundle");
      }
      
      // init PSF
      // --------------------------------------------
      
      if(psf_model=="GAUSSHERMITE")
	psf = PSF_p(new specex::GaussHermitePSF(gauss_hermite_deg));
      else if(psf_model=="GAUSSHERMITE2")
	psf = PSF_p(new specex::GaussHermite2PSF(gauss_hermite_deg,gauss_hermite_deg2));
      else if(psf_model=="HATHERMITE")
	psf = PSF_p(new specex::HatHermitePSF(gauss_hermite_deg));
      else if(psf_model=="HATMOFFAT")
	psf = PSF_p(new specex::HatMoffatPSF(gauss_hermite_deg));
      else if(psf_model=="DISKMOFFAT")
	psf = PSF_p(new specex::DiskMoffatPSF(gauss_hermite_deg));
      else
	SPECEX_ERROR("don't know this psf model");
      
      psf->ccd_image_n_cols = image.n_cols();
      psf->ccd_image_n_rows = image.n_rows();
      psf->hSizeX = half_size_x;
      psf->hSizeY = half_size_y;
      psf->FiberTraces.clear();
    }else{
      read_psf(psf,input_psf_filename);

      // erase traces outside fiber range
      /*
	for (specex::TraceSet::iterator it = psf->FiberTraces.begin(); it != psf->FiberTraces.end() ; ) {
	if((it->first<first_fiber)||(it->first>last_fiber)) {
	SPECEX_INFO("erase trace of fiber " << it->first);
	psf->FiberTraces.erase(it++);
	}else {
	++it;
	}
	}
      */

      
      if(per_fiber) {
	SPECEX_INFO("Use as many parameters along X as there are fibers (solve all at once, ie a single bundle, and overwrite legendre_deg_x)");
	int nfibers = psf->FiberTraces.size();
	spectro->number_of_fiber_bundles_per_ccd=nfibers;
	spectro->number_of_fibers_per_bundle=nfibers;
	SPECEX_DEBUG("number of fibers  = " << nfibers);
	legendre_deg_x = nfibers - 1;
      }else{
	spectro->AutoConfigure(psf->FiberTraces);
      }
      
    }
    

    
    // init PSF fitter
    // -------------------------------------------- 
    PSF_Fitter fitter(psf,image,weight,rdnoise);
    
    fitter.polynomial_degree_along_x    = legendre_deg_x;
    fitter.polynomial_degree_along_wave = legendre_deg_wave;
    fitter.psf->psf_error               = psf_error;
    fitter.corefootprint_weight_boost   = psf_core_wscale;
    fitter.write_tmp_results            = write_tmp_results;
    
#ifdef EXTERNAL_TAIL
    fitter.scheduled_fit_of_psf_tail    = fit_psf_tails;
#endif

#ifdef CONTINUUM
    fitter.scheduled_fit_of_continuum   = fit_continuum;
#endif

    fitter.scheduled_fit_with_weight_model  = use_variance_model;
    
    fitter.scheduled_fit_of_traces      = fit_traces;
    fitter.scheduled_fit_of_sigmas      = fit_sigmas;
    fitter.scheduled_fit_of_psf         = fit_psf;
    fitter.direct_simultaneous_fit      = (input_psf_filename != "");
    
    fitter.psf->gain = 1; // images are already in electrons
    fitter.psf->readout_noise = 0; // readnoise is a property of image, not PSF
        
    if(header.find("CAMERA") != header.end()) {
      fitter.psf->camera_id = header["CAMERA"];
      SPECEX_INFO("CAMERA = " << fitter.psf->camera_id );
    }else{
      SPECEX_WARNING("CAMERA Id not found in header");
    }
    
    fitter.priors = priors;
    
    
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
      }
      

      // now check mask ?
      if(psf->ParamsOfBundles[bundle].fiber_min < first_fiber) {
	psf->ParamsOfBundles[bundle].fiber_min = first_fiber;
	SPECEX_INFO("restricting fiber range first fiber = " << first_fiber);
      }
      if(psf->ParamsOfBundles[bundle].fiber_max > last_fiber) {
	psf->ParamsOfBundles[bundle].fiber_max = last_fiber;
	SPECEX_INFO("restricting fiber range last fiber = " << last_fiber);
      }
      
      
      
      // load traces in PSF 
      // --------------------------------------------
      for(int fiber=psf->ParamsOfBundles[bundle].fiber_min; fiber<=psf->ParamsOfBundles[bundle].fiber_max; fiber++) {
	if(traceset.find(fiber) != traceset.end())
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

      
      int margin = -psf->hSizeY+1; // we need to include spots that contribute to the image signal
      ymin+=margin;
      ymax-=margin;
      

      SPECEX_INFO("valid y(=rows) range = " << ymin << " " << ymax);

      vector<Spot_p> spots;
      allocate_spots_of_bundle(spots,lamp_lines_filename,psf->FiberTraces,bundle,psf->ParamsOfBundles[bundle].fiber_min,psf->ParamsOfBundles[bundle].fiber_max,ymin,ymax,min_wavelength,max_wavelength);
      SPECEX_INFO("number of spots = " << spots.size());
      
      if(write_tmp_results)
	write_spots_xml(spots,"spots-init.xml");

      /* DEBUGGING
      {
	for(int fiber=psf->ParamsOfBundles[bundle].fiber_min; fiber <= psf->ParamsOfBundles[bundle].fiber_max ;fiber ++) {
	  cout << "DEBUGGING check trace in fiber " << fiber << endl;
	  cout<< "X_vs_W xmin,xmax,deg " << psf->GetTrace(fiber).X_vs_W.xmin << " " << psf->GetTrace(fiber).X_vs_W.xmax << " " << psf->GetTrace(fiber).X_vs_W.deg << endl;
	  cout<< "Y_vs_W xmin,xmax,deg " << psf->GetTrace(fiber).Y_vs_W.xmin << " " << psf->GetTrace(fiber).Y_vs_W.xmax << " " << psf->GetTrace(fiber).Y_vs_W.deg << endl;
	  cout<< "X_vs_Y xmin,xmax,deg " << psf->GetTrace(fiber).X_vs_Y.xmin << " " << psf->GetTrace(fiber).X_vs_Y.xmax << " " << psf->GetTrace(fiber).X_vs_Y.deg << endl;
	  cout<< "W_vs_Y xmin,xmax,deg " << psf->GetTrace(fiber).W_vs_Y.xmin << " " << psf->GetTrace(fiber).W_vs_Y.xmax << " " << psf->GetTrace(fiber).W_vs_Y.deg << endl;
	}
	exit(12);
      }
      */
      

      // starting fit
      // --------------------------------------------
      bool init_psf = (input_psf_filename=="");
      fitter.FitEverything(spots,init_psf);

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
      write_psf_fits(fitter.psf,output_fits_filename,&fitted_spots);
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


#include <specex_serialization_implement.h>


