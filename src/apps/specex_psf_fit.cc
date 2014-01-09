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
#include <specex_gauss_hermite_psf.h>
#include <specex_serialization.h>



using namespace std;
using namespace specex;

namespace popts = boost::program_options;




int main ( int argc, char *argv[] ) {
  

  

  // default arguments
  // --------------------------------------------
  string psf_model = "GAUSSHERMITE";
  string spectrograph_name = "BOSS";
  int    first_fiber_bundle=1;
  int    last_fiber_bundle=1;
  int    first_fiber=0;
  int    last_fiber=100000;
  int    half_size=4;
  
  string arc_image_filename="";
  string xy_trace_fits_name="";
  string wy_trace_fits_name="";
  string lamp_lines_filename="";  
  double min_wavelength = 0;
  double max_wavelength = 1e6;
  
  int gauss_hermite_deg = 3;
  double gauss_hermite_sigma = 1.1;
  int legendre_deg_wave = 4;
  int legendre_deg_x = 1;
  
  double psf_error = 0.01;
  
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
    ( "half_size", popts::value<int>( &half_size ), "half size of PSF stamp (full size is 2*half_size+1)")
    ( "psfmodel", popts::value<string>( &psf_model ), "PSF model, default is GAUSSHERMITE")
    ( "positions", "fit positions of each spot individually after global fit for debugging")
    ( "verbose,v", "turn on verbose mode" )
    ( "lamplines", popts::value<string>( &lamp_lines_filename ), "lamp lines ASCII file name (def. is $IDLSPEC2D_DIR/opfiles/lamplines.par)" )
    ( "core", "dump core files when harp exception is thrown" )
    ( "gauss_hermite_deg",  popts::value<int>( &gauss_hermite_deg ), "degree of Hermite polynomials (same for x and y, only if GAUSSHERMITE psf)")
    ( "gauss_hermite_sigma",  popts::value<double>( &gauss_hermite_sigma ), "sigma of Gauss-Hermite PSF (same for x and y, only if GAUSSHERMITE psf)")
    ( "legendre_deg_wave",  popts::value<int>( &legendre_deg_wave ), "degree of Legendre polynomials along wavelength (can be reduced if missing data)")
    ( "legendre_deg_x",  popts::value<int>( &legendre_deg_x ), "degree of Legendre polynomials along x_ccd (can be reduced if missing data)")
    ( "psf_error",  popts::value<double>( &psf_error ), "psf fractional uncertainty (default is 0.01, for weights in the fit)")
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
    // weight is 0 or 1
    for(int j=0;j<weight.Ny();j++) {
      for(int i=0;i<weight.Nx();i++) {
	if(weight(i,j)>0)
	  weight(i,j)=1;
	else
	  weight(i,j)=0;
      }
    }
    
    // init PSF
    // --------------------------------------------
    specex::PSF_p psf;

    if(psf_model=="GAUSSHERMITE") {
      psf = PSF_p(new specex::GaussHermitePSF(gauss_hermite_deg));
      boost::static_pointer_cast<specex::GaussHermitePSF>(psf)->sigma = gauss_hermite_sigma;
      
      //} else if(psf_model=="GAUSSIAN") {
      //psf = new GaussPSF();
    }else {
      SPECEX_ERROR("don't know this psf model");
    }
    psf->ccd_image_n_cols = image.n_cols();
    psf->ccd_image_n_rows = image.n_rows();
    psf->hSizeX = half_size;
    psf->hSizeY = half_size;
    
        
    

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
    PSF_Fitter fitter(psf,image,weight);
    
    fitter.polynomial_degree_along_x    = legendre_deg_x;
    fitter.polynomial_degree_along_wave = legendre_deg_wave;
    fitter.psf_error                    = psf_error;

    fitter.gain = 1; // images are already in electrons
    // fitter.readout_noise = 2; // b1, evaluated using spatial variance in sdProc-b1-00108382.fits[1:4000,1:600]

    // compute mean readout noise
    fitter.readout_noise = 2;
    if(image_infos.find("RDNOISE0") != image_infos.end()) fitter.readout_noise = atof(image_infos["RDNOISE0"].c_str());
    
    
    for(int j=0;j<weight.Ny();j++) {
      for(int i=0;i<weight.Nx();i++) {
	if(weight(i,j)>0)
	  weight(i,j)=1/square(fitter.readout_noise);
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
    
    
    // loop on fiber bundles 
    // -------------------------------------------- 
    for(int bundle = first_fiber_bundle; bundle <= last_fiber_bundle ; bundle ++) {
    
      // allocate bundle in PSF if necessary
      if(psf->ParamsOfBundles.find(bundle)==psf->ParamsOfBundles.end()) {
	psf->ParamsOfBundles[bundle] = specex::PSF_Params();
	psf->ParamsOfBundles[bundle].fiber_min = spectro->number_of_fibers_per_bundle*bundle;
	psf->ParamsOfBundles[bundle].fiber_max = psf->ParamsOfBundles[bundle].fiber_min+spectro->number_of_fibers_per_bundle-1; // included
	
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
      int ymin = 696+3; // range of usable r CCD coordinates, hard coded for now
      int ymax = 3516-3; // range of usable r CCD coordinates, hard coded for now
      vector<Spot_p> spots;
      allocate_spots_of_bundle(spots,*spectro,lamp_lines_filename,traceset,bundle,psf->ParamsOfBundles[bundle].fiber_min,psf->ParamsOfBundles[bundle].fiber_max,ymin,ymax,min_wavelength,max_wavelength);
      SPECEX_INFO("number of spots = " << spots.size());
      
      // starting fit
      // --------------------------------------------
      fitter.FitEverything(spots,true);

      {
	// writing spots as xml
	std::ofstream os("spots.xml");
	boost::archive::xml_oarchive xml_oa ( os );
	xml_oa << BOOST_SERIALIZATION_NVP(spots);
	os.close();
	SPECEX_INFO("wrote spots in " << "spots.xml");
      }

      if(fit_individual_spots_position) // for debugging
	fitter.FitIndividualSpotPositions(spots);
    } // end of loop on bundles
    
    write_psf_xml(fitter.psf,"psf.xml");
    //write_psf_fits_image(fitter.psf,"psf.fits",500,2000,1,4);


    
    
    
    
    
  // ending
  // --------------------------------------------
  }catch(harp::exception e) {
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



