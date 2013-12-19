#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <map>

#include <boost/program_options.hpp>

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
#include <specex_gauss_hermite_analytic_psf.h>
#include <specex_serialization.h>



using namespace std;
using namespace specex;

namespace popts = boost::program_options;




int main ( int argc, char *argv[] ) {
  

  

  // default arguments
  // --------------------------------------------
  string psf_model = "GAUSSHERMITE";
  string spectrograph_name = "BOSS";
  int    fiber_bundle=1;
  
  string arc_image_filename="";
  string xy_trace_fits_name="";
  string wy_trace_fits_name="";
  string lamp_lines_filename="";  
  double min_wavelength = 0;
  double max_wavelength = 1e6;
  
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
    ( "bundle,b", popts::value<int>( &fiber_bundle ), "fiber bundle")
    ( "psfmodel", popts::value<string>( &psf_model ), "PSF model, default is GAUSSHERMITE")
    ( "positions", "fit positions of each spot individually after global fit for debugging")
    ( "verbose,v", "turn on verbose mode" )
    ( "lamplines", popts::value<string>( &lamp_lines_filename ), "lamp lines ASCII file name (def. is $IDLSPEC2D_DIR/opfiles/lamplines.par)" )
    ( "core", "dump core files when harp exception is thrown" )
    //( "out", popts::value<string>( &outfile ), "output image file" )
    ;

  popts::variables_map vm;
  popts::store(popts::command_line_parser( argc, argv ).options(desc).run(), vm);
  popts::notify(vm);
  
  if ( ( argc < 2 ) || vm.count( "help" ) || ( ! vm.count( "arc" ) )  || ( ! vm.count( "xy" ) ) || ( ! vm.count( "wy" ) ) ) {
    cerr << endl;
    cerr << desc << endl;
    cerr << "example:" << endl;
    cerr << argv[0] << " --arc sdProc-b1-00108382.fits --xy redux/spFlat-b1-00108381.fits.gz --wy redux/spArc-b1-00108382.fits.gz -v" << endl;
    return EXIT_FAILURE;
  }
  
  specex_set_verbose(vm.count("verbose")>0);
  specex_set_dump_core(vm.count("core")>0);
  
  bool fit_individual_spots_position = vm.count("positions");
  
  SPECEX_INFO("using lamp lines file " << lamp_lines_filename); 

  
  try {
    
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
    int i=406; int j=787;
    cout << "DEBUG160 w = " << weight(i,j) << " " << i << " " << j << endl;
    // write_new_fits_image("weight.fits",weight); exit(12);
    
    // init PSF
    // --------------------------------------------
    specex::PSF_p psf;

    if(psf_model=="GAUSSHERMITE") {
      psf = PSF_p(new specex::GaussHermitePSF(4));
      //} else if(psf_model=="GAUSSIAN") {
      //psf = new GaussPSF();
    }else {
      SPECEX_ERROR("don't know this psf model");
    }
    
    
    // define spectrograph (this should be improved)
    // --------------------------------------------
    Spectrograph *spectro = 0;
    if(spectrograph_name == "BOSS") {
      spectro = new BOSS_Spectrograph();
    }else{
      SPECEX_ERROR("unknown spectrograph");
    }
    

    // read trace set derived in BOSS pipeline (this should be improved)
    // --------------------------------------------    
    TraceSet traceset;
    if(spectrograph_name == "BOSS") {
      int xy_trace_hdu = 1; // in the spFlat file , warning : this HDU numbering starts at 0 (to check)
      int wy_trace_hdu = 2; // int the spArc file , warning : this HDU numbering starts at 0 (to check)
      read_BOSS_traceset_in_fits(traceset,wy_trace_fits_name,wy_trace_hdu,xy_trace_fits_name,xy_trace_hdu);
    }
    
    // loading arc lamp spots
    // --------------------------------------------  
    int ymin = 696+3; // range of usable r CCD coordinates, hard coded for now
    int ymax = 3516-3; // range of usable r CCD coordinates, hard coded for now
    vector<Spot> spots;
    allocate_spots_of_bundle(spots,*spectro,lamp_lines_filename,traceset,fiber_bundle,ymin,ymax,min_wavelength,max_wavelength);
    SPECEX_INFO("number of spots = " << spots.size());
    
    write_spots_list(spots,"initial_spots.list");
    
    // load traces in PSF
    // --------------------------------------------
    psf->FiberTraces.clear();  
    for(int fiber=spectro->number_of_fibers_per_bundle*fiber_bundle; fiber<spectro->number_of_fibers_per_bundle*(fiber_bundle+1); fiber++) {
      psf->FiberTraces[fiber]=traceset[fiber];
    } 
    
    // init PSF fitter
    // -------------------------------------------- 
    PSF_Fitter fitter(psf,image,weight);
    
    fitter.gain = 1; // images are already in electrons
    // fitter.readout_noise = 2; // b1, evaluated using spatial variance in sdProc-b1-00108382.fits[1:4000,1:600]
    fitter.readout_noise = 2.5; // need to modify this image.KeyVal("RDNOISE0"); // WARNING use val from first amp here
    fitter.flatfield_error = 0.02; // 2%
    
    for(int j=0;j<weight.Ny();j++) {
      for(int i=0;i<weight.Nx();i++) {
	if(weight(i,j)>0)
	  weight(i,j)=1/square(fitter.readout_noise);
	else
	  weight(i,j)=0;
      }
    }

    fitter.mask.Clear();

    SPECEX_INFO(
		"PSF '" << psf_model << "' stamp size = " 
		<< psf->hSizeX << "x" << psf->hSizeY << " npar(loc) = " 
		<< psf->FixedCoordNPar() << " npar(glob) = " << psf->VaryingCoordNPar()
		);
    
    
    // starting fit
    // --------------------------------------------
    bool init_psf = true;
    fitter.verbose = true;
    std::vector<Spot*> spot_array = array_of_pointer(spots);
    fitter.FitEverything(spot_array,init_psf);
    
    write_psf_xml(fitter.psf,"psf.xml");
    write_psf_fits(fitter.psf,"psf.fits",500,2000,4);

    if(fit_individual_spots_position) // for debugging
      fitter.FitIndividualSpotPositions(spot_array);
    
    
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



