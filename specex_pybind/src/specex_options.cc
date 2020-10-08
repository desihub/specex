#include <string>
#include <iostream>
#include <specex_options.h>
#include <boost/program_options.hpp>

using namespace std;

namespace popts = boost::program_options;

int specex::Options::parse(int argc, char *argv[] ){

  if(getenv("SPECEXDATA"))
    lamp_lines_filename =
      string(getenv("SPECEXDATA"))+"/specex_linelist_desi.txt";

  // reading arguments
  // --------------------------------------------

  desc.add_options()
    ( "help,h", "display usage information" )
    ( "arc,a", popts::value<string>( &arc_image_filename ), "arc preprocessed fits image file name (mandatory)" )
    ( "in-psf", popts::value<string>( &input_psf_filename ), " input psf file name (fits or xml, can contain only the traces, mandatory)")  
    ( "out-psf", popts::value<string>( &output_fits_filename ), " output psf fits file name (mandatory)")  
    ( "only-trace-fit", "only fit the trace coordinates (default is trace,sigma,psf)")
    ( "no-trace-fit", "do not fit the trace coordinates (default is trace,sigma,psf)")
    ( "no-sigma-fit", "do not fit the gaussian sigma (default is trace,sigma,psf)")
    ( "no-psf-fit", "do not fit the psf (default is trace,sigma,psf)")
    ( "flux-hdu", popts::value<int>( &flux_hdu )->default_value(1), " flux hdu in input arc fits, for unusual data format")
    ( "ivar-hdu", popts::value<int>( &ivar_hdu )->default_value(2), " ivar hdu in input arc fits, for unusual data format")
    ( "mask-hdu", popts::value<int>( &mask_hdu )->default_value(3), " mask hdu in input arc fits, for unusual data format")
    ( "header-hdu", popts::value<int>( &header_hdu )->default_value(1), " header hdu in input arc fits, for unusual data format")

    // ( "trace", popts::value<string>( &trace_filename ), "fits image file name with image extensions XTRACE and YTRACE for legendre polynomial of wavelength (mandatory)" )
    
    ( "xcoord-hdu", popts::value<int>( &xtrace_hdu )->default_value(-1), "hdu of x trace legendre polynomial of wavelength (default is extension XTRACE)" )
    ( "ycoord-hdu", popts::value<int>( &ytrace_hdu )->default_value(-1), "hdu of y trace legendre polynomial of wavelength (default is extension YTRACE)" )
    ( "first-bundle", popts::value<int>( &first_fiber_bundle )->default_value(0), "first fiber bundle to fit")
    ( "last-bundle", popts::value<int>( &last_fiber_bundle )->default_value(0), "last fiber bundle to fit")
    ( "single-bundle", "fit a single bundle with all fibers")
    ( "first-fiber", popts::value<int>( &first_fiber )->default_value(0), "first fiber (must be in bundle)")
    ( "last-fiber", popts::value<int>( &last_fiber )->default_value(100000), "last fiber (must be in bundle)")
    ( "half-size-x", popts::value<int>( &half_size_x )->default_value(8), "half size of PSF stamp (full size is 2*half_size+1)")
    ( "half-size-y", popts::value<int>( &half_size_y )->default_value(5), "half size of PSF stamp (full size is 2*half_size+1)")
    ( "psfmodel", popts::value<string>( &psf_model )->default_value("GAUSSHERMITE"), "PSF model, default is GAUSSHERMITE")
    ( "positions", "fit positions of each spot individually after global fit for debugging")
    ( "verbose,v", "turn on verbose mode (deprecated, true by default)" )
    ( "quiet", "no info message, only warning" )
    ( "debug", "turn on debug mode" )
    ( "trace-prior-deg", popts::value<int>( &trace_prior_deg )->default_value(0) , "force equal trace coeff in bundle starting at this degree")
    ( "lamplines", popts::value<string>( &lamp_lines_filename ), "lamp lines ASCII file name (def. is $SPECEXDATA/specex_linelist_desi.txt)" )
    ( "core", "dump core files when harp exception is thrown" )
    ( "gauss-hermite-deg",  popts::value<int>( &gauss_hermite_deg )->default_value(6), "degree of Hermite polynomials (same for x and y, only if GAUSSHERMITE psf)")
    ("gauss-hermite-deg2",  popts::value<int>( &gauss_hermite_deg2 )->default_value(2), "degree of Hermite polynomials (same for x and y, only if GAUSSHERMITE2 psf)")
    //( "gauss_hermite_sigma",  popts::value<double>( &gauss_hermite_sigma ), "sigma of Gauss-Hermite PSF (same for x and y, only if GAUSSHERMITE psf)")
    ( "legendre-deg-wave",  popts::value<int>( &legendre_deg_wave )->default_value(3), "degree of Legendre polynomials along wavelength (can be reduced if missing data)")
    ( "legendre-deg-x",  popts::value<int>( &legendre_deg_x )->default_value(1), "degree of Legendre polynomials along x_ccd (can be reduced if missing data)")
    ( "trace-deg-wave",  popts::value<int>( &trace_deg_wave )->default_value(6), "degree of Legendre polynomials along wavelength for fit of traces")
    ( "trace-deg-x",  popts::value<int>( &trace_deg_x )->default_value(6), "degree of Legendre polynomials along x_ccd for fit of traces")
    ( "psf-error",  popts::value<double>( &psf_error )->default_value(0), "psf fractional uncertainty (default is 0.01, for weights in the fit)")
    ( "psf-core-wscale",  popts::value<double>( &psf_core_wscale ), "scale up the weight of pixels in 5x5 PSF core")
    ( "broken-fibers", popts::value<string>( &broken_fibers_string ), "broken fibers (comma separated list)")
    
#ifdef EXTERNAL_TAIL
    ( "fit-psf-tails", "unable fit of psf tails")
#endif
#ifdef CONTINUUM
    ( "fit-continuum", "unable fit of continuum")
#endif
    ( "variance-model", "refit at the end with a model of the variance to avoid Poisson noise bias")
    
    ( "out-psf-xml", popts::value<string>( &output_xml_filename ), " output psf xml file name")
   
    ( "out-spots", popts::value<string>( &output_spots_filename ), " output spots file name")  
    ( "prior", popts::value< vector<string> >( &argurment_priors )->multitoken(), " gaussian prior on a param : 'name' value error")  
    ( "tmp_results", " write tmp results")  
    ( "nlines", popts::value<int>( &max_number_of_lines )->default_value(200)," maximum number of emission lines used in fit (uses a algorithm to select best ones based on S/N and line coverage")  
    //( "out", popts::value<string>( &outfile ), "output image file" )
    ;

  try {    
    popts::store(popts::command_line_parser( argc, argv ).options(desc).run(), vm);
    popts::notify(vm);
  }catch(std::exception& e) {
    cerr << "error in arguments : " << e.what() << endl;
    cerr << "try --help for options" << endl;
    return EXIT_FAILURE;
  }
  
  if ( ( argc < 2 ) || vm.count( "help" ) ) {
      cerr << endl;
      cerr << desc << endl;      
      return EXIT_FAILURE;
  }
  if ( ! vm.count( "arc" ) ) {
    cerr << endl;
    cerr << "missing --arc , try --help for options" << endl;      
    return EXIT_FAILURE;
  }
  if ( ! vm.count( "in-psf" ) ) {
    cerr << endl;
    cerr << "missing --in-psf , try --help for options" << endl;      
    return EXIT_FAILURE;
  }    
  if(lamp_lines_filename == "") {
    cerr << endl;
    cerr << "missing lamp_lines_filename either define env. variable SPECEXDATA or use option --lamplines" << endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
  
}

