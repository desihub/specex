#include <iostream>
#include <fstream>
#include <ctime>

#include <boost/program_options.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>

#include <harp.hpp>

#include "specex_message.h"
#include "specex_fits.h"
#include "specex_psf.h"
#include "specex_psf_io.h"
#include "specex_gauss_hermite_psf.h"
#include "specex_image_data.h"
#include "specex_serialization.h"

using namespace std;
using namespace specex;

namespace popts = boost::program_options;

int main(int argc, char *argv[]) {


  string psf_filename = "";
  string output_fits_image_filename = "";
  int fiber = 0;
  double wavelength = 6000;
  int oversampling = 1;
  
  int    half_size_x=0;
  int    half_size_y=0;
  
  // reading arguments
  // --------------------------------------------
  popts::options_description desc ( "Allowed Options" );
  desc.add_options()
    ( "help,h", "display usage information" )
    ( "psf", popts::value<string>( &psf_filename ), "psf xml or fits filename" )
    ( "out", popts::value<string>( &output_fits_image_filename ), " output fits image file name")
    ( "fiber", popts::value<int>( &fiber ), "")
    ( "wavelength", popts::value<double>( &wavelength ), "")
    ( "oversampling", popts::value<int>( &oversampling ), "")
    ( "half_size_x", popts::value<int>( &half_size_x ), "half size of PSF stamp (full size is 2*half_size+1)")
    ( "half_size_y", popts::value<int>( &half_size_y ), "half size of PSF stamp (full size is 2*half_size+1)")
    ( "debug", "debug mode" )
    ;
  
  popts::variables_map vm;
  
  try {
  
    
    popts::store(popts::command_line_parser( argc, argv ).options(desc).run(), vm);
    popts::notify(vm);
    
    if ( ( argc < 2 ) || vm.count( "help" ) || ( ! vm.count( "psf" ) )  || ( ! vm.count( "out" ) ) 
	 || ( ! vm.count( "fiber" ) )  || ( ! vm.count( "wavelength" ) ) 
	 ) {
      cerr << endl;
      cerr << desc << endl;
     return EXIT_FAILURE;
    }
  }catch(std::exception) {
    cerr << "error in arguments" << endl;
    cerr << endl;
    cerr << desc << endl;
    return EXIT_FAILURE;
  }
  try {

    specex_set_verbose(true);
    specex_set_debug(vm.count("debug")>0);
    
    specex::PSF_p psf;
    specex::read_psf(psf,psf_filename);
    
  
    if(half_size_x>0) psf->hSizeX = half_size_x;
    if(half_size_y>0) psf->hSizeY = half_size_y;
    
    cout << "wavelength = " << wavelength << " fiber = " << fiber << endl;
    write_psf_fits_image(psf,output_fits_image_filename,fiber,wavelength,oversampling);
    
  }catch(harp::exception e) {
    cerr << "FATAL ERROR (harp) " << e.what() << endl;
    return EXIT_FAILURE;
  }catch(std::exception) {
    cerr << "unknown error" << endl;
    return EXIT_FAILURE;
  }
  
  

  return EXIT_SUCCESS;
}
