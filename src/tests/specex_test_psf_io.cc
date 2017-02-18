#include <iostream>
#include <fstream>
#include <ctime>

#include <harp.hpp>

#include "specex_message.h"
#include "specex_psf.h"
#include "specex_psf_io.h"
#include "specex_serialization.h"
#include <boost/program_options.hpp>

namespace popts = boost::program_options;
using namespace std;

int main(int argc, char *argv[]) {
  
  
  string filename="";
  popts::options_description desc ( "Allowed Options" );
  desc.add_options()
    ( "help,h", "display usage information" )
    ( "in,i", popts::value<string>( &filename ), "input psf filename" );

  popts::variables_map vm;
  try {
    popts::store(popts::command_line_parser( argc, argv ).options(desc).run(), vm);
    popts::notify(vm);
  }catch(std::exception e) {
    cerr << "error in arguments" << endl;
    return EXIT_FAILURE;
  }
  if(filename=="") {
    cerr << "need input filename" << endl;
    return EXIT_FAILURE;
  }
  specex_set_verbose(true);
  specex::PSF_p psf;
  specex::read_psf(psf,filename);
  specex::write_psf_xml(psf,"toto.xml");
  specex::write_psf_fits(psf,"toto.fits");
  
  return EXIT_SUCCESS;
}
