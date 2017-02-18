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
  
  
  string input_filename="";
  string output_filename="";
  popts::options_description desc ( "Allowed Options" );
  desc.add_options()
    ( "help,h", "display usage information" )
    ( "in,i", popts::value<string>( &input_filename ), "input psf filename (xml or fits)" )
  ( "out,o", popts::value<string>( &output_filename ), "output psf filename (xml or fits)" );

  popts::variables_map vm;
  try {
    popts::store(popts::command_line_parser( argc, argv ).options(desc).run(), vm);
    popts::notify(vm);
  }catch(std::exception e) {
    cerr << "error in arguments" << endl;
    return EXIT_FAILURE;
  }
  if(input_filename=="") {
    cerr << "need input filename" << endl;
    return EXIT_FAILURE;
  }
  if(output_filename=="") {
    cerr << "need input filename" << endl;
    return EXIT_FAILURE;
  }
  
  specex_set_verbose(true);
  specex::PSF_p psf;
  specex::read_psf(psf,input_filename);
  specex::write_psf(psf,output_filename);
  
  
  return EXIT_SUCCESS;
}
