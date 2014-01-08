#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <map>

#include <boost/program_options.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>

#include <harp.hpp>

#include <specex_spot.h>
#include <specex_spot_array.h>
#include <specex_message.h>

using namespace std;
using namespace specex;

namespace popts = boost::program_options;


int main ( int argc, char *argv[] ) {
  
  string spots_xml_filename="";
  string spots_list_filename="";
  
  
  // reading arguments
  // --------------------------------------------
  popts::options_description desc ( "Allowed Options" );
  desc.add_options()
    ( "help,h", "display usage information" )
    ( "in", popts::value<string>( &spots_xml_filename ), "" )
    ( "out", popts::value<string>( &spots_list_filename ), "")
    ;

  popts::variables_map vm;
  popts::store(popts::command_line_parser( argc, argv ).options(desc).run(), vm);
  popts::notify(vm);
  
  if ( ( argc < 2 ) || vm.count( "help" ) || ( ! vm.count( "in" ) ) || ( ! vm.count( "out" ) ) ) {
    cerr << endl;
    cerr << desc << endl;
    cerr << "example:" << endl;
    cerr << argv[0] << " --in spots.xml --out spots.list" << endl;
    return EXIT_FAILURE;
  }
  
  
  try {

    
    // read spots
    // --------------------------------------------
    vector<specex::Spot_p> spots;
    {
      std::ifstream is(spots_xml_filename.c_str());
      boost::archive::xml_iarchive xml_ia ( is );
      xml_ia >> BOOST_SERIALIZATION_NVP(spots);
      is.close();
    } 
    
    SPECEX_INFO("number of spots = " << spots.size());
    
    write_spots_list(spots,spots_list_filename);
    

  // ending
  // --------------------------------------------
  }catch(harp::exception e) {
    cerr << "FATAL ERROR (harp) " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}



