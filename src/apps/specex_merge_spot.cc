#include <iostream>
#include <fstream>
#include <ctime>

#include <boost/program_options.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>

#include "harp.hpp"

#include "specex_message.h"
#include "specex_spot.h"
#include <specex_spot_array.h>
#include "specex_serialization.h"

using namespace std;
namespace popts = boost::program_options;

void usage(const char* pg) {

  cerr << pg << " file1.xml file2.xml ... --out ofile.xml" << endl;
  
  exit(12);
}


int main( int argc, char *argv[] ) {
  
  vector<string> input_xml_filenames;
  string output_xml_filename="";
  
  specex_set_verbose(true);

  if(argc<3) usage(argv[0]);
  
  
  for(int a=1;a<argc;a++) {

    if(strcmp(argv[a],"--out")==0) {
      output_xml_filename = argv[++a];
      continue;
    }
    
    if(argv[a][0]=='-') usage(argv[0]);
    
    input_xml_filenames.push_back(argv[a]);
    
  }
  if( output_xml_filename == "" ) usage(argv[0]);
  if( input_xml_filenames.size() == 0) usage(argv[0]);

  vector<specex::Spot_p> spots;
  for(size_t i=0;i<input_xml_filenames.size();i++) {
    
    vector<specex::Spot_p> spots_i;
    std::ifstream is(input_xml_filenames[i].c_str());
    boost::archive::xml_iarchive xml_ia ( is );
    xml_ia >> BOOST_SERIALIZATION_NVP(spots_i);
    is.close();

    for(size_t s=0;s<spots_i.size();s++)
      spots.push_back(spots_i[s]);
  }

  {
    
    std::ofstream os(output_xml_filename.c_str());
    boost::archive::xml_oarchive xml_oa ( os );
    xml_oa << BOOST_SERIALIZATION_NVP(spots);
    os.close();
    SPECEX_INFO("wrote spots in " << output_xml_filename);
  }
  
  return EXIT_SUCCESS;
}
