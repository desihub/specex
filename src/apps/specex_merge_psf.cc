#include <iostream>
#include <fstream>
#include <ctime>

#include <boost/program_options.hpp>

#include "harp.hpp"

#include "specex_message.h"
#include "specex_psf.h"
#include "specex_psf_io.h"
#include "specex_serialization.h"

using namespace std;
namespace popts = boost::program_options;

void usage(const char* pg) {

  cerr << pg << " file1.xml file2.xml ... --out-fits ofile.fits --out-xml ofile.xml" << endl;
  
  exit(12);
}


int main( int argc, char *argv[] ) {
  
  vector<string> input_xml_filenames;
  string output_xml_filename="";
  string output_fits_filename="";
  
  specex_set_verbose(true);

  if(argc<3) usage(argv[0]);
  
  
  for(int a=1;a<argc;a++) {

    if(strcmp(argv[a],"--out-fits")==0) {
      output_fits_filename = argv[++a];
      continue;
    }
    if(strcmp(argv[a],"--out-xml")==0) {
      output_xml_filename = argv[++a];
      continue;
    }
    
    if(argv[a][0]=='-') usage(argv[0]);
    
    input_xml_filenames.push_back(argv[a]);
    
  }
  if( output_fits_filename == "" && output_xml_filename == "" ) usage(argv[0]);
  if( input_xml_filenames.size() == 0) usage(argv[0]);

  specex::PSF_p psf;
  specex::read_psf_xml(psf,input_xml_filenames[0]);
  
  for(size_t i=1;i<input_xml_filenames.size();i++) {
    
    specex::PSF_p psf_i;
    specex::read_psf_xml(psf_i,input_xml_filenames[i]);
    
    psf->Append(psf_i);
  }
  
  if(output_fits_filename!="")
    specex::write_psf_fits(psf,output_fits_filename);
  
  if(output_xml_filename!="")
    specex::write_psf_xml(psf,output_xml_filename);
  
  return EXIT_SUCCESS;
}
