#include <iostream>
#include <fstream>

#include "harp.hpp"

#include "specex_message.h"
#include "specex_psf.h"
#include "specex_psf_io.h"
#include "specex_gauss_hermite_analytic_psf.h"
#include "specex_serialization.h"

using namespace std;

int main() {
  
  specex_set_verbose(true);
  
  
  if(0){
    specex::PSF_p psf(new specex::GaussHermitePSF());
    std::ofstream os("toto.xml");
    boost::archive::xml_oarchive xml_oa ( os );
    xml_oa << BOOST_SERIALIZATION_NVP(psf);
    os.close();
  }
  
  if(0){
    std::ifstream is("toto.xml");
    boost::archive::xml_iarchive xml_ia ( is );
    specex::PSF_p psf;
    xml_ia >> BOOST_SERIALIZATION_NVP(psf);
    is.close();    
  }
  if(0){
    specex::PSF_p psf(new specex::GaussHermitePSF());
    specex::write_psf_xml(psf,"toto.xml");
  }
  if(0){
    specex::PSF_p psf;
    specex::read_psf_xml(psf,"toto.xml");
  }
  if(1) {
    specex::PSF_p psf;
    specex::read_psf_xml(psf,"psf.xml");
    specex::write_psf_xml(psf,"psf2.xml");

    specex::write_psf_fits(psf,"psf2.fits",500,2000,4);
    
  }
  return EXIT_SUCCESS;
}
