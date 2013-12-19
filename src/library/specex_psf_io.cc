#include <fstream>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>

#include <harp.hpp>

#include <specex_psf.h>
#include <specex_psf_io.h>
#include <specex_message.h>
#include <specex_image_data.h>
#include <specex_fits.h>

// also need #included <specex_serialisation.h>  before main() 

void specex::write_psf_fits_image(const specex::PSF_p psf, const string& filename, const double& x, const double& y, int oversampling) {
  
  harp::vector_double P=psf->FixedCoordParams(x,y);
  
  int nx = 2*psf->hSizeX*oversampling+1;
  int ny = 2*psf->hSizeY*oversampling+1;
  
  

  specex::image_data img(nx,ny);
  for(int j=0;j<ny;j++) {
    for(int i=0;i<nx;i++) {
      
      int ib = (i-nx/2)/oversampling;
      int jb = (j-ny/2)/oversampling;
      double dx = (i-nx/2)/double(oversampling)-ib;
      double dy = (j-ny/2)/double(oversampling)-jb;
      
      img(i,j)=psf->PSFValueWithParams(x-dx,y-dy,ib+int(x),jb+int(y),P,0,0);
    }
  }
  specex::write_new_fits_image(filename,img);
}



void specex::write_psf_xml(const specex::PSF_p psf, const std::string& filename) {
  
  std::ofstream os(filename.c_str());
  boost::archive::xml_oarchive xml_oa ( os );

  xml_oa << BOOST_SERIALIZATION_NVP(psf);
  
  os.close();

  SPECEX_INFO("wrote psf in " << filename);
}

void specex::read_psf_xml(specex::PSF_p& psf, const std::string& filename) {
  
  
  std::ifstream is(filename.c_str());
  
  boost::archive::xml_iarchive xml_ia ( is );

  xml_ia >> BOOST_SERIALIZATION_NVP(psf);
  
  is.close();

  SPECEX_INFO("read psf in " << filename);
  
}
