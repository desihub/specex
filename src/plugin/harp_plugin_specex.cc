#include <harp_plugin_specex.h>
#include <specex_psf.h>
#include <specex_psf_io.h>
#include <specex_serialization.h>

using namespace std;

// HARP/src/libharp/data/harp/psf.hpp

harp::specex_psf::specex_psf ( boost::property_tree::ptree const & props ) : harp::psf ( "specex", props ) {

  // read specex psf in xml format
  
  cout << "Reading a specex::PSF" << endl;
  
  string file_type = props.get("type","");
  string path = props.get("path","");
  double wavebin = props.get("wavebin",1);
  
  if(wavebin<=0) {
    HARP_THROW("unphysical wave bin");
  }

  if( file_type == "xml")
    read_psf_xml(actual_specex_psf,path);
  else if( file_type == "fits")
    read_psf_fits(actual_specex_psf,path);
  
  {
    // compute range of wavelength accessible for all fibers as defined by the measured arc lamp lines
    double wavemin=0;
    double wavemax=100000;
    for(std::map<int,specex::Trace>::const_iterator it = actual_specex_psf->FiberTraces.begin(); it != actual_specex_psf->FiberTraces.end(); ++it) {
      wavemin=max(wavemin,it->second.Y_vs_W.xmin);
      wavemax=min(wavemax,it->second.Y_vs_W.xmax);
    }
    
    int nlambda = int((wavemax-wavemin)/wavebin);
    lambda_.resize(nlambda);
    for(int i=0;i<nlambda;i++)
      lambda_(i)=wavemin+wavebin*i;
  }

  
  nspec_ = actual_specex_psf->FiberTraces.size(); // assume this is the number of fibers
  rows_  = actual_specex_psf->ccd_image_n_rows;
  cols_  = actual_specex_psf->ccd_image_n_cols;
  
  cout << "This PSF is a " << actual_specex_psf->Name() << endl;
}


void harp::specex_psf::response ( size_t spec_index, size_t lambda_index, size_t & x_offset, size_t & y_offset, harp::matrix_double & patch ) const {
  
  int nx=2*actual_specex_psf->hSizeX+1;
  int ny=2*actual_specex_psf->hSizeY+1;
  
  std::map<int,specex::Trace>::const_iterator it = actual_specex_psf->FiberTraces.find(int(spec_index));
  if(it==actual_specex_psf->FiberTraces.end()) HARP_THROW("specex_psf::response don't have PSF for this fiber/spec");
  const specex::Trace &trace = it->second;
  double x_center = trace.X_vs_W.Value(lambda_(lambda_index));
  double y_center = trace.Y_vs_W.Value(lambda_(lambda_index));
  int x_pic = int(floor(x_center));
  
  HARP_THROW("NEED TO FINISH THIS");


}


size_t harp::specex_psf::response_nnz_estimate ( ) const {
  HARP_THROW( "julien doesn't know what response_nnz_estimate means." );
}


// This export statement allows use of the plugin with serialization operations in HARP.

BOOST_CLASS_EXPORT(harp::specex_psf)

// Define the plugin creation function here.

harp::psf * harp::specex_psf_create ( boost::property_tree::ptree const & props ) {
  return new harp::specex_psf ( props );
}


// Define the initialize function that will be called by the plugin registry.  There must
// be only one such function per *.so file, but you can register multiple plugins within
// this function call.

void initialize ( void * registry ) {

  harp::plugin_registry * reg = static_cast < harp::plugin_registry * > ( registry );

  string const & version = harp::source_version();

  // register psf plugin
  reg->register_psf ( "specex", harp::specex_psf_create, version );
  
  return;
}

