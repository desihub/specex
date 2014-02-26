#include <harp_plugin_specex.h>
#include <specex_psf_io.h>

using namespace std;


harp::specex_psf::specex_psf ( boost::property_tree::ptree const & props ) : harp::psf ( "specex", props ) {

  // read specex psf in xml format
  
  cout << "Reading a specex::PSF" << endl;
  
  string file_type = props.get("type","");
  string path = props.get("path","");
  
  
  if( file_type == "xml")
    read_psf_xml(actual_specex_psf,path);
  else if( file_type == "fits")
    read_psf_fits(actual_specex_psf,path);
  
  nspec_ = actual_specex_psf->FiberTraces.size(); // assume this is the number of fibers
  nlambda_ = 0; // I don't have lambda
  rows_ = actual_specex_psf->ccd_image_n_rows;
  cols_ = actual_specex_psf->ccd_image_n_cols;

  
  cout << "This PSF is a " << actual_specex_psf->Name() << endl;
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

