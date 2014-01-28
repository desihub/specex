#include <iostream>
#include <cstdio>

#include <boost/program_options.hpp>

#include <harp.hpp>

#include <specex_fits.h>
#include <specex_image_data.h>
#include <specex_message.h>
#include <specex_linalg.h>


using namespace std;
//using namespace harp;
namespace popts = boost::program_options;


int main ( int argc, char *argv[] ) {

  specex_set_verbose(true);

  double tstart;
  double tstop;

  cout.precision ( 12 );
  cerr.precision ( 12 );

  string infile = "";
  string outfile = "toto.fits";
  
  
  // Parse commandline options
  
  popts::options_description desc ( "Allowed Options" );
  
  desc.add_options()
  ( "help,h", "display usage information" )
  ( "in", popts::value<string>( &infile ), "input image file" )
  ( "out", popts::value<string>( &outfile ), "output image file" )
  ;

  popts::variables_map vm;

  popts::store(popts::command_line_parser( argc, argv ).options(desc).run(), vm);
  
  popts::notify(vm);

  if ( ( argc < 2 ) || vm.count( "help" ) || ( ! vm.count( "in" ) ) ) {
    cerr << endl;
    cerr << desc << endl;
    return 0;
  }
  
  // reading input
  
  boost::property_tree::ptree props;
  props.put ( "format", "fits" );
  props.put ( "path", infile );
  
  // ./src/libharp-mpi/harp/image.hpp:  typedef boost::shared_ptr < harp::image > image_p;
  harp::image_p img ( harp::image::create ( props ) );
  size_t nr = img->n_rows();
  size_t nc = img->n_cols();
  
  // ./src/libharp/math/harp/linalg.hpp:  typedef boost::numeric::ublas::vector < double > vector_double;
  harp::vector_double vals(nc*nr);
  img->values(vals);
  
  vals *= 2;
  

  // now write image
  specex::write_new_fits_image(outfile,nr,nc,vals);
  
  specex::image_data img2(12,20); specex::zero(img2.data);
  for(int i=0;i<img2.Nx();i++) {
    img2(i,3)=1;
  }
  specex::write_new_fits_image("tata.fits",img2);
  

  return 0;
}
