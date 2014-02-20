#include <iostream>
#include <fstream>
#include <ctime>

#include <boost/program_options.hpp>

#include "harp.hpp"

#include "specex_message.h"
#include "specex_fits.h"
#include "specex_psf.h"
#include "specex_psf_io.h"
#include "specex_gauss_hermite_psf.h"
#include "specex_image_data.h"
#include "specex_serialization.h"
#include "specex_model_image.h"

using namespace std;
namespace popts = boost::program_options;

int main( int argc, char *argv[] ) {
  
  string input_psf_xml_filename="";
  string input_psf_fits_filename="";
  string output_image_fits_filename="";
  
  specex_set_verbose(true);
  
  popts::options_description desc ( "Allowed Options" );
  desc.add_options()
    ( "help,h", "display usage information" )
    ( "in_xml", popts::value<string>( & input_psf_xml_filename ), "input PSF xml file name" )
    ( "in_fits", popts::value<string>( & input_psf_fits_filename ), "input PSF fits file name" )
    ( "out", popts::value<string>( & output_image_fits_filename ), "output image fits file name" )
    ;
  
  popts::variables_map vm;
  try {
    popts::store(popts::command_line_parser( argc, argv ).options(desc).run(), vm);
    popts::notify(vm);
    
    if ( ( argc < 2 ) || vm.count( "help" ) || ( ! vm.count( "in_xml" ) && ( ! vm.count( "in_fits" ) )  || ( ! vm.count( "out" ) )  ) ) {
      cerr << endl;
      cerr << desc << endl;
      return EXIT_FAILURE;
    }
  }catch(std::exception e) {
    cerr << "error in arguments" << endl;
    cerr << endl;
    cerr << desc << endl;
    return EXIT_FAILURE;
  }
  
  specex::PSF_p psf;
  
  if(vm.count( "in_fits" ))
    SPECEX_ERROR("reading fits not implemented yet");
  
  specex::read_psf_xml(psf,input_psf_xml_filename);
  
  //specex::write_psf_fits(psf,"toto.fits");
  

  if(0){
    // testing coordinate system  
    double xc = 400;
    double yc = 2000;
    harp::vector_double P(psf->LocalNAllPar());
    P.clear();
    int i=floor(xc);
    //for(i=floor(xc)-3;i<=floor(xc)+3;i++)
    for(int j=floor(yc)-3;j<=floor(yc)+3;j++) 
      cout << i << " " << j << " " << psf->PSFValueWithParamsXY(400,2000,i,j,P,NULL,NULL) << endl;
    
    return 0;
  }
  
  /* 
  // write an image of the psf at a precise location
  write_psf_fits_image(psf,"psf-image.fits",400,2000,1,4);
  
  {
    specex::PSF_p ipsf(new specex::GaussHermitePSF());
    // read psf in specter format
    
    read_psf_fits(ipsf,"test-psf-for-specter2.fits");
    
    write_psf_fits_image(ipsf,"psf-image-bis.fits",400,2000,1,4);
  }

  */
  

  // now create a test image
  SPECEX_INFO("Generate an image with spots of flux=1 at each 50A starting at 3500A for each fiber addressed by PSF");
  
  /* dump things
  
  For fiber 20 at wavelength=4000 Angstroms:
  
  What is the CCD x,y position?
  What is the reduced x', y' in the ranges [-1,1]?
  What are the Legendre coefficients that contribute to GH(0,0)?
  After summing those multiplied by Legendre polynomials evaluated at x' and y', what is the final coefficient for the GH(0,0) term?
  
  I'll compare to those numbers and then we'll continue from there...
  */
  {
    int bundle = 1;
    int fiber = 20;
    
    specex::Trace& trace = psf->FiberTraces[fiber];
    double wave = 4000;
    double x_ccd = trace.X_vs_W.Value(wave);
    double y_ccd = trace.Y_vs_W.Value(wave);
    
    specex::Legendre2DPol_p pol = psf->ParamsOfBundles[bundle].AllParPolXW[0];
    
    double rx = 2*(x_ccd-pol->xmin)/(pol->xmax-pol->xmin)-1;
    double ry = 2*(y_ccd-pol->ymin)/(pol->ymax-pol->ymin)-1;
    
    SPECEX_INFO("wave x_ccd y_ccd " << wave << " " << x_ccd << " " << y_ccd);
    SPECEX_INFO("reduced x and y " << rx << " " << ry);
    
  }

  
  SPECEX_INFO("create list of spots ...");
  
  vector<specex::Spot_p> spots;
  
  double wave_step = 50;

  for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf->ParamsOfBundles.begin();
      bundle_it != psf->ParamsOfBundles.end(); ++bundle_it) {
    for(std::map<int,specex::Trace>::const_iterator fiber_it = psf->FiberTraces.begin(); 
	fiber_it != psf->FiberTraces.end(); ++fiber_it) {
      
      const specex::Trace& trace = fiber_it->second;
      double trace_wave_min = trace.X_vs_W.xmin;
      double trace_wave_max = trace.X_vs_W.xmax;
      double min_wave = wave_step*int(max(3500.,trace_wave_min)/wave_step+1);
      
      // loop on wavelength
      for(double wave = min_wave ; wave<= trace_wave_max; wave += wave_step) {
	
	specex::Spot_p spot(new specex::Spot());
	
	spot->xc = trace.X_vs_W.Value(wave);
	spot->yc = trace.Y_vs_W.Value(wave);
	spot->flux = 1;
	spot->wavelength = wave;
	spot->fiber = fiber_it->first;
	spot->fiber_bundle = bundle_it->first;
	spots.push_back(spot);
	
	//break; // DEBUG, only one spot per fiber
      }
    }
  }
  
  // remove continuum
  for(std::map<int,specex::PSF_Params>::iterator bundle_it = psf->ParamsOfBundles.begin();
      bundle_it != psf->ParamsOfBundles.end(); ++bundle_it) {
    bundle_it->second.ContinuumPol.coeff.clear();
  }


   // create image
  specex::image_data img(psf->ccd_image_n_cols,psf->ccd_image_n_rows);
  
  specex::image_data weight(psf->ccd_image_n_cols,psf->ccd_image_n_rows);
  for(size_t i=0;i<weight.data.size();i++) weight.data(i)=1;
  
  bool only_on_spots = false;
  bool only_psf_core = false;
  bool only_positive = false;

  SPECEX_INFO("computing  image ...");  
  specex::parallelized_compute_model_image(img,weight,psf,spots,only_on_spots,only_psf_core,only_positive);
  
  SPECEX_INFO("writing image ...");
  specex::write_new_fits_image(output_image_fits_filename,img);
  
  return EXIT_SUCCESS;
}
