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

  // create image
  specex::image_data img(psf->ccd_image_n_cols,psf->ccd_image_n_rows);
  

  // define frame
  int begin_i=psf->ccd_image_n_cols;
  int end_i=0;
  int begin_j=psf->ccd_image_n_rows;
  int end_j=0;
  
  // loop on fiber bundles
  for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf->ParamsOfBundles.begin();
      bundle_it != psf->ParamsOfBundles.end(); ++bundle_it) {
    
    
    const specex::PSF_Params& params_of_bundle = bundle_it->second;
    
    // loop on fibers in psf
    for(std::map<int,specex::Trace>::const_iterator fiber_it = psf->FiberTraces.begin(); 
	fiber_it != psf->FiberTraces.end(); ++fiber_it) {
      
      int fiber_id = fiber_it->first;
      if(fiber_id<params_of_bundle.fiber_min || fiber_id>params_of_bundle.fiber_max) continue;
      
      const specex::Trace& trace = fiber_it->second;
      
      int x0=trace.X_vs_W.Value(trace.X_vs_W.xmin);
      int x1=trace.X_vs_W.Value(trace.X_vs_W.xmax);
      int y0=trace.Y_vs_W.Value(trace.Y_vs_W.xmin);
      int y1=trace.Y_vs_W.Value(trace.Y_vs_W.xmax);
      
      

      begin_i=min(begin_i,min(x0,x1));
      end_i=max(end_i,max(x0,x1));
      begin_j=min(begin_j,min(y0,y1));
      end_j=max(end_j,max(y0,y1));
    }
  }
  begin_i=max(0,begin_i-psf->hSizeX);
  end_i=min(int(psf->ccd_image_n_cols),end_i+psf->hSizeX+1);
  begin_j=max(0,begin_j-psf->hSizeY);
  end_j=min(int(psf->ccd_image_n_rows),end_j+psf->hSizeY+1);
  
  SPECEX_INFO("frame x=" << begin_i << "," << end_i << " y=" << begin_j << "," << end_j);
  
  
  // loop on fiber bundles
  for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf->ParamsOfBundles.begin();
      bundle_it != psf->ParamsOfBundles.end(); ++bundle_it) {
    
    int bundle_id = bundle_it->first;
    const specex::PSF_Params& params_of_bundle = bundle_it->second;
    
    // loop on fibers in psf
    for(std::map<int,specex::Trace>::const_iterator fiber_it = psf->FiberTraces.begin(); 
	fiber_it != psf->FiberTraces.end(); ++fiber_it) {
      
      int fiber_id = fiber_it->first;
      if(fiber_id<params_of_bundle.fiber_min || fiber_id>params_of_bundle.fiber_max) continue;
      
      SPECEX_INFO("filling image for fiber "<< fiber_id);

      const specex::Trace& trace = fiber_it->second;
      
      double trace_wave_min = trace.X_vs_W.xmin;
      double trace_wave_max = trace.X_vs_W.xmax;
      double min_wave = 50*int(max(3500.,trace_wave_min)/50+1);
      
      // loop on wavelength
      for(double wave = min_wave ; wave<= trace_wave_max; wave += 50) {
	
	// find center of spot for this wavelength and fiber
	double x_ccd = trace.X_vs_W.Value(wave);
	double y_ccd = trace.Y_vs_W.Value(wave);
	
	// if(fiber_id==20) cout << wave << " " << x_ccd << " " << y_ccd << " " << fiber_id << endl;
	
	// get psf parameters at this location in ccd
	
	
	int x_begin = int(x_ccd)-psf->hSizeX; // included 
	int x_end   = x_begin+2*psf->hSizeX+1; // not included
	int y_begin = int(y_ccd)-psf->hSizeY; // included 
	int y_end   = y_begin+2*psf->hSizeY+1; // not included
	
	// to go faster, pre-compute PSF parameters
	
	harp::vector_double psf_params = psf->AllLocalParamsXW(x_ccd,wave,bundle_id);
	
#ifdef EXTERNAL_TAIL
	double tail_amp = psf->r_tail_amplitude(wave);
#endif
	// loop on pixels of PSF footprint
	for(int j=y_begin;j<y_end;++j)
	  for(int i=x_begin;i<x_end;++i)
	    img(i,j) += psf->PSFValueWithParamsXY(x_ccd,y_ccd,i,j,psf_params,0,0);
	
#ifdef EXTERNAL_TAIL
	for(int j=begin_j;j<end_j;j++)
	  for(int i=begin_i;i<end_i;i++)
	    img(i,j) += tail_amp*psf->TailProfile(i-x_ccd,j-y_ccd);
#endif	
	

      } // end of loop on wavelength
    } // end of loop on fibers
  } // end of loop on fiber bundles
  
  // write image
  SPECEX_INFO("writing image ...");
  specex::write_new_fits_image(output_image_fits_filename,img);
  
  

  return EXIT_SUCCESS;
}
