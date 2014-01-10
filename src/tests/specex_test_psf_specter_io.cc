#include <iostream>
#include <fstream>
#include <ctime>

#include "harp.hpp"

#include "specex_message.h"
#include "specex_fits.h"
#include "specex_psf.h"
#include "specex_psf_io.h"
#include "specex_gauss_hermite_psf.h"
#include "specex_image_data.h"
#include "specex_serialization.h"

using namespace std;

int main() {
  

  specex_set_verbose(true);
  
  specex::PSF_p psf;
  specex::read_psf_xml(psf,"psf.xml");
  
  if(0){
    // testing coordinate system  
    double xc = 400;
    double yc = 2000;
    harp::vector_double P(psf->LocalNPar());
    P *= 0;
    int i=floor(xc);
    int j = floor(yc);
    //for(i=floor(xc)-3;i<=floor(xc)+3;i++)
    for(int j=floor(yc)-3;j<=floor(yc)+3;j++) 
      cout << i << " " << j << " " << psf->PSFValueWithParamsXY(400,2000,i,j,P,NULL,NULL) << endl;
    
    return 0;
  }
  
  // write psf in specter format
  psf->WriteFits("test-psf-for-specter2.fits");

  // write an image of the psf at a precise location
  write_psf_fits_image(psf,"psf-image.fits",400,2000,1,4);
  
  {
    specex::PSF_p ipsf(new specex::GaussHermitePSF());
    // read psf in specter format
    
    ipsf->ReadFits("test-psf-for-specter2.fits");

    write_psf_fits_image(ipsf,"psf-image-bis.fits",400,2000,1,4);
  }

  
  

  // now create a test image
  cout << "Generate an image with spots of flux=1 at each 50A starting at 3500A for each fiber addressed by PSF" << endl;
  
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
    
    specex::Legendre2DPol& pol = psf->ParamsOfBundles[bundle].Polynomials[0];
    
    double rx = 2*(x_ccd-pol.xmin)/(pol.xmax-pol.xmin)-1;
    double ry = 2*(y_ccd-pol.ymin)/(pol.ymax-pol.ymin)-1;
    
    cout << "wave x_ccd y_ccd " << wave << " " << x_ccd << " " << y_ccd << endl;
    cout << "reduced x and y " << rx << " " << ry << endl;
    
  }

  // create image
  specex::image_data img(psf->ccd_image_n_cols,psf->ccd_image_n_rows);
  

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
      
      const specex::Trace& trace = fiber_it->second;
      
      double trace_wave_min = trace.X_vs_W.xmin;
      double trace_wave_max = trace.X_vs_W.xmax;
      double min_wave = 50*int(max(3500.,trace_wave_min)/50+1);
      
      // loop on wavelength
      for(double wave = min_wave ; wave<= trace_wave_max; wave += 50) {
	
	// find center of spot for this wavelength and fiber
	double x_ccd = trace.X_vs_W.Value(wave);
	double y_ccd = trace.Y_vs_W.Value(wave);
	
	if(fiber_id==20) cout << wave << " " << x_ccd << " " << y_ccd << " " << fiber_id << endl;
	
	// get psf parameters at this location in ccd
	
	
	int x_begin = int(x_ccd)-psf->hSizeX; // included 
	int x_end   = x_begin+2*psf->hSizeX+1; // not included
	int y_begin = int(y_ccd)-psf->hSizeY; // included 
	int y_end   = y_begin+2*psf->hSizeY+1; // not included
	
	// to go faster, pre-compute PSF parameters
	
	harp::vector_double psf_params = psf->LocalParamsXW(x_ccd,wave,bundle_id);
	
	// loop on pixels of PSF footprint
	for(int j=y_begin;j<y_end;++j)
	  for(int i=x_begin;i<x_end;++i)
	    img(i,j) += psf->PSFValueWithParamsXY(x_ccd,y_ccd,i,j,psf_params,0,0);
	
      } // end of loop on wavelength
    } // end of loop on fibers
  } // end of loop on fiber bundles
  
  // write image
  specex::write_new_fits_image("test-ccd-image-for-specter2.fits",img);
  
  

  return EXIT_SUCCESS;
}
