#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <map>

#include <boost/program_options.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>

#include <harp.hpp>

#include <specex_message.h>
#include <specex_psf.h>
#include <specex_trace.h>
#include <specex_spot.h>
#include <specex_spot_array.h>
#include <specex_fits.h>
#include <specex_gauss_hermite_psf.h>
#include <specex_serialization.h>
#include <specex_image_data.h>
#include <specex_stamp.h>



using namespace std;
using namespace specex;

namespace popts = boost::program_options;




int main ( int argc, char *argv[] ) {
  

  

  // default arguments
  // --------------------------------------------
  string spectrograph_name = "BOSS";
  int    first_fiber_bundle=1;
  int    last_fiber_bundle=1;
  
  string psf_xml_filename="";
  string spots_xml_filename="";
  string input_fits_image_filename="";
  string output_fits_image_filename="";  
  double psf_error=0;
  double readout_noise=2;
  
  // reading arguments
  // --------------------------------------------
  popts::options_description desc ( "Allowed Options" );
  desc.add_options()
    ( "help,h", "display usage information" )
    ( "psf", popts::value<string>( &psf_xml_filename ), "psf xml filename" )
    ( "spots", popts::value<string>( &spots_xml_filename ), "spots xml filename" )
    ( "in", popts::value<string>( &input_fits_image_filename ), " input fits image file name" )
    ( "out", popts::value<string>( &output_fits_image_filename ), " output fits image file name")
    ( "readout-noise", popts::value<double>( &readout_noise ), " readout noise for pull")
    ( "psf-error", popts::value<double>( &psf_error ), " psf relative error for pull")
    ( "verbose,v", "turn on verbose mode" )
    ( "core", "dump core files when harp exception is thrown" )
    ;
  
  popts::variables_map vm;
  
  try {
  
    
    popts::store(popts::command_line_parser( argc, argv ).options(desc).run(), vm);
    popts::notify(vm);
    
    if ( ( argc < 2 ) || vm.count( "help" ) || ( ! vm.count( "psf" ) )  || ( ! vm.count( "spots" ) ) 
	 || ( ! vm.count( "in" ) )  ) {
      cerr << endl;
      cerr << desc << endl;
      cerr << "example:" << endl;
      cerr << argv[0] << " --psf psf.xml  --spots spots.xml --in sdProc-b1-00108382.fits --out res.fits" << endl;
      return EXIT_FAILURE;
    }
  }catch(std::exception e) {
    cerr << "error in arguments" << endl;
    cerr << endl;
    cerr << desc << endl;
    return EXIT_FAILURE;
  }
  try {

    specex_set_verbose(vm.count("verbose")>0);
    specex_set_dump_core(vm.count("core")>0);
    
  

    
    
    // open image
    // --------------------------------------------
    image_data image,weight;
    read_fits_images(input_fits_image_filename,image,weight);
    // weight is 0 or 1
    for(int j=0;j<weight.Ny();j++) {
      for(int i=0;i<weight.Nx();i++) {
	if(weight(i,j)>0)
	  weight(i,j)=1;
	else
	  weight(i,j)=0;
      }
    }
    
    // read PSF
    // --------------------------------------------
    specex::PSF_p psf;
    {
      std::ifstream is(psf_xml_filename.c_str());
      boost::archive::xml_iarchive xml_ia ( is );
      xml_ia >> BOOST_SERIALIZATION_NVP(psf);
      is.close();
    }
    
    
    // read spots
    // --------------------------------------------
    vector<specex::Spot_p> spots;
    {
      std::ifstream is(spots_xml_filename.c_str());
      boost::archive::xml_iarchive xml_ia ( is );
      xml_ia >> BOOST_SERIALIZATION_NVP(spots);
      is.close();
    } 
    
    SPECEX_INFO("number of spots = " << spots.size());
    
    // create output images
    // --------------------------------------------
    image_data model(image.n_cols(),image.n_rows());
    image_data data_in_stamp(image.n_cols(),image.n_rows());
    image_data variance(image.n_cols(),image.n_rows());
    
    model.data *= 0;
    data_in_stamp.data *= 0;
    variance.data *= 0;

    for(size_t s=0;s<spots.size();s++) {
      specex::Spot_p spot = spots[s];
      
      Stamp stamp(image);
      psf->StampLimits(spot->xc,spot->yc,stamp.begin_i,stamp.end_i,stamp.begin_j,stamp.end_j);
      stamp.begin_i = max(0,stamp.begin_i);
      stamp.end_i   = min(stamp.Parent_n_cols(),stamp.end_i);
      stamp.begin_j = max(0,stamp.begin_j);
      stamp.end_j   = min(stamp.Parent_n_rows(),stamp.end_j);

      harp::vector_double params = psf->LocalParamsXW(spot->xc,spot->wavelength,spot->fiber_bundle);
            
      for (int j=stamp.begin_j; j <stamp.end_j; ++j) {  
	for (int i=stamp.begin_i ; i < stamp.end_i; ++i) {

	  if(weight(i,j)<=0) continue;

	  double val =  spot->flux*psf->PSFValueWithParamsXY(spot->xc,spot->yc, i, j, params, 0, 0);
	  model(i,j) += val;
	  data_in_stamp(i,j) = image(i,j);
	  variance(i,j) += square(readout_noise) + square(psf_error*val);
	  if(val>0) variance(i,j) += val;
	}
      } // end of loop on stamp pixels
    } // end of loop on spots

    // compute chi2
    // ---------------------------
    image_data residual = image; 
    residual.data -= model.data;

    double chi2 = 0;
    int ndata = 0;
    for(size_t i=0; i<residual.data.size() ;i++) {
      if(variance.data(i)>0) {
	chi2 += square(residual.data(i))/(variance.data(i));
	ndata ++;
      }
    }
    cout << "psf" << " hx=" << psf->hSizeX << " hy=" <<  psf->hSizeX << " lpar=" << psf->LocalNPar() << " gnpar=" << psf->BundleNPar(spots[0]->fiber_bundle) << endl;
    int npar = psf->BundleNPar(spots[0]->fiber_bundle)+spots.size()+psf->TracesNPar();
    int ndf  = ndata-npar; 
    cout << "chi2/ndf = " << chi2/ndf << " ndf=" << ndf << " ndata=" << ndata << endl;

    
    // write images of model and residuals
    // ---------------------------
    if(output_fits_image_filename!="") {
      fitsfile * fp;  
      harp::fits::create ( fp, output_fits_image_filename );
      harp::fits::img_append < double > ( fp, image.n_rows(), image.n_cols() );
      harp::fits::img_write ( fp, model.data );
      
      harp::fits::img_append < double > ( fp, image.n_rows(), image.n_cols() );
      harp::fits::img_write ( fp, residual.data );
      
      residual = data_in_stamp; 
      residual.data -= model.data;

      harp::fits::img_append < double > ( fp, image.n_rows(), image.n_cols() );
      harp::fits::img_write ( fp, residual.data );
      
      for(size_t i=0; i<residual.data.size() ;i++) {
	if(variance.data(i)>0) residual.data(i) /= sqrt(variance.data(i));
      }
    
      harp::fits::img_append < double > ( fp, image.n_rows(), image.n_cols() );
      harp::fits::img_write ( fp, residual.data );
      
      harp::fits::close ( fp );
    }
    
    

  // ending
  // --------------------------------------------
  }catch(harp::exception e) {
    cerr << "FATAL ERROR (harp) " << e.what() << endl;
    return EXIT_FAILURE;
  }catch(std::exception e) {
    cerr << "FATAL ERROR " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}



