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
    vector<specex::Spot_p> input_spots;
    {
      std::ifstream is(spots_xml_filename.c_str());
      boost::archive::xml_iarchive xml_ia ( is );
      xml_ia >> BOOST_SERIALIZATION_NVP(input_spots);
      is.close();
    } 
    
    SPECEX_INFO("number of input spots = " << input_spots.size());

    vector<specex::Spot_p> spots;
    for(size_t s=0;s<input_spots.size();s++) {
      specex::Spot_p& spot = input_spots[s];
      if(spot->status == 0) continue;
      spots.push_back(spot);
    }
    
    SPECEX_INFO("number of selected spots = " << spots.size());
    
    // create output images
    // --------------------------------------------
    // definition of fitted region of image
    // ----------------------------------------------------
    double xc_min=1e20;
    double xc_max=-1e20;
    double yc_min=1e20;
    double yc_max=-1e20;
    for(size_t s=0;s<spots.size();++s) {
      const specex::Spot_p& spot = spots[s];
      if(spot->xc<xc_min) xc_min=spot->xc;
      if(spot->xc>xc_max) xc_max=spot->xc;
      if(spot->yc<yc_min) yc_min=spot->yc;
      if(spot->yc>yc_max) yc_max=spot->yc;
    }
    Stamp global_stamp(image);
    global_stamp.begin_i = max(global_stamp.begin_i,int(xc_min+0.5)-psf->hSizeX);
    global_stamp.end_i   = min(global_stamp.end_i,int(xc_max+0.5)+psf->hSizeX+1); 
    global_stamp.begin_j = max(global_stamp.begin_j,int(yc_min+0.5)-psf->hSizeY);
    global_stamp.end_j   = min(global_stamp.end_j,int(yc_max+0.5)+psf->hSizeY+1); 
    // ----------------------------------------------------
    
    
    // create a list of stamps
    vector<specex::Stamp> spot_stamps;
    for(size_t s=0;s<spots.size();s++) {
      specex::Spot_p spot = spots[s];
      
      Stamp stamp(image);
      psf->StampLimits(spot->xc,spot->yc,stamp.begin_i,stamp.end_i,stamp.begin_j,stamp.end_j);
      stamp.begin_i = max(0,stamp.begin_i);
      stamp.end_i   = min(stamp.Parent_n_cols(),stamp.end_i);
      stamp.begin_j = max(0,stamp.begin_j);
      stamp.end_j   = min(stamp.Parent_n_rows(),stamp.end_j);
      spot_stamps.push_back(stamp);
    }

    
    // fill model image
    image_data model(image.n_cols(),image.n_rows());
    specex::zero(model.data);

#ifdef CONTINUUM

    const PSF_Params * psf_params = &(psf->ParamsOfBundles.find(spots[0]->fiber_bundle)->second);
    
    for (int j=global_stamp.begin_j; j <global_stamp.end_j; ++j) {  
      
      for(int fiber=psf_params->fiber_min;fiber<=psf_params->fiber_max;fiber++) {
	double x = psf->GetTrace(fiber).X_vs_Y.Value(double(j));
	double w = psf->GetTrace(fiber).W_vs_Y.Value(double(j));
	double continuum_flux = psf->ContinuumPol.Value(w);
	double expfact_for_continuum=continuum_flux/(2*M_PI*square(psf->continuum_sigma_x));
	for (int i=global_stamp.begin_i ; i <global_stamp.end_i; ++i) {
	  
	  if(weight(i,j)<=0) continue;
	  model(i,j) += expfact_for_continuum*exp(-0.5*square((i-x)/psf->continuum_sigma_x));
	} // end of loop on i
      } // end of loop on fiber
    } // end of loop on j

#endif  

    
    for(size_t s=0;s<spots.size();s++) {
      specex::Spot_p spot = spots[s];
      
      harp::vector_double params = psf->AllLocalParamsXW(spot->xc,spot->wavelength,spot->fiber_bundle);
      
#ifdef EXTERNAL_TAIL
      
      double r_tail_amplitude = psf->RTailAmplitudePol.Value(spot->wavelength);

      // first fill tails on stamp
      for (int j=global_stamp.begin_j; j <global_stamp.end_j; ++j) {  
	for (int i=global_stamp.begin_i ; i < global_stamp.end_i; ++i) {
	  if(weight(i,j)<=0) continue;
	  model(i,j) += spot->flux*psf->TailValueA(r_tail_amplitude,i-spot->xc,j-spot->yc);
	}
      } // end of loop on stamp pixels
      
#endif



      // then the core of the psf for this spot's stamp only
      const Stamp& spot_stamp = spot_stamps[s];
      for (int j=spot_stamp.begin_j; j <spot_stamp.end_j; ++j) {  
	for (int i=spot_stamp.begin_i ; i < spot_stamp.end_i; ++i) {
	  
	  if(weight(i,j)<=0) continue;
	  model(i,j) += spot->flux*psf->PSFValueWithParamsXY(spot->xc,spot->yc, i, j, params, 0, 0);
	}
      } // end of loop on stamp pixels
    } // end of loop on spots

    // now compute variance 
    image_data variance(image.n_cols(),image.n_rows());
    specex::zero(variance.data);
    image_data data_in_stamp(image.n_cols(),image.n_rows());
    specex::zero(data_in_stamp.data);
    
    

    for(size_t i=0; i<model.data.size() ;i++) {
      double flux = model.data(i);
      if(flux==0) continue; // never been on this pixel

      // use data here for the variance : 
      flux = max(0.,image.data(i));
      
      variance.data(i) = square(readout_noise) + square(psf_error*flux); // readout and psf error
      if(flux>0) variance.data(i) += flux; // Poisson noise
    } // end of loop on all pixels
    

    // compute chi2 (image and core)
    // ---------------------------
    image_data residual = image; 
    residual.data -= model.data;

    double chi2_image = 0;
    int ndata_image = 0;
    for(size_t i=0; i<residual.data.size() ;i++) {
      if(variance.data(i)>0) {
	chi2_image += square(residual.data(i))/(variance.data(i));
	ndata_image ++;
      }
    }
    double chi2_core = 0;
    int ndata_core = 0;
    for(size_t s=0;s<spots.size();s++) {
      const Stamp& spot_stamp = spot_stamps[s];
      for (int j=spot_stamp.begin_j; j <spot_stamp.end_j; ++j) {  
	for (int i=spot_stamp.begin_i ; i < spot_stamp.end_i; ++i) {
	  if(variance(i,j)>0) {
	    chi2_core += square(residual(i,j))/(variance(i,j));
	    ndata_core ++;
	  }
	  
	}
      } // end of loop on stamp pixels
    }
    cout << "psf" << " hx=" << psf->hSizeX << " hy=" <<  psf->hSizeY << " lpar=" << psf->LocalNAllPar() << " gnpar=" << psf->BundleNAllPar(spots[0]->fiber_bundle) << endl;
    int npar = psf->BundleNAllPar(spots[0]->fiber_bundle)+spots.size()+psf->TracesNPar();
    int ndf_image  = ndata_image-npar; 
    int ndf_core   = ndata_core-npar; 
    cout << "image chi2/ndf     = " << chi2_image/ndf_image << " ndf=" << ndf_image << " ndata=" << ndata_image << endl;
    cout << "spot core chi2/ndf = " << chi2_core/ndf_core << " ndf=" << ndf_core << " ndata=" << ndata_core << endl;
    
    
    // write images of model and residuals
    // ---------------------------
    if(output_fits_image_filename!="") {
      fitsfile * fp;  
      harp::fits::create ( fp, output_fits_image_filename );
      harp::fits::img_append < double > ( fp, image.n_rows(), image.n_cols() );
      harp::fits::img_write ( fp, model.data );
      harp::fits::key_write(fp,"WHAT","MODEL","");
  
      harp::fits::img_append < double > ( fp, image.n_rows(), image.n_cols() );
      harp::fits::img_write ( fp, residual.data );
      harp::fits::key_write(fp,"EXTNAME","RESALL","");
  
      residual = data_in_stamp; 
      residual.data -= model.data;

      harp::fits::img_append < double > ( fp, image.n_rows(), image.n_cols() );
      harp::fits::img_write ( fp, residual.data );
      harp::fits::key_write(fp,"EXTNAME","RES","");
  
      for(size_t i=0; i<residual.data.size() ;i++) {
	if(variance.data(i)>0) residual.data(i) /= sqrt(variance.data(i));
      }
    
      harp::fits::img_append < double > ( fp, image.n_rows(), image.n_cols() );
      harp::fits::img_write ( fp, residual.data );
      harp::fits::key_write(fp,"EXTNAME","PULL","");
  
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



