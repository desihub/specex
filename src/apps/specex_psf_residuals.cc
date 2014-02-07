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
#include <specex_hat_hermite_psf.h>
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
    Stamp global_stamp(image);
    global_stamp.begin_i = 10000;
    global_stamp.end_i   = 0;
    global_stamp.begin_j = 10000;
    global_stamp.end_j   = 0;
    
    vector<specex::Stamp> spot_stamps;
    for(size_t s=0;s<spots.size();s++) {
      specex::Spot_p spot = spots[s];
      
      Stamp stamp(image);
      psf->StampLimits(spot->xc,spot->yc,stamp.begin_i,stamp.end_i,stamp.begin_j,stamp.end_j);
      stamp.begin_i = max(0,stamp.begin_i);
      stamp.end_i   = min(stamp.Parent_n_cols(),stamp.end_i);
      stamp.begin_j = max(0,stamp.begin_j);
      stamp.end_j   = min(stamp.Parent_n_rows(),stamp.end_j);
      
      global_stamp.begin_i = min(global_stamp.begin_i,stamp.begin_i);
      global_stamp.end_i   = max(global_stamp.end_i,stamp.end_i);
      global_stamp.begin_j = min(global_stamp.begin_j,stamp.begin_j);
      global_stamp.end_j   = max(global_stamp.end_j,stamp.end_j);
      
      spot_stamps.push_back(stamp);
    }
    
    
    // fill model image
    image_data model(image.n_cols(),image.n_rows());
    image_data model_for_var(image.n_cols(),image.n_rows());
    image_data core_footprint(image.n_cols(),image.n_rows());
    
    const PSF_Params * psf_params = &(psf->ParamsOfBundles.find(spots[0]->fiber_bundle)->second);

#ifdef CONTINUUM
    
    for (int j=global_stamp.begin_j; j <global_stamp.end_j; ++j) {  

      int begin_i = max(global_stamp.begin_i,int(floor(psf->GetTrace(psf_params->fiber_min).X_vs_Y.Value(double(j))+0.5))-psf->hSizeX-1);
      int end_i   = min(global_stamp.end_i,int(floor(psf->GetTrace(psf_params->fiber_max).X_vs_Y.Value(double(j))+0.5))+psf->hSizeX+2);
      


      for(int fiber=psf_params->fiber_min;fiber<=psf_params->fiber_max;fiber++) {
	double x = psf->GetTrace(fiber).X_vs_Y.Value(double(j));
	double w = psf->GetTrace(fiber).W_vs_Y.Value(double(j));
	double continuum_flux = psf->ContinuumPol.Value(w);
	double expfact_for_continuum=continuum_flux/(2*M_PI*square(psf->continuum_sigma_x));
	if(expfact_for_continuum!=0) {
	  for (int i=begin_i ; i <end_i; ++i) {    
	    if(weight(i,j)<=0) continue;
	    double val = expfact_for_continuum*exp(-0.5*square((i-x)/psf->continuum_sigma_x));
	    model(i,j) += val;
	    if(val>0) model_for_var(i,j) += val;
	  } // end of loop on i
	}
      } // end of loop on fiber
    } // end of loop on j

#endif  
    
    
    
#ifdef EXTERNAL_TAIL    
    for(size_t s=0;s<spots.size();s++) {
      specex::Spot_p spot = spots[s];
      
      double r_tail_amplitude = spot->flux*psf->RTailAmplitudePol.Value(spot->wavelength);
      double r_tail_amplitude_for_var = max(r_tail_amplitude,1e-20);
      
      if(spot->flux<0) r_tail_amplitude = 0; // this is frozen_flux=0 in fitter

      // first fill tails on stamp
      for (int j=global_stamp.begin_j; j <global_stamp.end_j; ++j) {  
	
	int begin_i = max(global_stamp.begin_i,int(floor(psf->GetTrace(psf_params->fiber_min).X_vs_Y.Value(double(j))+0.5))-psf->hSizeX-1);
	int end_i   = min(global_stamp.end_i,int(floor(psf->GetTrace(psf_params->fiber_max).X_vs_Y.Value(double(j))+0.5))+psf->hSizeX+2);
	
	
	for (int i=begin_i ; i < end_i; ++i) {
	  if(weight(i,j)<=0) continue;
	  double prof = psf->TailProfile(i-spot->xc,j-spot->yc);
	  model(i,j) += r_tail_amplitude*prof;
	  model_for_var(i,j) += r_tail_amplitude_for_var*prof;
	  
	}
      } // end of loop on stamp pixels
    }
#endif
    

    
    for(size_t s=0;s<spots.size();s++) {
      specex::Spot_p spot = spots[s];
       
      harp::vector_double params = psf->AllLocalParamsXW(spot->xc,spot->wavelength,spot->fiber_bundle);
      
      
      // then the core of the psf for this spot's stamp only
      const Stamp& spot_stamp = spot_stamps[s];
      for (int j=spot_stamp.begin_j; j <spot_stamp.end_j; ++j) { 
	
	
	int begin_i = max(spot_stamp.begin_i,int(floor(psf->GetTrace(psf_params->fiber_min).X_vs_Y.Value(double(j))+0.5))-psf->hSizeX-1);
	int end_i   = min(spot_stamp.end_i,int(floor(psf->GetTrace(psf_params->fiber_max).X_vs_Y.Value(double(j))+0.5))+psf->hSizeX+2);
	
	for (int i=begin_i ; i < end_i; ++i) {
	  
	  if(weight(i,j)<=0) continue;
	  double val = spot->flux*psf->PSFValueWithParamsXY(spot->xc,spot->yc, i, j, params, 0, 0);
	  model(i,j) += val;
	  model_for_var(i,j) += max(1e-20,val);
	  core_footprint(i,j) = 1;
	}
      } // end of loop on stamp pixels
    } // end of loop on spots
    
    // now compute variance 
    image_data variance(image.n_cols(),image.n_rows());
    image_data variance2(image.n_cols(),image.n_rows());
    image_data data_in_stamp(image.n_cols(),image.n_rows());
    
    for(size_t i=0; i<model.data.size() ;i++) {
      double flux = model_for_var.data(i);
      if(flux==0) continue; // never been on this pixel
      
      data_in_stamp.data(i) = image.data(i);
      
      variance.data(i) = square(readout_noise) + square(psf_error*flux); // readout and psf error
      if(flux>0) variance.data(i) += flux; // Poisson noise
      
      // use data here for the variance : 
      flux = max(0.,image.data(i));
      
      variance2.data(i) = square(readout_noise) + square(psf_error*flux); // readout and psf error
      if(flux>0) variance2.data(i) += flux; // Poisson noise
      
    } // end of loop on all pixels
    

    // compute chi2 (image and core)
    // ---------------------------
    image_data residual = image; 
    residual.data -= model.data;
    
    double chi2_image = 0;
    int ndata_image = 0;
    double chi2_core = 0;
    double chi2_core2 = 0;
    double chi2_core_trimed = 0;
    int ndata_core = 0;
    int ndata_core_trimed = 0;
    int npix_trim = 20;
    for (int j=global_stamp.begin_j; j <global_stamp.end_j; ++j) {  

      int begin_i = max(global_stamp.begin_i,int(floor(psf->GetTrace(psf_params->fiber_min).X_vs_Y.Value(double(j))+0.5))-psf->hSizeX-1);
      int end_i   = min(global_stamp.end_i,int(floor(psf->GetTrace(psf_params->fiber_max).X_vs_Y.Value(double(j))+0.5))+psf->hSizeX+2);

      for (int i=begin_i ; i <end_i; ++i) {
	if(variance(i,j)>0) {

	  // variance(i,j) = 1; // debug

	  double dchi2 = square(residual(i,j))/(variance(i,j));
	  
	  chi2_image += dchi2;
	  ndata_image ++;
	
	  if(core_footprint(i,j)>0) {
	    chi2_core += dchi2;
	    chi2_core2 += square(residual(i,j))/(variance2(i,j));
	    ndata_core ++;
	    if(i>global_stamp.begin_i+npix_trim && i<global_stamp.end_i-npix_trim) {
	      chi2_core_trimed += dchi2;
	      ndata_core_trimed ++;
	    }
	  }
	}
      }
    }
    
    cout << "psf" << " hx=" << psf->hSizeX << " hy=" <<  psf->hSizeY << " lpar=" << psf->LocalNAllPar() << " gnpar=" << psf->BundleNAllPar(spots[0]->fiber_bundle) << endl;
    int npar = psf->BundleNAllPar(spots[0]->fiber_bundle)+spots.size(); 
    // psf->TracesNPar();
    npar -= psf_params->AllParPolXW[0]->coeff.size(); // not in last fit
    npar -= psf_params->AllParPolXW[1]->coeff.size(); // not in last fit
    cout << "npar in last fit = " << npar << endl;

    if(0) { // not included in final fit 
#ifdef EXTERNAL_TAIL
      npar += psf->RTailAmplitudePol.coeff.size();
#endif
#ifdef CONTINUUM
      npar += psf->ContinuumPol.coeff.size();
#endif
    }


    int ndf_image  = ndata_image-npar; 
    int ndf_core   = ndata_core-npar;
    int ndf_core_trimed   = ndata_core_trimed-npar;
    cout << "readout noise      = " << readout_noise << endl;
    cout << "assumed psf error  = " << psf_error << endl;
    cout << "number of spots    = " << spots.size() << endl;
    cout << "image chi2/ndf     = " << chi2_image << "/" << ndf_image << " = " << chi2_image/ndf_image << endl;
    cout << "spot core chi2/ndf (using model in var) = " << chi2_core << "/" << ndf_core << " =" << chi2_core/ndf_core << endl;
    cout << "spot core chi2/ndf (using data in var)  = " << chi2_core2 << "/" << ndf_core << " =" << chi2_core2/ndf_core << endl;
    if(ndf_core_trimed>0)
      cout << "spot core chi2/ndf (using model in var, trimed,  " << npix_trim << " pixs)  = " << chi2_core_trimed << "/" << ndf_core_trimed << " =" << chi2_core2/ndf_core_trimed << endl;
    cout << "ndata in image = " << ndata_image << endl;
    cout << "ndata in core  = " << ndata_core << endl;
    
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



