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
#include <specex_model_image.h>



using namespace std;
using namespace specex;

namespace popts = boost::program_options;




int main ( int argc, char *argv[] ) {
  

  

  // default arguments
  // --------------------------------------------
  string spectrograph_name = "BOSS";
  
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
    Stamp global_stamp = compute_stamp(image,psf,spots,0,0);
    
    // fill model image
    image_data model(image.n_cols(),image.n_rows());
    
    image_data core_footprint(image.n_cols(),image.n_rows());
    int hsize=2;   
    for(size_t s=0;s<spots.size();s++) {     
      specex::Spot_p spot = spots[s];
      
      int i = int(floor(spot->xc+0.5));
      int j = int(floor(spot->yc+0.5));
      int begin_i = max(0,i-hsize);
      int end_i = min(int(image.n_cols()),i+hsize+1);
      int begin_j = max(0,j-hsize);
      int end_j = min(int(image.n_rows()),j+hsize+1);
      for(j=begin_j;j<end_j;j++)
	for(i=begin_i;i<end_i;i++)
	  core_footprint(i,j)=1;
    }
    
    bool only_on_spots = false;
    bool only_psf_core = false;
    bool only_positive = false;
    
    parallelized_compute_model_image(model,weight,psf,spots,only_on_spots,only_psf_core,only_positive,0,0);
    
    const image_data& model_for_var = model;


    // now compute variance 
    image_data variance(image.n_cols(),image.n_rows());
    image_data variance2(image.n_cols(),image.n_rows());
    image_data data_in_stamp(image.n_cols(),image.n_rows());
    
    for(size_t i=0; i<model.data.size() ;i++) {
      double flux = model_for_var.data(i);
      if(flux==0) continue; // never been on this pixel
      
      data_in_stamp.data(i) = image.data(i);
      
      variance.data(i) = square(readout_noise);
      if(flux>0)
	variance.data(i) += (flux + square(psf_error*flux));
	  
      // use data here for the variance : 
      flux = max(0.,image.data(i));
      
      variance2.data(i) = square(readout_noise);
      if(flux>0) 
	variance2.data(i) += (flux + square(psf_error*flux));
      
    } // end of loop on all pixels
    

    // compute chi2 (image and core)
    // ---------------------------
    image_data residual = image; 

    // zero if no model
    for(size_t i=0;i<residual.data.size();i++)
      if(model.data(i)==0) 
	residual.data(i)=0;
      else
	residual.data(i) -= model.data(i);
    
    
    double chi2_image = 0;
    int ndata_image = 0;
    double chi2_core = 0;
    double chi2_core2 = 0;
    int ndata_core = 0;
    for (int j=global_stamp.begin_j; j <global_stamp.end_j; ++j) {  

      
      
      int begin_i = global_stamp.begin_i;
      int end_i   = global_stamp.end_i;

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
	  }
	}
      }
    }
    
    cout << "psf" << " hx=" << psf->hSizeX << " hy=" <<  psf->hSizeY << " lpar=" << psf->LocalNAllPar() << " gnpar=" << psf->BundleNAllPar(spots[0]->fiber_bundle) << endl;
    
    int npar = 0;
    for(std::map<int,PSF_Params>::iterator bundle_it = psf->ParamsOfBundles.begin(); bundle_it != psf->ParamsOfBundles.end(); bundle_it++) {

      npar += psf->BundleNAllPar(bundle_it->first)+spots.size(); 
      // psf->TracesNPar();
      npar -= bundle_it->second.AllParPolXW[0]->coeff.size(); // not in last fit
      npar -= bundle_it->second.AllParPolXW[1]->coeff.size(); // not in last fit
    }
    cout << "npar in last fit = " << npar << endl;
    
    int ndf_image  = ndata_image-npar; 
    int ndf_core   = ndata_core-npar;
    cout << "readout noise      = " << readout_noise << endl;
    cout << "assumed psf error  = " << psf_error << endl;
    cout << "number of spots    = " << spots.size() << endl;
    cout << "image chi2/ndf     = " << chi2_image << "/" << ndf_image << " = " << chi2_image/ndf_image << endl;
    cout << "spot core 5x5 chi2/ndf (using model in var) = " << chi2_core << "/" << ndf_core << " =" << chi2_core/ndf_core << endl;
    cout << "spot core 5x5 chi2/ndf (using data in var)  = " << chi2_core2 << "/" << ndf_core << " =" << chi2_core2/ndf_core << endl;
    cout << "ndata in image = " << ndata_image << endl;
    cout << "ndata in core  = " << ndata_core << endl;
    
    // write images of model and residuals
    // ---------------------------
    if(output_fits_image_filename!="") {
      fitsfile * fp;  
      harp::fits::create ( fp, output_fits_image_filename );
      
      harp::fits::img_append < double > ( fp, image.n_rows(), image.n_cols() );
      harp::fits::img_write ( fp, image.data );
      harp::fits::key_write(fp,"WHAT","DATA","");
      
      harp::fits::img_append < double > ( fp, image.n_rows(), image.n_cols() );
      harp::fits::img_write ( fp, model.data );
      harp::fits::key_write(fp,"EXTNAME","MODEL","");
  
      harp::fits::img_append < double > ( fp, image.n_rows(), image.n_cols() );
      harp::fits::img_write ( fp, residual.data );
      harp::fits::key_write(fp,"EXTNAME","RESIDUAL","");
        
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
  }
  /*
  catch(std::exception e) {
    cerr << "FATAL ERROR " << e.what() << endl;
    return EXIT_FAILURE;
  }
  */
  return EXIT_SUCCESS;
}



