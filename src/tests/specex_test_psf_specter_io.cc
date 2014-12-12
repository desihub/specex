#include <iostream>
#include <fstream>
#include <ctime>

#include <boost/program_options.hpp>

#include <harp.hpp>

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
  string output_psf_fits_filename="";
  
  specex_set_verbose(true);
  
  popts::options_description desc ( "Allowed Options" );
  desc.add_options()
    ( "help,h", "display usage information" )
    ( "in_xml", popts::value<string>( & input_psf_xml_filename ), "input PSF xml file name" )
    ( "in_fits", popts::value<string>( & input_psf_fits_filename ), "input PSF fits file name" )
    ( "out_fits_image", popts::value<string>( & output_image_fits_filename ), "output image fits file name" )
    ( "out_fits_psf", popts::value<string>( & output_psf_fits_filename ), "output psf fits file name" )
    ;
  
  popts::variables_map vm;
  try {
    popts::store(popts::command_line_parser( argc, argv ).options(desc).run(), vm);
    popts::notify(vm);
    
    if ( ( argc < 2 ) || vm.count( "help" ) || ( ! vm.count( "in_xml" ) && ( ! vm.count( "in_fits" ) ) )  ||  ( ! vm.count( "out_fits_image") && ( ! vm.count( "out_fits_psf" ) ) ) ) {
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
  
  /*
  if(psf->Name() !=  "GaussHermite2PSF") {
    cout << "code only works with GAUSSHERMITE2PSF" << endl;
    return 0;
  }
  */

  
  if(vm.count( "out_fits_psf" ))
  specex::write_psf_fits(psf,output_psf_fits_filename);
  

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
  
  
  // now create a test image
  SPECEX_INFO("Generate an image with spots of flux=1 at each 50A starting at 3500A for each fiber addressed by PSF");
  
  
  
  SPECEX_INFO("create list of spots ...");
  
  vector<specex::Spot_p> spots;
  
  double wave_step = 50;

  for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf->ParamsOfBundles.begin();
      bundle_it != psf->ParamsOfBundles.end(); ++bundle_it) {
    for(int fiber = bundle_it->second.fiber_min;
	fiber <= bundle_it->second.fiber_max; ++fiber) {
      
      const specex::Trace& trace = psf->GetTrace(fiber);
      double trace_wave_min = trace.X_vs_W.xmin;
      double trace_wave_max = trace.X_vs_W.xmax;
      double min_wave = 3600; //wave_step*int(max(3500.,trace_wave_min)/wave_step+1);
      
      SPECEX_INFO("bundle " << bundle_it->first << " fiber " << fiber << " wave range = " << trace_wave_min << " " << trace_wave_max);
      SPECEX_INFO("bundle " << bundle_it->first << " fiber " << fiber << " coefs = " << trace.X_vs_W.coeff);
      
      

      // loop on wavelength
      for(double wave = min_wave ; wave<= trace_wave_max; wave += wave_step) {
	
	specex::Spot_p spot(new specex::Spot());
	
	spot->xc = trace.X_vs_W.Value(wave);
	spot->yc = trace.Y_vs_W.Value(wave);
	spot->flux = 1;
	spot->wavelength = wave;
	spot->fiber = fiber;
	spot->fiber_bundle = bundle_it->first;
	spots.push_back(spot);
	SPECEX_INFO("fiber=" << spot->fiber << " wave=" << spot->wavelength << " x=" << spot->xc << " y=" << spot->yc);

	//break; // DEBUG, only one spot per fiber
      }
      
      //exit(12); // DEBUG
    }
  }
  
#ifdef CONTINUUM
  SPECEX_INFO("remove continuum");
  for(std::map<int,specex::PSF_Params>::iterator bundle_it = psf->ParamsOfBundles.begin();
      bundle_it != psf->ParamsOfBundles.end(); ++bundle_it) {
    bundle_it->second.ContinuumPol.coeff.clear();
  }
#endif

  if(0) {
    SPECEX_INFO("remove psf tail");
    for(std::map<int,specex::PSF_Params>::iterator bundle_it = psf->ParamsOfBundles.begin();
	bundle_it != psf->ParamsOfBundles.end(); ++bundle_it) {
      bundle_it->second.AllParPolXW[psf->ParamIndex("TAILAMP")]->coeff.clear();
    }
  }

   // create image
  specex::image_data img(psf->ccd_image_n_cols,psf->ccd_image_n_rows);
  
  specex::image_data weight(psf->ccd_image_n_cols,psf->ccd_image_n_rows);
  for(size_t i=0;i<weight.data.size();i++) weight.data(i)=1;
  
  bool only_on_spots = false;
  bool only_psf_core = false;
  bool only_positive = false;

  if(output_image_fits_filename!="") {
    SPECEX_INFO("computing  image ...");  
    specex::parallelized_compute_model_image(img,weight,psf,spots,only_on_spots,only_psf_core,only_positive,0,0);
    //specex::compute_model_image(img,weight,psf,spots,only_on_spots,only_psf_core,only_positive,-1,-1,12);
  
    SPECEX_INFO("writing image ...");
    fitsfile * fp;  
    harp::fits::create ( fp, output_image_fits_filename );
    harp::fits::img_append < double > ( fp, img.n_rows(), img.n_cols() );
    harp::fits::img_write ( fp, img.data, false );
    
    int status=0;
    char comment[800];
    fits_write_comment(fp,"generateds by specex_test_psf_specter_io",&status);
    sprintf(comment,"from %s",input_psf_xml_filename.c_str());
    fits_write_comment(fp,comment,&status);
    sprintf(comment,"spots generated every %fA, starting at %fA",wave_step,spots[0]->wavelength);
    fits_write_comment(fp,comment,&status);
    fits_write_comment(fp,"flux in each spot is 1 electron",&status);
    harp::fits::close ( fp );
    


  }else{
    SPECEX_INFO("no full image");  
  }
  
  { // one spot
    vector<specex::Spot_p> one_spot;
    specex::Spot_p spot = spots[0];
    spot->wavelength = (double)((int)spot->wavelength); 
    one_spot.push_back(spot);
    
    img    = specex::image_data(spot->xc+50,spot->yc+50);
    weight = specex::image_data(spot->xc+50,spot->yc+50);
    for(size_t i=0;i<weight.data.size();i++) weight.data(i)=1;
    
    char filename[1000];
    
    specex::PSF_Params& params_of_bundle = psf->ParamsOfBundles.find(spot->fiber_bundle)->second;
    
    /*
    if(psf->HasParam("GHNSIG")) {
      params_of_bundle.AllParPolXW[psf->ParamIndex("GHNSIG")]->coeff.clear();
      params_of_bundle.AllParPolXW[psf->ParamIndex("GHNSIG")]->coeff(0)=1000; // don't want to cut off core
      cout << "set GHNSIG to very large value" << endl;
    }
    */
    
    harp::vector_double psf_params = psf->AllLocalParamsFW(spot->fiber,spot->wavelength,spot->fiber_bundle);
    
    double scal2    = 0;
    if(psf->Name() ==  "GaussHermite2PSF") {
      scal2    = psf_params[psf->ParamIndex("GH2-0-0")];
      cout << "norm of first  GH PSF = " << (1-scal2) << endl;
      cout << "norm of second GH PSF = " <<  scal2 << endl;
    }
    
    vector<string> names = psf->DefaultParamNames();
    for(size_t k=0;k<psf_params.size(); k++) {
      cout << k << " " << names[k] << " " << psf_params(k) << endl;
    }
    
    img.data.clear();
    specex::parallelized_compute_model_image(img,weight,psf,one_spot,only_on_spots,only_psf_core,only_positive,50,50,spot->fiber_bundle);
    specex::image_data all_contributions_img = img;
    
    
    params_of_bundle.AllParPolXW[psf->ParamIndex("TAILAMP")]->coeff.clear();
    img.data.clear();
    specex::parallelized_compute_model_image(img,weight,psf,one_spot,only_on_spots,only_psf_core,only_positive,50,50,spot->fiber_bundle);
    specex::image_data without_tails_img = img;
    specex::image_data only_tails_img = all_contributions_img;
    only_tails_img.data -= without_tails_img.data;
    
    for(size_t p=0;p<params_of_bundle.AllParPolXW.size();p++) {
      if(names[p].find("GH2-")!=names[p].npos) {
	cout << "clearing " << names[p] << endl;
	params_of_bundle.AllParPolXW[p]->coeff.clear();
      }
    }
    
    sprintf(filename,"single_spot_fiber_%d_wave_%d.fits",spot->fiber,(int)spot->wavelength);
    specex::write_new_fits_image(filename,all_contributions_img);
    sprintf(filename,"single_spot_fiber_%d_wave_%d_without_tails.fits",spot->fiber,(int)spot->wavelength);
    specex::write_new_fits_image(filename,without_tails_img);

    //sprintf(filename,"single_spot_fiber_%d_wave_%d_only_tails.fits",spot->fiber,(int)spot->wavelength);
    //specex::write_new_fits_image(filename,only_tails_img);
    
    if(psf->Name() ==  "GaussHermite2PSF") {
      img.data.clear();
      specex::parallelized_compute_model_image(img,weight,psf,one_spot,only_on_spots,only_psf_core,only_positive,50,50,spot->fiber_bundle);
      specex::image_data only_core_gaussian_img = img;
    
      // need to compute naked gaussian to fix amplitude offset
      for(size_t p=0;p<params_of_bundle.AllParPolXW.size();p++) {
	if(names[p].find("GH-")!=names[p].npos) {
	  cout << "clearing " << names[p] << endl;
	  params_of_bundle.AllParPolXW[p]->coeff.clear();
	}
      }
      
      img.data.clear();
      specex::parallelized_compute_model_image(img,weight,psf,one_spot,only_on_spots,only_psf_core,only_positive,50,50,spot->fiber_bundle);
      specex::image_data only_core_naked_gaussian_img = img;
      only_core_gaussian_img.data -= scal2*only_core_naked_gaussian_img.data;
      
      
      
      specex::image_data only_second_gaussian_img = without_tails_img;
      only_second_gaussian_img.data -= only_core_gaussian_img.data;
      
    
      specex::image_data zero = all_contributions_img;
      zero.data -= only_core_gaussian_img.data;
      zero.data -= only_second_gaussian_img.data;
      zero.data -= only_tails_img.data;
      
      sprintf(filename,"single_spot_fiber_%d_wave_%d_only_core_gaussian.fits",spot->fiber,(int)spot->wavelength);
      specex::write_new_fits_image(filename,only_core_gaussian_img);
      sprintf(filename,"single_spot_fiber_%d_wave_%d_only_second_gaussian.fits",spot->fiber,(int)spot->wavelength);
      specex::write_new_fits_image(filename,only_second_gaussian_img);
      
      sprintf(filename,"zero.fits");
      specex::write_new_fits_image(filename,zero);
    
    
    
    }
  } // end of one spot

  return EXIT_SUCCESS;
}
