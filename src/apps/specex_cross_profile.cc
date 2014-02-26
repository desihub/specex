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
  double wave=6000;
  string psf_xml_filename="";
  string spots_xml_filename="";
  string input_fits_image_filename="";
  string output_filename="";  
  double psf_error=0;
  double readout_noise=2;
  
  // reading arguments
  // --------------------------------------------
  popts::options_description desc ( "Allowed Options" );
  desc.add_options()
    ( "help,h", "display usage information" )
    ( "psf", popts::value<string>( &psf_xml_filename ), "psf xml filename" )
    ( "spots", popts::value<string>( &spots_xml_filename ), "spots xml filename" )
    ( "wave", popts::value<double>( &wave ), "" )
    ( "in", popts::value<string>( &input_fits_image_filename ), " input fits image file name" )
    ( "out", popts::value<string>( &output_filename ), " output file name")
    ( "readout-noise", popts::value<double>( &readout_noise ), " readout noise for pull")
    ( "psf-error", popts::value<double>( &psf_error ), " psf relative error for pull")
    ( "core", "dump core files when harp exception is thrown" )
    ;
  
  popts::variables_map vm;
  
  try {
  
    
    popts::store(popts::command_line_parser( argc, argv ).options(desc).run(), vm);
    popts::notify(vm);
    
    if ( ( argc < 2 ) || vm.count( "help" ) || ( ! vm.count( "psf" ) )  || ( ! vm.count( "spots" ) ) 
	 || ( ! vm.count( "in" ) ) || ( ! vm.count( "wave" ) ) || ( ! vm.count( "out" ) )  ) {
      cerr << endl;
      cerr << desc << endl;
      cerr << "example:" << endl;
      cerr << argv[0] << " --psf psf.xml  --spots spots.xml --in sdProc-b1-00108382.fits --wave 5460 --out cross_profile.list" << endl;
      return EXIT_FAILURE;
    }
  }catch(std::exception e) {
    cerr << "error in arguments" << endl;
    cerr << endl;
    cerr << desc << endl;
    return EXIT_FAILURE;
  }
  try {

    specex_set_verbose(true);
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


    double spot_wave = 0;
    vector<specex::Spot_p> spots;
    {  
      double dist = 1e12;
      for(size_t s=0;s<input_spots.size();s++) {
	specex::Spot_p spot = input_spots[s];
	if(spot->status == 0) continue;
	if(fabs(spot->wavelength -wave)<dist) {spot_wave = spot->wavelength;dist=fabs(spot->wavelength -wave);}
	spots.push_back(spot);
      }
    }
    
    SPECEX_INFO("wavelength of closest spot = " << spot_wave << " for input wave " << wave);
    
    int first_fiber = 10000;
    int last_fiber  = 0;
    
    for(std::map<int,PSF_Params>::const_iterator it = psf->ParamsOfBundles.begin();
	it != psf->ParamsOfBundles.end(); ++it) {
      first_fiber = min(first_fiber,it->second.fiber_min);
      last_fiber  = max(last_fiber,it->second.fiber_max);
    }
    // create output images
    // --------------------------------------------
    
    double begin_x = psf->Xccd(first_fiber,spot_wave)-psf->hSizeX;
    double end_x = psf->Xccd(last_fiber,spot_wave)+psf->hSizeX+1;
    int j = int(psf->Yccd((first_fiber+last_fiber)/2,spot_wave));
    
    SPECEX_INFO("selected wave = " << spot_wave);
    SPECEX_INFO("y ccd = " << j);
    SPECEX_INFO("x ccd range " << begin_x << " " << end_x);
    
    ofstream os(output_filename.c_str());
    os << "# x :" << endl;
    os << "# y :" << endl;
    os << "# w :" << endl;
    os << "# data :" << endl;
    os << "# model : core+tail+cont" << endl;
    os << "# core :" << endl;
    os << "# tail :" << endl;
    os << "# cont :" << endl;
    os << "#end" << endl;
    
    vector<harp::vector_double> psf_params;
    for(size_t s=0;s<spots.size();s++) {
      specex::Spot_p& spot = spots[s];
      psf_params.push_back(psf->AllLocalParamsFW(spot->fiber,spot->wavelength,spot->fiber_bundle));
    }

    

    
    
   
    for(int i=begin_x ; i<end_x; i++) {
 
      os << i << " " << j << " " << spot_wave << " ";
      
      double data = 0;
      double val   = 0;
      double val_core   = 0;
      double val_tail   = 0;
      double val_cont   = 0;
    
      if(weight(i,j)=0) continue;
      data = image(i,j);
      
      

#ifdef CONTINUUM
      for(std::map<int,PSF_Params>::const_iterator it = psf->ParamsOfBundles.begin();
	  it != psf->ParamsOfBundles.end(); ++it) {


	const PSF_Params &params_of_bundle = it->second;
	double expfact_for_continuum=1./(2*M_PI*square(params_of_bundle.continuum_sigma_x));

	
	for(int fiber = it->second.fiber_min ; fiber <= it->second.fiber_max; fiber ++) {
	  double x_center = psf->Xccd(fiber,spot_wave);
	  double continuum_val = params_of_bundle.ContinuumPol.Value(wave)*expfact_for_continuum*exp(-0.5*square((i-x_center)/params_of_bundle.continuum_sigma_x));
	  val += continuum_val;
	  val_cont += continuum_val;
	}
      }
#endif

      
      for(size_t s=0;s<spots.size();s++) {
	specex::Spot_p& spot = spots[s];

	bool in_core = fabs(j-spot->yc)<psf->hSizeY+1;
	double val = spot->flux*psf->PSFValueWithParamsXY(spot->xc, spot->yc,i, j,
							  psf_params[s], 0, 0, in_core, 0);
	  
	val += val;
	val_core += val;
      }
	
      os << " " << data << " " << val << " " << val_core << " " << val_tail << " " << val_cont;
      os << endl;
    }
    
  os.close();
  

  
  // ending
  // --------------------------------------------
  }catch(harp::exception e) {
    cerr << "FATAL ERROR (harp) " << e.what() << endl;
    return EXIT_FAILURE;
  }
    /*
      }catch(std::exception e) {
    cerr << "FATAL ERROR " << e.what() << endl;
    return EXIT_FAILURE;
    }*/

  return EXIT_SUCCESS;
}



