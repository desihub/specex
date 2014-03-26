#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <map>

#include <boost/program_options.hpp>

#include <boost/archive/xml_oarchive.hpp>

#include <harp.hpp>

#include <specex_message.h>
#include <specex_linalg.h>
#include <specex_psf.h>
#include <specex_trace.h>
#include <specex_fits.h>
#include <specex_boss_io.h>
#include <specex_psf_fitter.h>
#include <specex_psf_io.h>
#include <specex_serialization.h>

using namespace std;

namespace popts = boost::program_options;

#define _GNU_SOURCE 1
#ifndef __USE_GNU
#define __USE_GNU
#endif
#include <fenv.h>


int main ( int argc, char *argv[] ) {
  
  // to crash when NaN
  feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
  

  // default arguments
  // --------------------------------------------
  string input_science_image_filename = "";
  string input_psf_xml_filename = "";
  string input_skyline_filename = "";
  string plugmap_filename = "";
  string output_psf_xml_filename = "";
  string output_psf_fits_filename = "";
  
  string spectrograph_name = "BOSS";
  

  
  // reading arguments
  // --------------------------------------------
  popts::options_description desc ( "Allowed Options" );
  desc.add_options()
    ( "help,h", "display usage information" )
    ( "in", popts::value<string>( &input_science_image_filename ), "fiber flat pre-reduced fits image file name (mandatory), ex:  sdProc-b1-00108382.fits" )
    ( "psf", popts::value<string>( &input_psf_xml_filename ), "input psf xml filename" )
    ( "plugmap", popts::value<string>( &plugmap_filename ), "plugmap filename (ex: plPlugMapM-3647-55173-01.par)" )
    ( "lines", popts::value<string>( &input_skyline_filename ), "input sky line filename (skylines.dat)" )
    ( "xml", popts::value<string>( &output_psf_xml_filename), "output psf xml filename" )
    ( "fits", popts::value<string>( &output_psf_fits_filename), "output psf fits filename" )
    ;

  popts::variables_map vm;
  
  try {
    
    popts::store(popts::command_line_parser( argc, argv ).options(desc).run(), vm);
    popts::notify(vm);
    
    if ( ( argc < 2 ) || vm.count( "help" ) || ( ! vm.count( "in" ) )  || ( ! vm.count( "psf" ) ) || ( ! vm.count( "lines" ) ) || ( ! vm.count( "xml" ) ) ) {
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

  specex_set_verbose(true);
  specex_set_dump_core(true);

  // read PSF
  // --------------------------------------------
  specex::PSF_p psf;
  {
    std::ifstream is(input_psf_xml_filename.c_str());
    boost::archive::xml_iarchive xml_ia ( is );
    xml_ia >> BOOST_SERIALIZATION_NVP(psf);
    is.close();
  }
  int saved_hsizex = psf->hSizeX;
  int saved_hsizey = psf->hSizeY;
  psf->hSizeX=3;
  psf->hSizeY=3;
  
  
  
  SPECEX_INFO("read sky lines in " << input_skyline_filename);
  vector<double> sky_lines;
  {
    ifstream file (input_skyline_filename.c_str());
    double wave;
    while(!file.eof()){
      std::string str;
      std::getline(file, str);
      if(str[0] == '#') continue;
      std::stringstream ss(str);
      if(!(ss >> wave)) continue;
      sky_lines.push_back(wave);
    }
    file.close();
  }
  SPECEX_INFO("done");
  for(size_t w=0;w<sky_lines.size();w++) {
    SPECEX_INFO("Using sky line " << sky_lines[w] << " A"); 
  }

  

  char band_id;
  int spectro_id;
  sscanf(psf->camera_id.c_str(),"%c%d",&band_id,&spectro_id);
  SPECEX_INFO("Camera " << band_id << spectro_id);
  SPECEX_INFO("Spectrograph number #" << spectro_id);
  SPECEX_INFO("Reading PlugMap in " << plugmap_filename);
  
  vector<int> sky_fibers;
  {
    ifstream file (plugmap_filename.c_str());
    double wave;
    while(!file.eof()){
      std::string str;
      std::getline(file, str);
      if(str.find("PLUGMAPOBJ")==str.npos) continue;
      if(str.find("SKY")==str.npos) continue;
      //cout << str << endl;
      std::stringstream ss(str);
      int count=0;
      string token;
      for(int c=0;c<24;c++) { 
	if(!(ss >> token)) SPECEX_ERROR("Error in reading " <<plugmap_filename);
      }
      int spectro, boss_fiber_id;
      if(!(ss >> spectro)) SPECEX_ERROR("Error in reading " <<plugmap_filename);
      if(!(ss >> boss_fiber_id)) SPECEX_ERROR("Error in reading " <<plugmap_filename);
      if(spectro != spectro_id) continue;
      int fiber = boss_fiber_id-1;
      if(spectro_id==2) fiber -= 500;
      sky_fibers.push_back(fiber);
    }
    file.close();
  }
  SPECEX_INFO("done");
  cout << "Using sky fibers ";
  for(size_t f=0;f<sky_fibers.size();f++) cout << " " << sky_fibers[f];
  cout << endl;
  
 

  // allocate spots
  // --------------------------------------------
  vector<specex::Spot_p> spots;

  //for(std::map<int,specex::PSF_Params>::iterator bundle_it = psf->ParamsOfBundles.begin(); bundle_it != psf->ParamsOfBundles.end(); bundle_it++) {
  //const specex::PSF_Params& params_of_bundle = bundle_it->second;
    // for(int fiber=params_of_bundle.fiber_min; fiber<=params_of_bundle.fiber_max; fiber++) {

  for(size_t f=0;f<sky_fibers.size();f++) {
    int fiber = sky_fibers[f];
    
    int bundle = -1;
    for(std::map<int,specex::PSF_Params>::iterator bundle_it = psf->ParamsOfBundles.begin(); bundle_it != psf->ParamsOfBundles.end(); bundle_it++) {
      if(bundle_it->second.fiber_min<=fiber && bundle_it->second.fiber_max>=fiber){
	bundle = bundle_it->second.bundle_id;
	break;
      }
    }
    if(bundle==-1) SPECEX_ERROR("Coudn't find bundle of fiber " << fiber);
    
    
    const specex::Trace& trace = psf->GetTrace(fiber);
    if(trace.Off()) continue;
    
    double wmin = trace.X_vs_W.xmin;
    double wmax = trace.X_vs_W.xmax;
    
    for(size_t w=0;w<sky_lines.size();w++) {
      
      if(sky_lines[w]<wmin || sky_lines[w]>wmax) continue;
      
      specex::Spot_p spot = specex::Spot_p(new specex::Spot());
      spot->wavelength = sky_lines[w];
      
      spot->xc         = trace.X_vs_W.Value(spot->wavelength);
      spot->yc         = trace.Y_vs_W.Value(spot->wavelength);
      spot->flux       = 0;
      spot->fiber      = fiber;
      spot->fiber_bundle = bundle;
      spot->eflux = 99;
      spot->initial_xc = spot->xc;
      spot->initial_yc = spot->yc;
      spot->initial_flux = spot->flux;
      spots.push_back(spot);
    }
  }
  
  // read image
  // --------------------------------------------
  specex::image_data image,weight;
  read_fits_images(input_science_image_filename,image,weight);
  
  // load fitter
  // --------------------------------------------
  specex::PSF_Fitter fitter(psf,image,weight);
  fitter.include_signal_in_weight = false;
  
  bool ok = true;
  
  for(std::map<int,specex::PSF_Params>::iterator bundle_it = psf->ParamsOfBundles.begin(); bundle_it != psf->ParamsOfBundles.end(); bundle_it++) {
    const specex::PSF_Params& params_of_bundle = bundle_it->second;
    vector<specex::Spot_p> selected_spots;
    for(size_t s=0;s<spots.size();s++) {
      if (spots[s]->fiber_bundle == params_of_bundle.bundle_id)
	selected_spots.push_back(spots[s]);
    }
    if(selected_spots.size()==0) continue;
    
    fitter.SelectFiberBundle(params_of_bundle.bundle_id); 
    ok = fitter.FitIndividualSpotFluxes(selected_spots);
    
  }
  
  return 0;

  // write PSF
  psf->hSizeX = saved_hsizex;
  psf->hSizeY = saved_hsizey;
  
  if(output_psf_xml_filename !="")
    write_psf_xml(psf, output_psf_xml_filename);
  if(output_psf_fits_filename !="" && psf->Name() == "GaussHermite2PSF") {
    write_psf_fits(psf,output_psf_fits_filename);
  }
 
}
