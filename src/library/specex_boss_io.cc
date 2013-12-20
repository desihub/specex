#include <iostream>

#include "harp.hpp"

#include "specex_trace.h"

#include "specex_fits.h"
#include "specex_boss_io.h"
#include "specex_spectrograph.h"
#include "specex_message.h"


using namespace std;


void specex::read_BOSS_singleset_in_fits(specex::TraceSet& traceset,
					 const std::string& filename, int hdu_number, specex::TraceSetType ttype,
					 bool verbose) {
 
  
  specex::BOSS_Spectrograph spectro;
  int BOSS_NUMBER_OF_FIBERS_PER_CCD = spectro.number_of_fiber_bundles_per_ccd*spectro.number_of_fibers_per_bundle;

  specex::FitsTable table(filename,hdu_number,verbose);
  
  // test formatting
  if(table.data.empty()) SPECEX_ERROR("BOSS IO error : empty data in table");

  double yjumplo=-12;
  double yjumphi=-12;
  double yjumpval=-12;
  
  std::vector<specex::FitsTableEntry>& data = table.data[0];
  int data_size = data.size();
  if(! (data_size == 4 || data_size == 7)) {cerr << "ERROR specex::Trace::read_sdss_fits data has " << data.size() << " rows and I expect 4 or 7" << endl; HARP_THROW("BOSS IO error"); }
  if(data[0].string_val != "legendre") {cerr << "ERROR specex::Trace::read_sdss_fits func is '" << data[0].string_val << "' and only legendre is implemented" << endl; HARP_THROW("BOSS IO error");}
  if(data[1].double_vals.size()  != 1 ) {cerr << "ERROR specex::Trace::read_sdss_fits col. 2 has wrong size (not 1)" << endl; HARP_THROW("BOSS IO error");}
  if(data[2].double_vals.size()  != 1 ) {cerr << "ERROR specex::Trace::read_sdss_fits col. 3 has wrong size (not 1)" << endl; HARP_THROW("BOSS IO error");}
  if(data_size==7) {
    if(data[4].double_vals.size()  != 1 ) {cerr << "ERROR specex::Trace::read_sdss_fits col. 5 has wrong size (not 1)" << endl; HARP_THROW("BOSS IO error");}
    if(data[5].double_vals.size()  != 1 ) {cerr << "ERROR specex::Trace::read_sdss_fits col. 6 has wrong size (not 1)" << endl; HARP_THROW("BOSS IO error");}
    if(data[6].double_vals.size()  != 1 ) {cerr << "ERROR specex::Trace::read_sdss_fits col. 7 has wrong size (not 1)" << endl; HARP_THROW("BOSS IO error");}
  
    
    
    for(int i=4;i<7;i++) {
      const string& colname = table.columns[i].name;
      specex::FitsTableEntry& entry = data[i];
      if(colname == "XJUMPLO") {
	yjumplo = entry.double_vals(0);
      }else if(colname == "XJUMPHI") {
	yjumphi = entry.double_vals(0); 
      }else if(colname == "XJUMPVAL") {
	yjumpval = entry.double_vals(0);
      }else{
	SPECEX_ERROR(" specex::Trace::read_sdss_fits dont' know what to do with column named '" << colname << "'");
      }
    }
    SPECEX_INFO(" specex::Trace::read_sdss_fits, has y-jump : lo,hi,val = " << yjumplo << "," << yjumphi << "," << yjumpval);
    if(yjumplo==-12 || yjumphi==-12 || yjumpval==-12) {
      SPECEX_ERROR("specex::Trace::read_sdss_fits, missing something");
    }
  }


  int ntotcoeff = data[3].double_vals.size();
  int ncoeff_per_fiber = ntotcoeff/NUMBER_OF_FIBERS_PER_CCD;
  if(ntotcoeff%NUMBER_OF_FIBERS_PER_CCD!=0) {
    cerr << "ERROR specex::Trace::read_sdss_fits col. 3 has bizarre size " << ntotcoeff << endl; 
    HARP_THROW("BOSS IO error");
  }
  
  
  if(traceset.size()!=NUMBER_OF_FIBERS_PER_CCD) traceset.resize(NUMBER_OF_FIBERS_PER_CCD);
  
  const double& xmin = data[1].double_vals(0);
  const double& xmax = data[2].double_vals(0);
  const harp::vector_double & coefficients = data[3].double_vals;
  
  int index=0;
  for(int fiber=0;fiber<NUMBER_OF_FIBERS_PER_CCD; fiber++) {

    specex::Trace& trace = traceset[fiber];
    trace.fiber=fiber;
    
    trace.yjumplo = yjumplo;
    trace.yjumphi = yjumphi;
    trace.yjumpval = yjumpval;
    
    specex::Legendre1DPol* pol = 0;
    switch(ttype) {
    case specex::WY : pol = &trace.lW_vs_Y; break;
    case  specex::XY : pol = &trace.X_vs_Y; break;
    case  specex::YW : pol = &trace.Y_vs_lW; break;
    case  specex::XW : pol = &trace.X_vs_lW; break;
    default : {cerr << "unknown TraceSetType " << ttype << endl; HARP_THROW("BOSS IO error"); };
    }
    
    pol->xmin = xmin;
    pol->xmax = xmax;
    pol->deg = ncoeff_per_fiber-1;
    pol->coeff.resize(ncoeff_per_fiber);
      
    for(int i=0;i<ncoeff_per_fiber;i++,index++) 
      pol->coeff(i) = coefficients(index);
  }
  
  if(verbose) cout <<"INFO specex::TraceSet::ReadSDSS_SingleSet_Fits successful" << endl;
  
}



void specex::read_BOSS_traceset_in_fits(
					specex::TraceSet& traceset,
					const std::string& wave_vs_y_filename, int wave_vs_y_hdu_number,
					const std::string& x_vs_y_filename, int x_vs_y_hdu_number,
					bool verbose
					) {
  
  if(verbose) cout << "INFO Spec2DTraceSet::ReadSDSS_FullSet_Fits starting" << endl;
  if(verbose) cout << "INFO Reading " << wave_vs_y_filename << "[" << wave_vs_y_hdu_number << "]" << endl;
  if(verbose) cout << "INFO Reading " << x_vs_y_filename << "[" << x_vs_y_hdu_number << "]" << endl;
  
  
  specex::read_BOSS_singleset_in_fits(traceset,wave_vs_y_filename,wave_vs_y_hdu_number,specex::WY,verbose);
  specex::read_BOSS_singleset_in_fits(traceset,x_vs_y_filename,x_vs_y_hdu_number,specex::XY,verbose);
  
  for(int fiber=0;fiber<traceset.size();fiber++) {
    
    specex::Trace& trace = traceset[fiber];
    // first need to invert lW_vs_Y
    trace.Y_vs_lW = trace.lW_vs_Y.Invert(2); 
    
    // now compose
    trace.X_vs_lW = composed_pol(trace.X_vs_Y,trace.Y_vs_lW);   
  }
  if(verbose) cout << "INFO specex::read_BOSS_traceset_in_fits done" << endl;
  

}

void specex::read_BOSS_keywords(PSF_p psf, const std::string& arc_image_filename) {
  
  fitsfile * fp= 0;
  harp::fits::open_read (fp,arc_image_filename);
  try{ harp::fits::key_read (fp,"EXPOSURE",psf->arc_exposure_id); }catch(...) { SPECEX_WARNING("could not read EXPOSURE key in " << arc_image_filename);}
  try{ harp::fits::key_read (fp,"PLATEID",psf->plate_id); }catch(...) { SPECEX_WARNING("could not read PLATEID key in " << arc_image_filename);}
  try{ harp::fits::key_read (fp,"MJD",psf->mjd); }catch(...) { SPECEX_WARNING("could not read MJD key in " << arc_image_filename);}
  try{ harp::fits::key_read (fp,"CAMERAS",psf->camera_id); }catch(...) { SPECEX_WARNING("could not read CAMERAS key in " << arc_image_filename);}
  harp::fits::close(fp);
}
