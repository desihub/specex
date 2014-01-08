#include <iostream>

#include "harp.hpp"

#include "specex_trace.h"

#include "specex_fits.h"
#include "specex_boss_io.h"
#include "specex_spectrograph.h"
#include "specex_message.h"


using namespace std;



std::vector<specex::Legendre1DPol> specex::read_BOSS_singleset_in_fits(const std::string& filename, int hdu_number, bool verbose) {
 
  
  specex::FitsTable table(filename,hdu_number,verbose);
  
  // SPECEX_INFO("EXIT FOR DEBUG"); exit(12);
  
  if(table.data.empty()) SPECEX_ERROR("BOSS IO error : empty data in table");

  

  double yjumplo=-12;
  double yjumphi=-12;
  double yjumpval=-12;
  
  
  std::vector<specex::FitsTableEntry>& data = table.data[0];
  int data_size = data.size();
  if(! (data_size == 4 || data_size == 7)) {cerr << "ERROR specex::Trace::read_sdss_fits data has " << data.size() << " rows and I expect 4 or 7" << endl; HARP_THROW("BOSS IO error"); }
  
  specex::FitsTableEntry& func_entry  = data[table.columns["FUNC"].col];
  specex::FitsTableEntry& xmin_entry  = data[table.columns["XMIN"].col];
  specex::FitsTableEntry& xmax_entry  = data[table.columns["XMAX"].col];
  specex::FitsTableEntry& coeff_entry = data[table.columns["COEFF"].col];
  
  
  if(func_entry.string_val != "legendre")  SPECEX_ERROR("specex::Trace::read_sdss_fits func is '" << func_entry.string_val << "' and only legendre is implemented");
  if(xmin_entry.double_vals.size()  != 1 ) SPECEX_ERROR("specex::Trace::read_sdss_fits xmin col.has wrong size (not 1)");
  if(xmax_entry.double_vals.size()  != 1 ) SPECEX_ERROR("specex::Trace::read_sdss_fits xmax col.has wrong size (not 1)");

  
  if(table.columns.find("XJUMPVAL") != table.columns.end()) {
    SPECEX_WARNING("JUMPS IN TRACES, NOT IMPLEMENTED YET!");
    yjumplo  = data[table.columns["XJUMPLO"].col].double_vals(0);
    yjumphi  = data[table.columns["XJUMPHI"].col].double_vals(0);
    yjumpval = data[table.columns["XJUMPVAL"].col].double_vals(0);
    SPECEX_INFO(" specex::Trace::read_sdss_fits, has y-jump : lo,hi,val = " << yjumplo << "," << yjumphi << "," << yjumpval);
  }
  
  if(table.columns["COEFF"].dimension.size()!=2) SPECEX_ERROR("specex::Trace::read_sdss_fits expects a 2 dim coeff entry");
  
  int ncoeff_per_fiber = table.columns["COEFF"].dimension[0];
  int nfibers          = table.columns["COEFF"].dimension[1];
  
  SPECEX_INFO("number of coefficients " << ncoeff_per_fiber << "x" << nfibers);
    
  const double& xmin = xmin_entry.double_vals(0);
  const double& xmax = xmax_entry.double_vals(0);
  const harp::vector_double & coefficients = coeff_entry.double_vals;
  
  std::vector<specex::Legendre1DPol> pols;
  int index = 0;
  for(int fiber=0;fiber<nfibers; fiber++) {

    specex::Legendre1DPol pol;

#warning needs to add jumps to legendre pol
    //trace.yjumplo  = yjumplo;
    //trace.yjumphi  = yjumphi;
    //trace.yjumpval = yjumpval;
    
    pol.xmin = xmin;
    pol.xmax = xmax;
    pol.deg = ncoeff_per_fiber-1;
    pol.coeff.resize(ncoeff_per_fiber);
          
    for(int i=0;i<ncoeff_per_fiber;i++,index++) 
      pol.coeff(i) = coefficients(index);

    pols.push_back(pol);

  }
  
  if(verbose) cout <<"INFO specex::TraceSet::ReadSDSS_SingleSet_Fits successful" << endl;
  return pols;
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
  
  
  std::vector<specex::Legendre1DPol> lW_vs_Y_pols = specex::read_BOSS_singleset_in_fits(wave_vs_y_filename,wave_vs_y_hdu_number,verbose);
  std::vector<specex::Legendre1DPol> X_vs_Y_pols  = specex::read_BOSS_singleset_in_fits(x_vs_y_filename,x_vs_y_hdu_number,verbose);
  
  if(lW_vs_Y_pols.size() !=  X_vs_Y_pols.size()) SPECEX_ERROR("lW_vs_Y and X_vs_Y dont' have same number of fibers");
  
  traceset.resize(X_vs_Y_pols.size());
  for(int fiber=0;fiber<int(lW_vs_Y_pols.size());fiber++) {
    
    specex::Trace& trace = traceset[fiber];
    trace.fiber = fiber;
    
    // need to convert lW into W
    { 
      specex::Legendre1DPol& lW_vs_Y_pol = lW_vs_Y_pols[fiber];
      int n = 100;
      harp::vector_double w(n);
      harp::vector_double y(n);

      
      for(int i=0;i<n;i++) {
	y(i)= lW_vs_Y_pol.xmin + ((lW_vs_Y_pol.xmax-lW_vs_Y_pol.xmin)*i)/(n-1);
	w(i)=pow(10.,lW_vs_Y_pol.Value(y(i)));
      }
      trace.W_vs_Y.xmin = lW_vs_Y_pol.xmin;
      trace.W_vs_Y.xmax = lW_vs_Y_pol.xmax;
      trace.W_vs_Y.deg  = lW_vs_Y_pol.deg+1; // add 1 to minimise errors
      trace.W_vs_Y.coeff.resize(trace.W_vs_Y.deg+1);     
      trace.W_vs_Y.Fit(y,w,0,false);
      
      trace.Y_vs_W.xmin = pow(10.,lW_vs_Y_pol.Value(lW_vs_Y_pol.xmin));
      trace.Y_vs_W.xmax = pow(10.,lW_vs_Y_pol.Value(lW_vs_Y_pol.xmax));
      trace.Y_vs_W.deg  = lW_vs_Y_pol.deg+1; // add 1 to minimise errors
      trace.Y_vs_W.coeff.resize(trace.Y_vs_W.deg+1);     
      trace.Y_vs_W.Fit(w,y,0,false);
      
    }
    
    trace.X_vs_Y = X_vs_Y_pols[fiber];
    
    // now compose
    trace.X_vs_W = composed_pol(trace.X_vs_Y,trace.Y_vs_W);   
  }
  if(verbose) cout << "INFO specex::read_BOSS_traceset_in_fits done" << endl;
  

}

void specex::read_BOSS_keywords(PSF_p psf, const std::string& arc_image_filename, std::map<std::string,std::string>& infos) {
  
  fitsfile * fp= 0;
  harp::fits::open_read (fp,arc_image_filename);
  try{ harp::fits::key_read (fp,"EXPOSURE",psf->arc_exposure_id); }catch(...) { SPECEX_WARNING("could not read EXPOSURE key in " << arc_image_filename);}
  try{ harp::fits::key_read (fp,"PLATEID",psf->plate_id); }catch(...) { SPECEX_WARNING("could not read PLATEID key in " << arc_image_filename);}
  try{ harp::fits::key_read (fp,"MJD",psf->mjd); }catch(...) { SPECEX_WARNING("could not read MJD key in " << arc_image_filename);}
  try{ harp::fits::key_read (fp,"CAMERAS",psf->camera_id); }catch(...) { SPECEX_WARNING("could not read CAMERAS key in " << arc_image_filename);}

  infos.clear();
  string keys[]={"RDNOISE0","RDNOISE1","RDNOISE2","RDNOISE3","EXPTIME"};
  for(int k=0;k<5;k++) {
    harp::fits::key_read (fp,keys[k],infos[keys[k]]);
  }

  harp::fits::close(fp);
}
