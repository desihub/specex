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
#include <specex_brent.h>

using namespace std;

namespace popts = boost::program_options;

#define _GNU_SOURCE 1
#ifndef __USE_GNU
#define __USE_GNU
#endif
#include <fenv.h>

namespace specex {
  
  class SpotC : public Spot {
  public :
    double cont; // continuum
    SpotC() : Spot() {
      cont = 0;
    }
  };
  BOOST_SERIALIZATION_SHARED_PTR(SpotC)  
  typedef boost::shared_ptr < specex::SpotC > SpotC_p;
  typedef boost::weak_ptr < specex::SpotC > SpotC_wp;
};



specex::PSF_p psf;
specex::image_data image,weight;

specex::Legendre2DPol dxpol;
specex::Legendre2DPol dypol;
harp::matrix_double A;
harp::vector_double B;
harp::vector_double H;
harp::vector_double pos_der;
harp::vector_double chi2pdf_of_spots;
bool include_tails = true;


struct BrentBox { // handler to data and params needed to compute chi2 in brent routine
  vector<specex::SpotC_p> spots_to_fit;
  bool fit_dx,fit_flux,fit_cont;
  double chi2_0;
  BrentBox(vector<specex::SpotC_p> i_spots_to_fit, bool i_fit_dx, bool i_fit_flux, bool i_fit_cont, const double& i_chi2_0) :
    spots_to_fit(i_spots_to_fit),
    fit_dx(i_fit_dx),
    fit_flux(i_fit_flux),
    fit_cont(i_fit_cont),
    chi2_0(i_chi2_0)
  {}
};

double computechi2ab(vector<specex::SpotC_p> spots_to_fit, bool fit_flux, bool fit_cont, bool fit_dx);
double compute_chi2_for_a_given_step(const double &current_step, BrentBox* bbox) {

  int npardx = 0;
  int npardy = 0;
  int first_flux=0;
  int first_cont=0;
  int npar = 0;
  if(bbox->fit_dx) { 
    npardx=dxpol.coeff.size(); 
    npardy=dypol.coeff.size(); 
    npar += (npardx+npardy);
  }
  if(bbox->fit_flux) {
    first_flux = npar;
    npar += bbox->spots_to_fit.size();
  }
  if(bbox->fit_cont) {
    first_cont = npar;
    npar += bbox->spots_to_fit.size();
  }
    

  harp::vector_double sB=current_step*B;
  if(bbox->fit_dx) {
    ublas::noalias(dxpol.coeff) += ublas::project(sB,ublas::range(0,npardx));
    ublas::noalias(dypol.coeff) += ublas::project(sB,ublas::range(npardx,npardx+npardy));
  }
  if(bbox->fit_flux) {
    for(size_t s=0;s<bbox->spots_to_fit.size();s++) {
      bbox->spots_to_fit[s]->flux += sB(first_flux+s);
    }
  }
  if(bbox->fit_cont) {
    for(size_t s=0;s<bbox->spots_to_fit.size();s++) {
      bbox->spots_to_fit[s]->cont += sB(first_cont+s);
    }
  }

  double chi2=computechi2ab(bbox->spots_to_fit,false,false,false);
  SPECEX_INFO("brent step = " << current_step << " chi2= " << chi2 << " dchi2=" << bbox->chi2_0-chi2);
  // rewind
  if(bbox->fit_dx) {
    ublas::noalias(dxpol.coeff) -= ublas::project(sB,ublas::range(0,npardx));
    ublas::noalias(dypol.coeff) -= ublas::project(sB,ublas::range(npardx,npardx+npardy));
  }
  if(bbox->fit_flux) {
    for(size_t s=0;s<bbox->spots_to_fit.size();s++) {
      bbox->spots_to_fit[s]->flux -= sB(first_flux+s);
    }
  }
  if(bbox->fit_cont) {
    for(size_t s=0;s<bbox->spots_to_fit.size();s++) {
      bbox->spots_to_fit[s]->cont -= sB(first_cont+s);
    }
  }

  return chi2;
  
}

void robustify(vector<specex::SpotC_p>& spots) {
  SPECEX_INFO("robustify");

  if(chi2pdf_of_spots.size() != spots.size()) SPECEX_ERROR("spots and chi2pdf_of_spots don't have same size");

  int n_discarded=0;

  double mean_chi2pdf=0;
  double rms_chi2pdf=0;
  int nspots=0;
  for(size_t s=0;s<spots.size();s++) {
    if(chi2pdf_of_spots(s)>0) {
      mean_chi2pdf += chi2pdf_of_spots(s);
      rms_chi2pdf += chi2pdf_of_spots(s)*chi2pdf_of_spots(s);
      nspots++;
    }
  }
  if(nspots==0)  SPECEX_ERROR("0 spots with pos. chi2pdf");
  mean_chi2pdf/=nspots;
  rms_chi2pdf=sqrt(rms_chi2pdf/nspots-mean_chi2pdf);
  SPECEX_INFO("mean chi2/ndf = " << mean_chi2pdf << " rms=" << rms_chi2pdf);
  double nsig=4;
  
    
  vector<specex::SpotC_p> tmp_spots;
  for(size_t s=0;s<spots.size();s++) {
    if(chi2pdf_of_spots(s)>mean_chi2pdf+nsig*rms_chi2pdf || spots[s]->flux<=0) { 
      SPECEX_INFO("discarding sky line at fiber " << spots[s]->fiber << " " << spots[s]->wavelength << " because chi2pdf=" << chi2pdf_of_spots(s) << " flux=" << spots[s]->flux);
      continue;
    }
    tmp_spots.push_back(spots[s]);
  }
  n_discarded = spots.size()-tmp_spots.size();
  spots = tmp_spots;
  

  SPECEX_INFO("keeping " << spots.size() << " sky lines, " << n_discarded << " discarded");  
  
}

double computechi2ab(vector<specex::SpotC_p> spots_to_fit, bool fit_flux, bool fit_cont, bool fit_dx) {

  
  
  int npardx = 0;
  int npardy = 0;
  int npar   = 0;
  int first_flux = 0;
  int first_cont = 0;
  if(fit_dx) { npardx=dxpol.coeff.size(); npardy=dypol.coeff.size(); npar += (npardx+npardy);}
  if(fit_flux) { first_flux = npar ; npar += spots_to_fit.size();}
  if(fit_cont) { first_cont = npar ; npar += spots_to_fit.size();}
  
  
  
  //cout << " B.size() H.size() npar = " << B.size() << " " << H.size() << " " << npar << endl;
  //cout << "fit_dx ... = " << fit_dx << " " << fit_flux << " " << fit_cont << endl;
  //cout << "spots_to_fit.size()  = " << spots_to_fit.size() << endl;
  
  
  if(fit_dx || fit_flux || fit_cont) {
    if(int(B.size()) != npar) {
      A.resize(npar,npar);
      B.resize(npar);
      H.resize(npar);
    }
    A.clear();
    B.clear();
    if(fit_dx) pos_der.resize(2);
    //cont_pos_der.resize(2);
  }
  
  double chi2=0;
  harp::vector_double contprof(2*psf->hSizeX+1);
  if(chi2pdf_of_spots.size() != spots_to_fit.size()) chi2pdf_of_spots.resize(spots_to_fit.size());
  for(size_t s=0;s<spots_to_fit.size();s++) {
    
    specex::SpotC_p spot = spots_to_fit[s];
    double& spot_chi2 = chi2pdf_of_spots(s);
    spot_chi2=0;

    int begin_i = max(0,int(floor(spot->xc))-psf->hSizeX);
    int end_i   = min(int(image.n_cols()),int(floor(spot->xc))+psf->hSizeX+1);
    int begin_j = max(0,int(floor(spot->yc))-psf->hSizeY);
    int end_j   = min(int(image.n_rows()),int(floor(spot->yc))+psf->hSizeY+1);
    
    harp::vector_double psf_params   = psf->AllLocalParamsXW(spot->xc,spot->wavelength,spot->fiber_bundle);
    harp::vector_double dx_monomials,dy_monomials;
    if(fit_dx) {
      dx_monomials = dxpol.Monomials(spot->xc,spot->wavelength);
      dy_monomials = dypol.Monomials(spot->xc,spot->wavelength);
    }

    double x = spot->xc+dxpol.Value(spot->xc,spot->wavelength);
    double y = spot->yc+dypol.Value(spot->xc,spot->wavelength);
    
    /* compute a cross-dispersion profile of PSF */
    contprof.clear();
    for(int i=begin_i;i<end_i;i++) {
      contprof(i-begin_i)=psf->PSFValueWithParamsXY(x,y,i,int(y),psf_params,NULL,NULL,true,include_tails);
    }


    int npix=0;
    for(int j=begin_j;j<end_j;j++) {
      for(int i=begin_i;i<end_i;i++) {
	double w=weight(i,j);
	if(w==0) continue;
	double res = image(i,j);
	if(fit_dx || fit_flux || fit_cont) H.clear();
	
	double prof = 0;
	if(fit_dx) {
	  pos_der.clear();
	  prof = psf->PSFValueWithParamsXY(x,y,i,j,psf_params,&pos_der,NULL,true,include_tails);
	}else{
	  prof = psf->PSFValueWithParamsXY(x,y,i,j,psf_params,NULL,NULL,true,include_tails);
	}
	
	{
	  res -= spot->flux*prof;
	  res -= spot->cont*contprof(i-begin_i);
	  
	  if(fit_flux) H(first_flux+s) += prof;
	  if(fit_cont) H(first_cont+s) += contprof(i-begin_i);
	  
	  if(fit_dx) {
	    ublas::noalias(ublas::project(H,ublas::range(0,npardx))) += (pos_der(0) * spot->flux)*dx_monomials;
	    ublas::noalias(ublas::project(H,ublas::range(npardx,npardx+npardy))) += (pos_der(1) * spot->flux)*dy_monomials;
	  } 
	}
	
	double dchi2 = w*res*res;
	chi2 += dchi2;
	spot_chi2 += dchi2;
	npix++;
	

	if(fit_dx || fit_flux) {
	  specex::syr(w,H,A);
	  specex::axpy(w*res,H,B);
	}

      } // end of loop on i
    } // end of loop on j
    
    if(npix>1) {
      spot_chi2 /= (npix-1);
    }else{
      spot_chi2=0;
    }

  } // end of loop on spot
  
  return chi2;
}

double fit(vector<specex::SpotC_p>& spots_to_fit, bool fit_flux, bool fit_cont, bool fit_dx) {
  int npar = 0;
  int npardx=dxpol.coeff.size();
  int npardy=dypol.coeff.size();
  if(fit_dx) npar += (npardx+npardy);
  if(fit_flux) npar += spots_to_fit.size();
  if(fit_cont) npar += spots_to_fit.size();
  
  double previous_chi2 = computechi2ab(spots_to_fit,fit_flux,fit_cont,fit_dx);
  int status = specex::cholesky_solve(A,B);
  if(status!=0) SPECEX_ERROR("cholesky solve failed with status " << status);
  int first_flux=0;
  int first_cont=0;
  if(fit_dx) {
    ublas::noalias(dxpol.coeff) += ublas::project(B,ublas::range(0,npardx));
    ublas::noalias(dypol.coeff) += ublas::project(B,ublas::range(npardx,npardx+npardy));
    first_flux = npardx+npardy;
  }
  if(fit_flux) {
    for(size_t s=0;s<spots_to_fit.size();s++) {
      spots_to_fit[s]->flux += B(first_flux+s);
      first_cont++;
    }
  }
  if(fit_cont) {
    for(size_t s=0;s<spots_to_fit.size();s++) {
      spots_to_fit[s]->cont += B(first_cont+s);
    }
  }
  double chi2  = computechi2ab(spots_to_fit,false,false,false);

  double dchi2 = previous_chi2-chi2;
  if(dchi2<0 && fabs(dchi2)>10 && fit_dx) {
    SPECEX_INFO("STARTING BRENT");
    // rewinding
    if(fit_dx) {
      ublas::noalias(dxpol.coeff) -= ublas::project(B,ublas::range(0,npardx));
      ublas::noalias(dypol.coeff) -= ublas::project(B,ublas::range(npardx,npardx+npardy));
    }
    if(fit_flux) {
      for(size_t s=0;s<spots_to_fit.size();s++) {
	spots_to_fit[s]->flux -= B(first_flux+s);
      }
    }
    if(fit_cont) {
      for(size_t s=0;s<spots_to_fit.size();s++) {
	spots_to_fit[s]->cont -= B(first_cont+s);
      }
    }
    double min_step=-0.05;
    double max_step=1.05;
    double prefered_step=1;
    int status=0;
    double best_chi2=0;
    double brent_precision = 0.01;
    
    BrentBox bbox(spots_to_fit,fit_dx,fit_flux,fit_cont,previous_chi2);
    
    
    double best_step = brent((AnalyticFunction*)(compute_chi2_for_a_given_step),
			     min_step,prefered_step,max_step,
			     brent_precision,&bbox,best_chi2,status,100);
    chi2=best_chi2;
    
    B *= best_step;
    
    if(fit_dx) {
      ublas::noalias(dxpol.coeff) += ublas::project(B,ublas::range(0,npardx));
      ublas::noalias(dypol.coeff) += ublas::project(B,ublas::range(npardx,npardx+npardy));
    }
    if(fit_flux) {
      for(size_t s=0;s<spots_to_fit.size();s++) {
	spots_to_fit[s]->flux += B(first_flux+s);
      }
    }
    if(fit_cont) {
      for(size_t s=0;s<spots_to_fit.size();s++) {
	spots_to_fit[s]->cont += B(first_cont+s);
      }
    } 
  }
  return chi2;
}




int main ( int argc, char *argv[] ) {
  
  // to crash when NaN
  feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
  
  A.resize(0,0);
  B.resize(0);
  

  // default arguments
  // --------------------------------------------
  string input_science_image_filename = "";
  string input_psf_xml_filename = "";
  string input_skyline_filename = "";
  string plugmap_filename = "";
  string output_psf_xml_filename = "";
  string output_psf_fits_filename = "";
  string output_list_filename = "";
  string output_residual_fits_image_filename = "";
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
    ( "list", popts::value<string>( &output_list_filename), "output list of corrections for some wave/fibers" )
    ( "res", popts::value<string>( &output_residual_fits_image_filename), "output residuals for monitoring" )
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
  {
    std::ifstream is(input_psf_xml_filename.c_str());
    boost::archive::xml_iarchive xml_ia ( is );
    xml_ia >> BOOST_SERIALIZATION_NVP(psf);
    is.close();
  }
  // modify PSF params to get accurate position derivatives
  // first copy params to restore them afterwards

  std::map<int,double> GHNSIG_values;
  for(std::map<int,specex::PSF_Params>::iterator it=psf->ParamsOfBundles.begin(); it !=psf->ParamsOfBundles.end(); ++it) {
    specex::PSF_Params& psf_params = it->second;
    for(int p=0;p<psf->LocalNAllPar();p++) {
      const string& name = psf->ParamName(p);
      if(name=="GHNSIG") {
	SPECEX_INFO("Set GHNSIG to 100000");
	GHNSIG_values[it->first]=psf_params.AllParPolXW[p]->coeff(0);
	psf_params.AllParPolXW[p]->coeff(0) = 100000;
	for(size_t i=1;i<psf_params.AllParPolXW[p]->coeff.size();i++)
	  psf_params.AllParPolXW[p]->coeff(i)=0;
      }
    }
  }
      

  int saved_hsizex = psf->hSizeX;
  int saved_hsizey = psf->hSizeY;
  psf->hSizeX=3;
  psf->hSizeY=4;
  
  
  
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
  
  vector<int> selected_fibers;
  {
    ifstream file (plugmap_filename.c_str());
    while(!file.eof()){
      std::string str;
      std::getline(file, str);
      if(str.find("PLUGMAPOBJ")==str.npos) continue;
      bool is_obj = (str.find("OBJECT")!=str.npos);
      if(!is_obj) continue; // not a fiber with a target
      
      bool is_sky = (str.find("SKY")!=str.npos);
      //bool is_star = (str.find("SPECTROPHOTO_STD")!=str.npos);
      bool is_gal = (str.find("GALAXY")!=str.npos);
      bool is_qso = (str.find("QSO")!=str.npos);
      
      //if(!is_sky) continue; // only on sky fibers
      
      if((!is_sky) && (!is_gal) && (!is_qso)) continue; 
      
      //cout << str << endl;
      std::stringstream ss(str);
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
      selected_fibers.push_back(fiber);
    }
    file.close();
  }
  SPECEX_INFO("done");
  cout << "Using fibers ";
  for(size_t f=0;f<selected_fibers.size();f++) cout << " " << selected_fibers[f];
  cout << endl;
  
 

  // allocate spots
  // --------------------------------------------
  vector<specex::SpotC_p> spots;
  spots.clear();

  //for(std::map<int,specex::PSF_Params>::iterator bundle_it = psf->ParamsOfBundles.begin(); bundle_it != psf->ParamsOfBundles.end(); bundle_it++) {
  //const specex::PSF_Params& params_of_bundle = bundle_it->second;
    // for(int fiber=params_of_bundle.fiber_min; fiber<=params_of_bundle.fiber_max; fiber++) {
  
  double wavemin=1e7;
  double wavemax=0;
  
  for(size_t f=0;f<selected_fibers.size();f++) {
    int fiber = selected_fibers[f];
    
    const specex::PSF_Params* bundle_params = 0;
    int bundle = -1;
    for(std::map<int,specex::PSF_Params>::iterator bundle_it = psf->ParamsOfBundles.begin(); bundle_it != psf->ParamsOfBundles.end(); bundle_it++) {
      if(bundle_it->second.fiber_min<=fiber && bundle_it->second.fiber_max>=fiber){
	bundle = bundle_it->second.bundle_id;
	bundle_params = &(bundle_it->second);
	break;
      }
    }
    if(bundle==-1) SPECEX_ERROR("Coudn't find bundle of fiber " << fiber);
    
    // check chi2
    if(bundle_params->ndata_in_core<=0) {
      SPECEX_INFO("Bundle of fiber " << fiber << " has ndata_in_core=" << bundle_params->ndata_in_core);
      continue;
    }
    if(bundle_params->chi2_in_core/bundle_params->ndata_in_core > 20) {
      SPECEX_INFO("Discard fiber " << fiber << " with bundle with core chi2/ndata = " << bundle_params->chi2_in_core/bundle_params->ndata_in_core);
      continue;
    }
   

    const specex::Trace& trace = psf->GetTrace(fiber);
    if(trace.Off()) continue;
    
    double wmin = trace.X_vs_W.xmin;
    double wmax = trace.X_vs_W.xmax;
    wavemin = min(wmin,wavemin);
    wavemax = max(wmax,wavemax);
    

    for(size_t w=0;w<sky_lines.size();w++) {
      
      if(sky_lines[w]<wmin || sky_lines[w]>wmax) continue;
      
      specex::SpotC_p spot = specex::SpotC_p(new specex::SpotC());
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
  read_fits_images(input_science_image_filename,image,weight);
  
  if(1) {
    // add few% of signal in weight to account for PSF errors
    double psf_error=0.01;
    for(size_t i=0;i<weight.data.size();i++) {
      double &w=weight.data(i);
      const double &f=image.data(i);
      if(w>0 && f>(20./sqrt(w))) {
	w = 1./(1./w + psf_error*psf_error*f*f);
      }
    }
  }
  


  
  {
    specex::PSF_Fitter fitter(psf,image,weight);
    fitter.include_signal_in_weight = false;
    
    bool ok = true;
    
    for(std::map<int,specex::PSF_Params>::iterator bundle_it = psf->ParamsOfBundles.begin(); bundle_it != psf->ParamsOfBundles.end(); bundle_it++) {
      const specex::PSF_Params& params_of_bundle = bundle_it->second;
      vector<specex::SpotC_p> selected_spots;
      vector<specex::Spot_p> selected_spots_as_spots;
      for(size_t s=0;s<spots.size();s++) {
	if (spots[s]->fiber_bundle == params_of_bundle.bundle_id) {
	  selected_spots.push_back(spots[s]);
	  selected_spots_as_spots.push_back(spots[s]);
	}
      }
      if(selected_spots.size()==0) continue;
      
      fitter.SelectFiberBundle(params_of_bundle.bundle_id); 
      ok = fitter.FitIndividualSpotFluxes(selected_spots_as_spots);
      
    }
    
    vector<specex::SpotC_p> tmp_spots;
    for(size_t s=0;s<spots.size();s++) {
      if(spots[s]->flux<=0 || spots[s]->eflux<=0 || (spots[s]->flux/spots[s]->eflux<20)) { 
	SPECEX_INFO("discarding sky line at fiber " <<  spots[s]->fiber << " " << spots[s]->wavelength);
	continue;
      }
      tmp_spots.push_back(spots[s]);
    }
    SPECEX_INFO("using " << tmp_spots.size() << " sky lines, " << spots.size()-tmp_spots.size() << " discarded");  
    spots = tmp_spots;
    
  }

  // allocate polynomials
  dxpol=specex::Legendre2DPol(2,0,image.n_cols(),2,wavemin,wavemax);
  dypol=specex::Legendre2DPol(2,0,image.n_cols(),2,wavemin,wavemax);
  


  // start actual fit

  double min_dchi2=max(0.1,0.0001*spots.size()*((2*psf->hSizeX+1)*(2*psf->hSizeY+1)));
  SPECEX_INFO("Min delta chi2 = " << min_dchi2 );
  
  int nmaxiter=20;
  bool need_robustification = true;
  bool fit_flux = false;
  bool fit_cont = false;
  bool fit_dx = false;
  for(int iter=0;iter<nmaxiter;iter++) {
    
    fit_flux = (iter%2==0); 
    fit_cont = fit_flux; 
    fit_dx   = (!fit_flux);
    
    /* TOO SLOW , can't run this
    if(iter>4) {
      fit_flux = true;
      fit_dx   = true;
      fit_cont = false;
    }
    */
    
    double previous_chi2 = computechi2ab(spots,false,false,false);
    
    double chi2=0;
    if(fit_dx) { // full fit
      chi2 = fit(spots,fit_flux,fit_cont,fit_dx);
    }else{ // spot by spot 
      for(size_t s=0;s<spots.size();s++) {
	vector<specex::SpotC_p> one_spot;
	one_spot.push_back(spots[s]);
	chi2 += fit(one_spot,fit_flux,fit_cont,fit_dx);
      }
    }
    double dchi2 = previous_chi2-chi2;
    
    SPECEX_INFO("#" << iter << " chi2= " << chi2 << " dchi2= " << dchi2 << " fit_dx=" << fit_dx << " fit_flux=" << fit_flux << " fit_cont=" << fit_cont << " nspots=" << spots.size() << " chi2/ndata=" << chi2/(spots.size()*(2*psf->hSizeX+1)*(2*psf->hSizeY+1)));

    
    

    bool has_just_robustified = false;

    if(dchi2<0 && fit_dx) {
      SPECEX_INFO("neg. dchi2, rewinding ");
      int running_index = 0;
      
      if(fit_dx) {
	int npardx=dxpol.coeff.size();
	int npardy=dypol.coeff.size();
	ublas::noalias(dxpol.coeff) -= ublas::project(B,ublas::range(0,npardx));
	ublas::noalias(dypol.coeff) -= ublas::project(B,ublas::range(npardx,npardx+npardy));
	running_index += npardx+npardy;
      }
      if(fit_flux) {
	for(size_t s=0;s<spots.size();s++,running_index++) {
	  spots[s]->flux -= B(running_index);
	}
      }
      if(fit_cont) {
	for(size_t s=0;s<spots.size();s++,running_index++) {
	  spots[s]->cont -= B(running_index);
	}
      }
      
      if(!need_robustification) 
	break;
      else {
	robustify(spots);
	need_robustification=false;
      }
    }
    
    if(iter>=3 && fit_dx && fabs(dchi2)<min_dchi2) {
      if(!need_robustification) {
	if(!has_just_robustified)
	  break;
      }else{
	robustify(spots);
	need_robustification=false;
	has_just_robustified=true;
	previous_chi2 = 1e12;
      }
    }
    previous_chi2 = chi2;
    
  } // end of loop on iterations
  
  // now apply thoses shifts to the PSF and dump a correction file

  ofstream *os=0;
  if(output_list_filename != "") {    
    os = new ofstream(output_list_filename.c_str());
    *os << "# fiber : " << endl;
    *os << "# wave : " << endl;
    *os << "# x : new coordinate" << endl;
    *os << "# y : new coordinate" << endl;
    *os << "# dx : correction that has been applied" << endl;
    *os << "# dy : correction that has been applied" << endl;
    *os << "# xt : trace fit of above" << endl;
    *os << "# yt : trace fit of above" << endl;
    *os << "#end " << endl;
  }
  
  double wmin=1e12;
  double wmax=0;
  for(std::map<int,specex::PSF_Params>::iterator bundle_it = psf->ParamsOfBundles.begin(); bundle_it != psf->ParamsOfBundles.end(); bundle_it++) {
    for(int fiber=bundle_it->second.fiber_min; fiber<=bundle_it->second.fiber_max; fiber++) {
      specex::Trace& trace = psf->GetTrace(fiber);
      
      wmin=min(wmin,trace.X_vs_W.xmin);
      wmax=max(wmax,trace.X_vs_W.xmax);
    }
  }
  for(std::map<int,specex::PSF_Params>::iterator bundle_it = psf->ParamsOfBundles.begin(); bundle_it != psf->ParamsOfBundles.end(); bundle_it++) {
    for(int fiber=bundle_it->second.fiber_min; fiber<=bundle_it->second.fiber_max; fiber++) {
      specex::Trace& trace = psf->GetTrace(fiber);
      
      int nwave=(trace.X_vs_W.deg+1)*2;
      harp::vector_double wave(nwave);
      harp::vector_double x(nwave);
      harp::vector_double y(nwave);
      
      for(int i=0;i<nwave;i++) {
	wave(i) = wmin+(i*(wmax-wmin))/(nwave-1);
	x(i)    = trace.X_vs_W.Value(wave(i));
	y(i)    = trace.Y_vs_W.Value(wave(i));
	double dx = dxpol.Value(x(i),wave(i));
	double dy = dypol.Value(x(i),wave(i));
	x(i)   += dx;
	y(i)   += dy;
	
	// if(os) *os << fiber << " " << wave(i) << " " << x(i) << " " << y(i) << " " << dx << " " << dy << endl;
	
      }    
      trace.X_vs_W.Fit(wave,x,0,false);
      trace.X_vs_Y.Fit(y,x,0,false);
      trace.Y_vs_W.Fit(wave,y,0,false);
      trace.W_vs_Y.Fit(y,wave,0,false);
      if(os) {
	for(int i=0;i<nwave;i++) {
	  *os << fiber << " " << wave(i) << " " << x(i) << " " << y(i) << " " << dxpol.Value(x(i),wave(i)) << " " << dypol.Value(x(i),wave(i)) << " " << trace.X_vs_W.Value(wave(i)) << " " << trace.Y_vs_W.Value(wave(i)) << endl;
	}
      }
    }
  }
  if(os) {
    os->close();
    SPECEX_INFO("wrote " << output_list_filename);
  }

  for(std::map<int,specex::PSF_Params>::iterator it=psf->ParamsOfBundles.begin(); it !=psf->ParamsOfBundles.end(); ++it) {
    specex::PSF_Params& psf_params = it->second;
    for(int p=0;p<psf->LocalNAllPar();p++) {
      const string& name = psf->ParamName(p);
      if(name=="GHNSIG") {
	psf_params.AllParPolXW[p]->coeff(0) = GHNSIG_values.find(it->first)->second;
      }
    }
  }
  // psf->ParamsOfBundles = saved_ParamsOfBundles; // this does not work because pointers!
  
  
  // now we want to refit gauss-hermite term of deg1 ...
  // TODO (complicated...)
  



  
  // write PSF
  psf->hSizeX = saved_hsizex;
  psf->hSizeY = saved_hsizey;
  
  
  if(output_psf_xml_filename !="")
    write_psf_xml(psf, output_psf_xml_filename);
  if(output_psf_fits_filename !="" && psf->Name() == "GaussHermite2PSF") {
    write_psf_fits(psf,output_psf_fits_filename);
  }
 

  

  // try to save residuals
  

  if(output_residual_fits_image_filename != "") {
    SPECEX_INFO("computing residuals");
    specex::image_data model_without_tails=image;
    model_without_tails.data *= 0;
    specex::image_data model_with_tails=image;
    model_with_tails.data *= 0;
    
    specex::image_data fibers=image;
    fibers.data *= 0;
    
    psf->hSizeY = 6;
    psf->hSizeX = 3;
    harp::vector_double contprof(2*psf->hSizeX+1);
    for(size_t s=0;s<spots.size();s++) {
      specex::SpotC_p spot = spots[s];

      double x = psf->GetTrace(spot->fiber).X_vs_W.Value(spot->wavelength);
      double y = psf->GetTrace(spot->fiber).Y_vs_W.Value(spot->wavelength);
      //double x = spot->xc+dxpol.Value(spot->xc,spot->wavelength);
      //double y = spot->yc+dypol.Value(spot->xc,spot->wavelength);
      

      int begin_i = max(0,int(floor(x))-psf->hSizeX);
      int end_i   = min(int(image.n_cols()),int(floor(x))+psf->hSizeX+1);
      int begin_j = max(0,int(floor(y))-psf->hSizeY);
      int end_j   = min(int(image.n_rows()),int(floor(y))+psf->hSizeY+1); 
      harp::vector_double psf_params   = psf->AllLocalParamsXW(x,spot->wavelength,spot->fiber_bundle);
      
      /* compute a cross-dispersion profile of PSF */
      contprof.clear();
      for(int i=begin_i;i<end_i;i++) {
	contprof(i-begin_i)=psf->PSFValueWithParamsXY(x,y,i,int(y),psf_params,NULL,NULL,true,include_tails);
      }
      
      for(int j=begin_j;j<end_j;j++) {
	for(int i=begin_i;i<end_i;i++) {
	  double w=weight(i,j);
	  if(w==0) continue;
	  double prof = psf->PSFValueWithParamsXY(x,y,i,j,psf_params,NULL,NULL,true,true); // with tail
	  model_with_tails(i,j) += spot->flux*prof;
	  model_with_tails(i,j) += spot->cont*contprof(i-begin_i);
	  
	  prof = psf->PSFValueWithParamsXY(x,y,i,j,psf_params,NULL,NULL,true,false); // without tail
	  model_without_tails(i,j) += spot->flux*prof;
	  model_without_tails(i,j) += spot->cont*contprof(i-begin_i);
	  fibers(i,j) = spot->fiber;
	}
      }
    }
    specex::image_data residuals=image;
    
    for(size_t i=0;i<residuals.data.size();i++) {
      if(fabs(model_without_tails.data(i))>0) 
	residuals.data(i)=image.data(i)-model_without_tails.data(i);
      else
	residuals.data(i)=0;
    }
    
    specex::image_data pull=residuals;
    
    for(size_t i=0;i<residuals.data.size();i++) {
      if(fabs(residuals.data(i))>0)
	pull.data(i)=residuals.data(i)*sqrt(weight.data(i));
    }
    
    
  
    fitsfile * fp;  
    harp::fits::create ( fp, output_residual_fits_image_filename);
      
    harp::fits::img_append < double > ( fp, image.n_rows(), image.n_cols() );
    harp::fits::img_write ( fp, image.data ,false);
    harp::fits::key_write(fp,"WHAT","DATA","");
    
    harp::fits::img_append < double > ( fp, image.n_rows(), image.n_cols() );
    harp::fits::img_write ( fp, weight.data ,false);
    harp::fits::key_write(fp,"EXTNAME","IVAR","");
    
    harp::fits::img_append < double > ( fp, image.n_rows(), image.n_cols() );
    harp::fits::img_write ( fp, model_with_tails.data ,false);
    harp::fits::key_write(fp,"EXTNAME","MODELWITHTAILS","");
    
    harp::fits::img_append < double > ( fp, image.n_rows(), image.n_cols() );
    harp::fits::img_write ( fp, model_without_tails.data ,false);
    harp::fits::key_write(fp,"EXTNAME","MODELWITHOUTTAILS","");
    
    harp::fits::img_append < double > ( fp, image.n_rows(), image.n_cols() );
    harp::fits::img_write ( fp, residuals.data ,false);
    harp::fits::key_write(fp,"EXTNAME","RESIDUAL","");
    
    harp::fits::img_append < double > ( fp, image.n_rows(), image.n_cols() );
    harp::fits::img_write ( fp, pull.data ,false);
    harp::fits::key_write(fp,"EXTNAME","PULL","");
    
    harp::fits::img_append < double > ( fp, fibers.n_rows(), image.n_cols() );
    harp::fits::img_write ( fp, fibers.data ,false);
    harp::fits::key_write(fp,"EXTNAME","FIBERS","");
    
    harp::fits::close ( fp );
    SPECEX_INFO("wrote residuals in " << output_residual_fits_image_filename);
    
  }
  
 

}
