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



specex::PSF_p psf;
specex::image_data image,weight;
vector<specex::Spot_p> spots;
specex::Legendre2DPol dxpol;
specex::Legendre2DPol dypol;
harp::matrix_double A;
harp::vector_double B;
harp::vector_double chi2pdf_of_spots;

void robustify() {
  SPECEX_INFO("robustify");

  if(chi2pdf_of_spots.size() != spots.size()) SPECEX_ERROR("spots and chi2pdf_of_spots don't have same size");

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
  double nsig=5;
  
  
  vector<specex::Spot_p> tmp_spots;
  for(size_t s=0;s<spots.size();s++) {
    if(chi2pdf_of_spots(s)>mean_chi2pdf+nsig*rms_chi2pdf) { 
      SPECEX_INFO("discarding sky line at fiber " << spots[s]->fiber << " " << spots[s]->wavelength << " because chi2pdf=" << chi2pdf_of_spots(s));
      continue;
    }
    tmp_spots.push_back(spots[s]);
  }

  SPECEX_INFO("keeping " << tmp_spots.size() << " sky lines, " << spots.size()-tmp_spots.size() << " discarded");  
  spots = tmp_spots;
}

double computechi2ab(bool fit_flux, bool fit_dx) {
  
  int npardx = 0;
  int npardy = 0;
  int npar   = 0;
  if(fit_dx) { npardx=dxpol.coeff.size(); npardy=dypol.coeff.size(); npar += (npardx+npardy);}
  if(fit_flux) { npar += spots.size();}
  
  harp::vector_double H;
  harp::vector_double pos_der;
  
  if(fit_dx || fit_flux) {
    A.resize(npar,npar); A.clear();
    B.resize(npar); B.clear();
    H.resize(npar); 
    pos_der.resize(2);
  }
  
  double chi2=0;
  
  if(chi2pdf_of_spots.size() != spots.size()) chi2pdf_of_spots.resize(spots.size());
  for(size_t s=0;s<spots.size();s++) {
    
    specex::Spot_p spot = spots[s];
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
    
    int npix=0;
    for(int j=begin_j;j<end_j;j++) {
      for(int i=begin_i;i<end_i;i++) {
	double w=weight(i,j);
	if(w==0) continue;
	double res = image(i,j);
	if(fit_dx || fit_flux) H.clear();
	
	double prof = 0;
	if(fit_dx) {
	  pos_der.clear();
	  prof = psf->PSFValueWithParamsXY(x,y,i,j,psf_params,&pos_der,NULL,true,true);
	}else{
	  prof = psf->PSFValueWithParamsXY(x,y,i,j,psf_params,NULL,NULL,true,true);
	}
	if(prof) {
	  res -= spot->flux*prof;
	  if(fit_flux) H(npardx+s) += prof;
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
  string output_list_filename = "";
  
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
    while(!file.eof()){
      std::string str;
      std::getline(file, str);
      if(str.find("PLUGMAPOBJ")==str.npos) continue;
      if(str.find("SKY")==str.npos) continue;
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
  spots.clear();

  //for(std::map<int,specex::PSF_Params>::iterator bundle_it = psf->ParamsOfBundles.begin(); bundle_it != psf->ParamsOfBundles.end(); bundle_it++) {
  //const specex::PSF_Params& params_of_bundle = bundle_it->second;
    // for(int fiber=params_of_bundle.fiber_min; fiber<=params_of_bundle.fiber_max; fiber++) {
  
  double wavemin=1e7;
  double wavemax=0;
  
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
    wavemin = min(wmin,wavemin);
    wavemax = max(wmax,wavemax);
    

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
  read_fits_images(input_science_image_filename,image,weight);
  
  
  {
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
    
    vector<specex::Spot_p> tmp_spots;
    for(size_t s=0;s<spots.size();s++) {
      if(spots[s]->flux<=0) { 
	SPECEX_WARNING("discarding sky line at fiber " <<  spots[s]->fiber << " " << spots[s]->wavelength);
	continue;
      }
      tmp_spots.push_back(spots[s]);
    }
    SPECEX_INFO("using " << tmp_spots.size() << " sky lines, " << spots.size()-tmp_spots.size() << " discarded");  
    spots = tmp_spots;
    
  }

  // allocate polynomials
  dxpol=specex::Legendre2DPol(1,0,image.n_cols(),1,wavemin,wavemax);
  dypol=specex::Legendre2DPol(1,0,image.n_cols(),1,wavemin,wavemax);
  


  // start actual fit
  int nmaxiter=20;
  int npardx=dxpol.coeff.size();
  int npardy=dypol.coeff.size();
  
  bool need_robustification = true;
  for(int iter=0;iter<nmaxiter;iter++) {
    
    bool fit_flux = (iter%2==0); 
    bool fit_dx   = !fit_flux;

    int npar = 0;
    if(fit_dx) npar += (npardx+npardy);
    if(fit_flux) npar += spots.size();
    
    SPECEX_INFO("fill matrix iter #" <<iter);    
    double previous_chi2 = computechi2ab(fit_flux,fit_dx);
    for(int p=0;p<npar;p++)
      if(A(p,p)==0) A(p,p)=1.e-12;
    
    SPECEX_INFO("solving iter #" <<iter);
    int status = specex::cholesky_solve(A,B);
    if(status!=0) SPECEX_ERROR("cholesky solve failed with status " << status);
    SPECEX_INFO("solving done");
    
    int first_flux=0;
    if(fit_dx) {
      ublas::noalias(dxpol.coeff) += ublas::project(B,ublas::range(0,npardx));
      ublas::noalias(dypol.coeff) += ublas::project(B,ublas::range(npardx,npardx+npardy));
      first_flux = npardx+npardy;
    }
    if(fit_flux) {
      for(size_t s=0;s<spots.size();s++) {
	spots[s]->flux += B(first_flux+s);
      }
     }
    
    SPECEX_INFO("compute chi2 iter #" <<iter);
    double chi2  = computechi2ab(false,false);
    double dchi2 =  previous_chi2-chi2;
    
    SPECEX_INFO("chi2= " << chi2 << " dchi2= " << dchi2);
    if(dchi2<0) {
      SPECEX_WARNING("neg. dchi2, rewinding ");
      if(fit_dx) {
	ublas::noalias(dxpol.coeff) -= ublas::project(B,ublas::range(0,npardx));
	ublas::noalias(dypol.coeff) -= ublas::project(B,ublas::range(npardx,npardx+npardy));
      }
      if(fit_flux) {
	for(size_t s=0;s<spots.size();s++) {
	  spots[s]->flux -= B(first_flux+s);
	}
      }
      if(!need_robustification) 
	break;
      else {
	robustify();
	need_robustification=false;
      }
    }
    
    if(iter>=3 && fit_dx && fabs(previous_chi2-chi2)<0.1) {
      if(!need_robustification) 
	break;
      else{
	robustify();
	need_robustification=false;
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
    *os << "# dx : " << endl;
    *os << "# dy : " << endl;
    *os << "#end " << endl;
  }
  for(std::map<int,specex::PSF_Params>::iterator bundle_it = psf->ParamsOfBundles.begin(); bundle_it != psf->ParamsOfBundles.end(); bundle_it++) {
    for(int fiber=bundle_it->second.fiber_min; fiber<=bundle_it->second.fiber_max; fiber++) {
      specex::Trace& trace = psf->GetTrace(fiber);
      
      double wmin=trace.X_vs_W.xmin;
      double wmax=trace.X_vs_W.xmax;
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
	
	if(os) *os << fiber << " " << wave(i) << " " << dx << " " << dy << endl;
	
      }    
      trace.X_vs_W.Fit(wave,x,0,false);
      trace.X_vs_Y.Fit(y,x,0,false);
      trace.Y_vs_W.Fit(wave,y,0,false);
      trace.W_vs_Y.Fit(y,wave,0,false);
    }
  }
  if(os) {
    os->close();
    SPECEX_INFO("wrote " << output_list_filename);
  }
  
  
  // write PSF
  psf->hSizeX = saved_hsizex;
  psf->hSizeY = saved_hsizey;
  
  if(output_psf_xml_filename !="")
    write_psf_xml(psf, output_psf_xml_filename);
  if(output_psf_fits_filename !="" && psf->Name() == "GaussHermite2PSF") {
    write_psf_fits(psf,output_psf_fits_filename);
  }
 
}
