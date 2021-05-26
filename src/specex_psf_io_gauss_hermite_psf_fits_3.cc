#include <fstream>
#include <string>
#include <unhrp.h>
#include <specex_fits.h>
#include <specex_trace.h>
#include <specex_gauss_hermite_psf.h>
#include <specex_string.h>

using namespace std ;

void read_gauss_hermite_psf_fits_version_3(specex::PSF_p& psf, fitsfile* fp, int hdu) {
  
  SPECEX_DEBUG("read_gauss_hermite_psf_fits_version_3 hdu=" << hdu);

  int status = 0;
  fits_movabs_hdu ( fp, hdu, NULL, &status ); harp::fits::check ( status );
  
  int GHDEGX;  harp::fits::key_read(fp,"GHDEGX",GHDEGX);
  int GHDEGY;  harp::fits::key_read(fp,"GHDEGY",GHDEGY);
  if(GHDEGX != GHDEGY) {
    SPECEX_ERROR("expect GHDEGX=GHDEGY");
  }

  // psf is now already created
  // psf = specex::PSF_p(new specex::GaussHermitePSF(GHDEGX));
  ((specex::GaussHermitePSF*)(&(*psf)))->SetDegree(GHDEGX);

  // read header keywords
  harp::fits::key_read(fp,"MJD",psf->mjd);
  harp::fits::key_read(fp,"PLATEID",psf->plate_id);
  harp::fits::key_read(fp,"CAMERA",psf->camera_id);
  harp::fits::key_read(fp,"ARCEXP",psf->arc_exposure_id);
  int nx; harp::fits::key_read(fp,"NPIX_X",nx); psf->ccd_image_n_cols = nx;
  int ny; harp::fits::key_read(fp,"NPIX_Y",ny); psf->ccd_image_n_rows = ny;
  harp::fits::key_read(fp,"HSIZEX",psf->hSizeX);
  harp::fits::key_read(fp,"HSIZEY",psf->hSizeY);
  int FIBERMIN; harp::fits::key_read(fp,"FIBERMIN",FIBERMIN);
  int FIBERMAX; harp::fits::key_read(fp,"FIBERMAX",FIBERMAX);
  double WAVEMIN; harp::fits::key_read(fp,"WAVEMIN",WAVEMIN);
  double WAVEMAX; harp::fits::key_read(fp,"WAVEMAX",WAVEMAX);
  int LEGDEG; harp::fits::key_read(fp,"LEGDEG",LEGDEG);
  

  // read table
  specex::FitsTable table;
  table.Read(fp); 
  int param_col = table.columns["PARAM"].col;
  int coeff_col = table.columns["COEFF"].col;
  int degx_col = table.columns["LEGDEGX"].col;
  int degw_col = table.columns["LEGDEGW"].col;
  
  std::vector<std::string> params;
  std::map<std::string,int> param_row;
  std::map<std::string,unhrp::vector_double > param_coeff;
  std::map<std::string,int > param_degx;
  std::map<std::string,int > param_degw;
  
  for(int i=0;i<table.data.size();i++) { 
    std::string pname=table.data[i][param_col].string_val;
    trim(pname);
    params.push_back(pname);
    param_row[pname]=i;
    param_coeff[pname]=table.data[i][coeff_col].double_vals;    
    param_degx[pname]=table.data[i][degx_col].int_vals[0];
    param_degw[pname]=table.data[i][degw_col].int_vals[0];
    SPECEX_DEBUG("read_gauss_hermite_psf " << i << " '" << pname << "' degx=" << param_degx[pname] << " degw=" << param_degw[pname]);
  }
  
  // find bundles
  vector<int> bundles;  
  for(int fiber=FIBERMIN;fiber<=FIBERMAX;fiber++) {
    int bundle = int(param_coeff["BUNDLE"][(LEGDEG+1)*fiber]);
    if(bundle<0) continue;
    if(std::find(bundles.begin(), bundles.end(), bundle) == bundles.end()) bundles.push_back(bundle);
  }
  SPECEX_DEBUG("Number of bundles = " << bundles.size());

  // loop on bundles and fill PSF parameters
  for(int b=0;b<bundles.size();b++) {
    int bundle=bundles[b];

    // here we can test the input bundle requirement
    
    int bundle_fibermin=10000;
    int bundle_fibermax=0;    
    vector<int> fibers_in_bundle;
    for(int fiber=FIBERMIN;fiber<=FIBERMAX;fiber++) {
      int fiber_bundle = int(param_coeff["BUNDLE"][(LEGDEG+1)*fiber]);
      if(fiber_bundle==bundle) {
	fibers_in_bundle.push_back(fiber);
	bundle_fibermin=min(fiber,bundle_fibermin);
	bundle_fibermax=max(fiber,bundle_fibermax);	
      }
    }
    if(fibers_in_bundle.size()==0) {
      SPECEX_WARNING("No fiber in bundle " << bundle << " ???");
      continue;
    }
    
    specex::PSF_Params bundle_params;
    bundle_params.bundle_id=bundle;
    bundle_params.fiber_min=bundle_fibermin;
    bundle_params.fiber_max=bundle_fibermax;
    
    double xmin=100000;
    double xmax=0;
    
    for(int fiber=bundle_fibermin;fiber<=bundle_fibermax;fiber++) {
      const specex::Trace& trace = psf->FiberTraces[fiber];
      for(double wave=WAVEMIN;wave<=WAVEMAX;wave+=10) {
	double x = trace.X_vs_W.Value(wave);
	xmin=min(xmin,x);
	xmax=max(xmax,x);	  
      }
    }
    SPECEX_DEBUG("bundle=" << bundle << " xmin xmax : " << xmin << " " << xmax);
    
#ifdef CONTINUUM
    // dealing with continuum
    bundle_params.ContinuumPol = specex::Legendre1DPol(param_degw["CONT"],WAVEMIN,WAVEMAX);
    for(int i=0;i<param_degw["CONT"]+1;i++) {
      bundle_params.ContinuumPol.coeff[i] = param_coeff["CONT"][i];
    }
#endif    

    // loop on parameters
    for(int pi=0;pi<(int)params.size();pi++) {      
      const string& pname=params[pi];
      if(pname=="GH-0-0") continue; // not used in c++ version
      if(pname=="CONT") continue; // not a local param
      if(pname=="BUNDLE") continue; // not a local param
      if(pname=="STATUS") continue; // not a local param
      
      SPECEX_DEBUG("Fitting pol of parameter " << pname << " in bundle " << bundle << " fibers in [" << bundle_fibermin << "," << bundle_fibermax << "]" );
      
      // now we need to fit a 2D legendre polynomial of X and wave
      // for the subset of fibers of this bundle

      int degx = param_degx[pname];
      int degw = param_degw[pname];
      
      specex::Pol_p pol(new specex::Pol(degx,xmin,xmax,degw,WAVEMIN,WAVEMAX));
      pol->name = pname;
      pol->Fill(true); // sparse or not sparse ????
      
      int npar = pol->Npar();
      unhrp::matrix_double A(npar,npar); A.clear();
      unhrp::vector_double B(npar); B.clear();
      
      int npoints=0;
      for(int fiber=bundle_fibermin;fiber<=bundle_fibermax;fiber++) {
	const specex::Trace& trace = psf->FiberTraces[fiber];
	specex::Legendre1DPol fiberpol(LEGDEG,WAVEMIN,WAVEMAX);
	for(int cj=0;cj<=LEGDEG;cj++)
	  fiberpol.coeff[cj]=param_coeff[pname][cj+(fiber-FIBERMIN)*(LEGDEG+1)];
	
	for(double wave=WAVEMIN;wave<WAVEMAX+0.01;wave+=(WAVEMAX-WAVEMIN)/(degw+1)) {
	  double x    = trace.X_vs_W.Value(wave);
	  double pval = fiberpol.Value(wave);
	  unhrp::vector_double der = pol->Monomials(x,wave);
	  specex::syr(1.,der,A); // A += der*der.transposed;
	  specex::axpy(pval,der,B); // B += pval*der;
	  npoints += 1;
	}
      }
      // now need to solve
      int status = specex::cholesky_solve(A,B);
      if(status != 0) 
	SPECEX_ERROR("Oups, failed to convert LegPol(fiber,wave) -> LegPol(x,wave) for bundle " << bundle << " and parameter " << pname << " " << "bundle fibermin,fibermax=" << bundle_fibermin << "," << bundle_fibermax << " degx,degw=" << degx << "," << degw << " " << "xmin,xmax=" << xmin << "," << xmax << " WAVEMIN,WAVEMAX=" << WAVEMIN << "," << WAVEMAX << " npoints=" << npoints );
      
      pol->coeff = B;
      bundle_params.AllParPolXW.push_back(pol);
       
    }    
    psf->ParamsOfBundles[bundle] = bundle_params;
  } 
}


