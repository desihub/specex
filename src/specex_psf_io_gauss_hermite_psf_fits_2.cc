#include <fstream>
#include <boost/algorithm/string.hpp>
#include <unhrp.h>
#include <specex_fits.h>
#include <specex_trace.h>
#include <specex_gauss_hermite_psf.h>
#include <specex_string.h>
#include <specex_unbst.h>

using namespace std ;

static void AddRow1(specex::FitsTable& table,const string& PARAM, double wavemin, double wavemax, unhrp::vector_double& coeff) {
  std::vector<specex::FitsTableEntry> row;
  {specex::FitsTableEntry entry; entry.string_val = PARAM; row.push_back(entry);}
  {specex::FitsTableEntry entry; entry.double_vals.resize(1); entry.double_vals[0] = wavemin; row.push_back(entry);}
  {specex::FitsTableEntry entry; entry.double_vals.resize(1); entry.double_vals[0] = wavemax; row.push_back(entry);}
  {specex::FitsTableEntry entry; entry.double_vals = coeff; row.push_back(entry);}
  table.data.push_back(row);
}

void read_gauss_hermite_psf_fits_version_2(specex::PSF_p& psf, fitsfile* fp, int hdu, int first_bundle, int last_bundle) {
  
  SPECEX_DEBUG("read_gauss_hermite_psf_fits_version_2 hdu=" << hdu << " first,last bundle=" << first_bundle << "," << last_bundle);

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
  
  harp::fits::key_read(fp,"MJD",psf->mjd);
  harp::fits::key_read(fp,"PLATEID",psf->plate_id);
  harp::fits::key_read(fp,"CAMERA",psf->camera_id);
  harp::fits::key_read(fp,"ARCEXP",psf->arc_exposure_id);
  
  int nx; harp::fits::key_read(fp,"NPIX_X",nx);
  psf->ccd_image_n_cols = nx;
  int ny; harp::fits::key_read(fp,"NPIX_Y",ny);
  psf->ccd_image_n_rows = ny;
  harp::fits::key_read(fp,"HSIZEX",psf->hSizeX);
  harp::fits::key_read(fp,"HSIZEY",psf->hSizeY);
  
  
  int bundlemin; harp::fits::key_read(fp,"BUNDLMIN",bundlemin);
  int bundlemax; harp::fits::key_read(fp,"BUNDLMAX",bundlemax);
  int fibermin; harp::fits::key_read(fp,"FIBERMIN",fibermin);
  int fibermax; harp::fits::key_read(fp,"FIBERMAX",fibermax);
  int nparams; harp::fits::key_read(fp,"NPARAMS",nparams);
  int legdeg; harp::fits::key_read(fp,"LEGDEG",legdeg);
  int legdeg_param_fit=0;
  int legdeg_xtrace_fit=0;
  int legdeg_ytrace_fit=0;
  try { 
    harp::fits::key_read(fp,"PALEGDEG",legdeg_param_fit);
    harp::fits::key_read(fp,"TXLEGDEG",legdeg_xtrace_fit);
    harp::fits::key_read(fp,"TYLEGDEG",legdeg_ytrace_fit);
  } catch(...) {
    SPECEX_DEBUG("No key PALEGDEG,TXLEGDEG or TYLEGDEG");
  }
  if(legdeg_param_fit>0) {
    SPECEX_DEBUG("Will use fitted legendre degrees"); 
  }else{
    legdeg_param_fit=legdeg-1;
    legdeg_xtrace_fit=legdeg-1;
    legdeg_ytrace_fit=legdeg-1;
  }
    
    
  SPECEX_DEBUG("legdeg_param_fit  = " << legdeg_param_fit);
  SPECEX_DEBUG("legdeg_xtrace_fit = " << legdeg_xtrace_fit);
  SPECEX_DEBUG("legdeg_ytrace_fit = " << legdeg_ytrace_fit);
  
  specex::FitsTable table;
  table.Read(fp);
  // 
  int param_col = table.columns["PARAM"].col;
  int wmin_col  = table.columns["WAVEMIN"].col;
  int wmax_col  = table.columns["WAVEMAX"].col;
  int coeff_col = table.columns["COEFF"].col;
  std::vector<std::string> params;
  std::map<std::string,int> param_row;
  std::map<std::string,double> param_wavemin;
  std::map<std::string,double> param_wavemax;
  std::map<std::string,unhrp::vector_double > param_coeff;
  for(int i=0;i<table.data.size();i++) { 
    std::string pname=table.data[i][param_col].string_val;
    trim(pname);
    SPECEX_DEBUG("read_gauss_hermite_psf " << i << " '" << pname << "'");
    params.push_back(pname);
    param_row[pname]=i;
    param_wavemin[pname]=table.data[i][wmin_col].double_vals[0];
    param_wavemax[pname]=table.data[i][wmax_col].double_vals[0];
    param_coeff[pname]=table.data[i][coeff_col].double_vals;    
  }

  
  if(first_bundle>=0 && last_bundle>=first_bundle && (bundlemin<first_bundle || bundlemax>last_bundle)) {
    SPECEX_DEBUG("Truncating input PSF bundles [" << bundlemin << "," << bundlemax << "] -> [" 
		 << first_bundle << "," << last_bundle << "]");
    int original_nfibers  = (fibermax-fibermin)+1;
    int original_nbundles = (bundlemax-bundlemin)+1;
    int nfibers_per_bundle = original_nfibers/original_nbundles; // not safe ...
    int new_fiber_begin = first_bundle*nfibers_per_bundle;
    int new_fiber_end   = (last_bundle+1)*nfibers_per_bundle;
    for(size_t s=0;s<params.size();s++) {
      string pname=params[s];
      //SPECEX_DEBUG("Truncated parameter " << pname);
      int ncoef_per_fiber = param_coeff[pname].size()/original_nfibers;
      param_coeff[pname] = specex::unbst::subrange(param_coeff[pname],new_fiber_begin*ncoef_per_fiber,new_fiber_end*ncoef_per_fiber);
    }
    fibermin=new_fiber_begin;
    fibermax=new_fiber_end-1;
    bundlemin=first_bundle;
    bundlemax=last_bundle;
    SPECEX_DEBUG("Truncated FIBERMIN " << fibermin);
    SPECEX_DEBUG("Truncated FIBERMAX " << fibermax);
    SPECEX_DEBUG("Truncated BUNDLEMIN " << bundlemin);
    SPECEX_DEBUG("Truncated BUNDLEMAX " << bundlemax);    
  }

  SPECEX_DEBUG("Reading and converting FiberTraces");
  // FiberTraces
  // check size makes sense
  int nfibers = fibermax-fibermin+1;
  int ncoeff  = nfibers*(legdeg+1);
  if( int(param_coeff["X"].size()) != ncoeff)
    SPECEX_ERROR("XCOEFF SIZE INCONSISTENT WITH LEGDEG AND NUMBER OF FIBERS");
  int index=0;
  for(int fiber=fibermin;fiber<=fibermax; fiber++,index++) {

    SPECEX_DEBUG("Fiber #" << fiber);

    specex::Trace trace;
    trace.fiber=fiber;
    trace.mask=0; //?
    trace.X_vs_W = specex::Legendre1DPol(legdeg_xtrace_fit,param_wavemin["X"],param_wavemax["X"]);
    trace.Y_vs_W = specex::Legendre1DPol(legdeg_ytrace_fit,param_wavemin["Y"],param_wavemax["Y"]);
    for(int i=0;i<legdeg_xtrace_fit+1;i++) {
      trace.X_vs_W.coeff[i] = param_coeff["X"][i+index*(legdeg+1)];
    }
    for(int i=0;i<legdeg_ytrace_fit+1;i++) {
      trace.Y_vs_W.coeff[i] = param_coeff["Y"][i+index*(legdeg+1)];
    }
    trace.X_vs_Y.deg  = trace.Y_vs_W.deg + 1; // add one for inversion
    trace.W_vs_Y.deg  = trace.Y_vs_W.deg + 1; // add one for inversion
    int npts=100;
    unhrp::vector_double w(npts);
    unhrp::vector_double x(npts);
    unhrp::vector_double y(npts);
    for(int i=0;i<npts;i++) {
      w[i]=param_wavemin["Y"]+((param_wavemax["Y"]-param_wavemin["Y"])/(npts-1))*i;
      x[i]=trace.X_vs_W.Value(w[i]);
      y[i]=trace.Y_vs_W.Value(w[i]);      
    }
    trace.X_vs_Y.Fit(y,x);
    trace.W_vs_Y.Fit(y,w);
    trace.synchronized = true;
    psf->FiberTraces[fiber] = trace;
  }

  // Now that's more complicated, need to loop on bundles
  int nbundles = (bundlemax-bundlemin)+1;
  int nfibers_per_bundle = nfibers/nbundles; // not safe ...
  SPECEX_DEBUG("Number of bundles= " << nbundles);
  SPECEX_DEBUG("Number of fibers per bundle= " << nfibers_per_bundle);
  
  for(int bundle = bundlemin ; bundle <= bundlemax ; bundle++) {
    
    
    //SPECEX_DEBUG("Fitting pols of bundle " << bundle);

    int bundle_fibermin=bundle*nfibers_per_bundle;
    int bundle_fibermax=(bundle+1)*nfibers_per_bundle-1;
    if(bundle_fibermin<fibermin || bundle_fibermax>fibermax) {
      SPECEX_ERROR("FIBERMIN,FIBERMAX,BUNDLEMIN,BUNDLEMAX ARE INCONSISTENT");
    }
    specex::PSF_Params bundle_params;
    bundle_params.bundle_id=bundle;
    bundle_params.fiber_min=bundle_fibermin;
    bundle_params.fiber_max=bundle_fibermax;
    
    double xmin=100000;
    double xmax=0;
    double wavemin = param_wavemin["X"]; // or any other param
    double wavemax = param_wavemax["X"]; // or any other param
    
    for(int fiber=bundle_fibermin;fiber<=bundle_fibermax;fiber++) {
      const specex::Trace& trace = psf->FiberTraces[fiber];
      for(double wave=wavemin;wave<=wavemax;wave+=10) {
	double x = trace.X_vs_W.Value(wave);
	xmin=min(xmin,x);
	xmax=max(xmax,x);	  
      }
    }
    SPECEX_DEBUG("bundle=" << bundle << " xmin xmax : " << xmin << " " << xmax);
#ifdef CONTINUUM
    // dealing with continuum
    bundle_params.ContinuumPol = specex::Legendre1DPol(legdeg,param_wavemin["CONT"],param_wavemax["CONT"]);
    for(int i=0;i<legdeg+1;i++) {
      bundle_params.ContinuumPol.coeff[i] = param_coeff["CONT"][i+(bundle_params.fiber_min-fibermin)*(legdeg+1)];
    }
#endif    
    

    // need to set : bundle_params.AllParPolXW  and bundle_params.FitParPolXW
    for(int i=0;i<(int)params.size();i++) {      
      const string& pname=params[i];
      
      SPECEX_DEBUG("Fitting pol of parameter " << pname << " in bundle " << bundle);

      if(pname=="X") continue; // trace
      if(pname=="Y") continue; // trace
      if(pname=="GH-0-0") continue; // not used in c++ version
      if(pname=="CONT") continue; // not a local param
      // now we need to fit a 2D legendre polynomial of X and wave
      // for the subset of fibers of this bundle

      // NOTE : we remove 1 deg along wave, because 1 was added in the c++/xml -> fits conversion X->fiber
      specex::Pol_p pol(new specex::Pol(nfibers_per_bundle-1,xmin,xmax,legdeg_param_fit,wavemin,wavemax));
      pol->name = pname;
      pol->Fill(true); // sparse or not sparse

      int npar = pol->Npar();
      unhrp::matrix_double A(npar,npar);
      unhrp::vector_double B(npar);
      A *= 0;
      //B *= 0;
      for(int i = 0; i < B.size(); i++) B[i] = 0;
      for(int fiber=bundle_fibermin;fiber<=bundle_fibermax;fiber++) {
	const specex::Trace& trace = psf->FiberTraces[fiber];
	specex::Legendre1DPol fiberpol(legdeg,param_wavemin[pname],param_wavemax[pname]);
	for(int i=0;i<=legdeg;i++)
	  fiberpol.coeff[i]=param_coeff[pname][i+(fiber-fibermin)*(legdeg+1)];
	for(double wave=wavemin;wave<=wavemax;wave+=(wavemax-wavemin)/(legdeg_param_fit+2)) {
	  double x    = trace.X_vs_W.Value(wave);
	  double pval = fiberpol.Value(wave);
	  unhrp::vector_double der = pol->Monomials(x,wave);
	  specex::syr(1.,der,A); // A += der*der.transposed;
	  specex::axpy(pval,der,B); // B += pval*der;
	}
      }
      // now need to solve
      int status = specex::cholesky_solve(A,B);
      if(status != 0) 
	SPECEX_ERROR("Failed to convert LegPol(fiber,wave) -> LegPol(x,wave) for bundle " << bundle << " and parameter " << pname);
      
      pol->coeff = B;
      //SPECEX_DEBUG("Fit of parameter " << pname << " in bundle " << bundle << " done");
      
      bundle_params.AllParPolXW.push_back(pol);
      
      
    }
    
    
    psf->ParamsOfBundles[bundle] = bundle_params;
  } 

  //SPECEX_ERROR("Finish implementation");
}


