#include <fstream>
#include <string>
#include <iostream>
#include <specex_pypsf.h>

using namespace std;

specex::image_data specex::PyPSF::get_trace(std::string axis){
  
  if(axis == "x"){
    return this->psf->pydata.coeff2d_x;
  }
  else{
    return this->psf->pydata.coeff2d_y;
  }
}

void specex::PyPSF::set_psf(
			    std::vector<std::string> &table_col0,
			    std::vector<double>      &table_col1,
			    std::vector<int>         &table_col2,
			    std::vector<int>         &table_col3
			    ){
 
  int status = 0;

  int GHDEGX = this->GHDEGX;
  int GHDEGY = this->GHDEGY;
  if(GHDEGX != GHDEGY) {
    SPECEX_ERROR("expect GHDEGX=GHDEGY");
  }
  
  ((specex::GaussHermitePSF*)(&(*psf)))->SetDegree(GHDEGX);

  this->psf->mjd              = this->mjd;
  this->psf->plate_id         = this->plate_id;
  this->psf->camera_id        = this->camera_id;
  this->psf->arc_exposure_id  = this->arc_exposure_id;
  this->psf->ccd_image_n_cols = this->NPIX_X;
  this->psf->ccd_image_n_rows = this->NPIX_Y;
  this->psf->hSizeX           = this->hSizeX;
  this->psf->hSizeY           = this->hSizeY;
  int FIBERMIN                = this->FIBERMIN;
  int FIBERMAX                = this->FIBERMAX;
  double WAVEMIN              = this->table_WAVEMIN;
  double WAVEMAX              = this->table_WAVEMAX;
  int LEGDEG                  = this->LEGDEG;

  std::vector<std::string> params;
  std::map<std::string,int> param_row;
  std::map<std::string,unbls::vector_double > param_coeff;
  std::map<std::string,int > param_degx;
  std::map<std::string,int > param_degw;
    
  int ncoeff_per_row = this->nfibers*this->ncoeff;
  
  for(int i=0; i < this->table_nrows;i++) { 
    std::string pname=table_col0[i];
    params.push_back(pname);
    param_row[pname]=i;
    for(int j=0; j < ncoeff_per_row; j++){
      param_coeff[pname] = unbls::vector_double(ncoeff_per_row);
      param_coeff[pname][j] = table_col1[i*ncoeff_per_row+j];
    }
    param_degx[pname]=table_col2[i];
    param_degw[pname]=table_col3[i];
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
      unbls::matrix_double A(npar,npar);
      unbls::zero(A);
      unbls::vector_double B(npar,0.);
      
      int npoints=0;
      for(int fiber=bundle_fibermin;fiber<=bundle_fibermax;fiber++) {
	const specex::Trace& trace = psf->FiberTraces[fiber];
	specex::Legendre1DPol fiberpol(LEGDEG,WAVEMIN,WAVEMAX);
	for(int cj=0;cj<=LEGDEG;cj++)
	  fiberpol.coeff[cj]=param_coeff[pname][cj+(fiber-FIBERMIN)*(LEGDEG+1)];
	
	for(double wave=WAVEMIN;wave<WAVEMAX+0.01;wave+=(WAVEMAX-WAVEMIN)/(degw+1)) {
	  double x    = trace.X_vs_W.Value(wave);
	  double pval = fiberpol.Value(wave);
	  unbls::vector_double der = pol->Monomials(x,wave);
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

void specex::PyPSF::init_traces(specex::PyOptions opts){

  SPECEX_INFO("Initializing a " << opts.psf_model << " PSF");
  this->psf = PSF_p(new specex::GaussHermitePSF(opts.gauss_hermite_deg));
    
  this->psf->hSizeX = opts.half_size_x;
  this->psf->hSizeY = opts.half_size_y;
  SPECEX_INFO("trace_deg_x=" << opts.trace_deg_x << " trace_deg_wave=" <<
	      opts.trace_deg_wave);
    
}

void specex::PyPSF::synchronize_traces(){
  SPECEX_INFO("synchronizing traces");
  specex::PSF_p psf = this->psf;
  for(std::map<int,specex::Trace>::iterator it=psf->FiberTraces.begin(); it!=psf->FiberTraces.end(); it++) {
    specex::Trace &trace = it->second;
    // X_vs_Y
    int deg     = trace.Y_vs_W.deg;
    double wmin = trace.Y_vs_W.xmin;
    double wmax = trace.Y_vs_W.xmax;
    
    int ddeg = 1; // add one degree for inversion
    unbls::vector_double ty(deg+ddeg+1);
    unbls::vector_double tx(deg+ddeg+1);
    for(int i=0;i<deg+ddeg+1;i++) {
      double wave=wmin+i*((wmax-wmin)/deg);
      ty[i]=trace.Y_vs_W.Value(wave);
      tx[i]=trace.X_vs_W.Value(wave);
    }
    trace.X_vs_Y = specex::Legendre1DPol(deg+ddeg,0,4000);
    trace.X_vs_Y.Fit(ty,tx,0,true);    
    trace.W_vs_Y = trace.Y_vs_W.Invert(ddeg);
    trace.synchronized=true;
    
    if(it==psf->FiberTraces.begin()) {
      SPECEX_INFO("X_vs_W deg=" << trace.X_vs_W.deg);
      SPECEX_INFO("Y_vs_W deg=" << trace.Y_vs_W.deg);
      SPECEX_INFO("X_vs_Y deg=" << trace.X_vs_Y.deg);
      SPECEX_INFO("W_vs_Y deg=" << trace.W_vs_Y.deg);      
    }    
  }
  
}

void specex::PyPSF::set_trace(py::array trace, int requested_deg, int isx) {

  auto trace_prop = trace.request();
  double *trace_vals = (double*) trace_prop.ptr;

  // assume both traces equal in size
  int ncoefs  = this->trace_ncoeff;
  int nfibers = this->nfibers;

  double WAVEMIN = this->trace_WAVEMIN;
  double WAVEMAX = this->trace_WAVEMAX;
  
  int requested_ncoefs = ncoefs;  
  
  if(requested_deg>0)  requested_ncoefs=requested_deg+1;
  
  for(int fiber=0;fiber<nfibers; fiber++) {
    // create or not the trace
    specex::Trace& trace = this->psf->FiberTraces[fiber];
    trace.fiber=fiber;
    trace.mask=0;
    trace.synchronized = false; // will do synchro after
    specex::Legendre1DPol pol(requested_ncoefs-1,WAVEMIN,WAVEMAX);
    for(int i=0;i<min(ncoefs,requested_ncoefs);i++)
      pol.coeff[i] = trace_vals[fiber*ncoefs + i]; //coeff(i,fiber);
    if(isx==1) trace.X_vs_W=pol;
    else trace.Y_vs_W=pol;
  }
  
  return ;
  
}

void specex::PyPSF::SetParamsOfBundle(){

  for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = this->psf->ParamsOfBundles.begin();
      bundle_it != this->psf->ParamsOfBundles.end(); ++bundle_it) {
      
    const specex::PSF_Params & params_of_bundle = bundle_it->second;
      
    int ndf = params_of_bundle.ndata - params_of_bundle.nparams;
    double chi2pdf = 0;
    
    if(ndf>0) chi2pdf = params_of_bundle.chi2/ndf;

    bundle_id.push_back(params_of_bundle.bundle_id);
    bundle_ndata.push_back(params_of_bundle.ndata);
    bundle_nparams.push_back(params_of_bundle.nparams);
    bundle_chi2pdf.push_back(chi2pdf);
    
  }
}
