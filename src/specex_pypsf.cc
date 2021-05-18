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
    unhrp::vector_double ty(deg+ddeg+1);
    unhrp::vector_double tx(deg+ddeg+1);
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

void specex::PyPSF::get_table(std::vector<std::string>         &table_string,
			      std::vector<std::vector<int>>    &table_int,
			      std::vector<std::vector<double>> &table_double){

  int status = 0;
  int nrows  = this->psf->pydata.table.data.size();
  int ncols  = this->psf->pydata.table.columns.size();

  cout << "nrows in get_tablestring = " << nrows << endl;
  cout << "ncols in get_tablestring = " << ncols << endl;
  
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
