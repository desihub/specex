#include <specex_psf.h>
#include <specex_psf_proc.h>
#include <specex_spot.h>
#include <specex_message.h>
#include <specex_image_data.h>
#include <specex_fits.h>
#include <specex_trace.h>

#include <specex_unhrp.h>

namespace unhrp = specex::unhrp;

static void _AddRow2(specex::FitsTable& table,const string& PARAM, unhrp::vector_double& coeff, int legdegx, int legdegw) {
  std::vector<specex::FitsTableEntry> row;
  {specex::FitsTableEntry entry; entry.string_val = PARAM; row.push_back(entry);}
  {specex::FitsTableEntry entry; entry.double_vals = coeff; row.push_back(entry);}
  {specex::FitsTableEntry entry; entry.int_vals.resize(1); entry.int_vals[0] = legdegx; row.push_back(entry);}
  {specex::FitsTableEntry entry; entry.int_vals.resize(1); entry.int_vals[0] = legdegw; row.push_back(entry);}
  table.data.push_back(row);
}


void _load_trace(specex::PSF_p psf, bool is_x) {
  
  int ncoeff=0;
  double WAVEMIN=0;
  double WAVEMAX=0;
  for(std::map<int,specex::Trace>::const_iterator it=psf->FiberTraces.begin(); it!=psf->FiberTraces.end(); it++) {
    const specex::Legendre1DPol *pol=0;
    if(is_x) pol = &(it->second.X_vs_W);
    else pol = &(it->second.Y_vs_W);
    if((pol->xmin)==0) continue;
    ncoeff=max(ncoeff,int(pol->coeff.size()));
    if(WAVEMIN==0) WAVEMIN=pol->xmin;
    else if(WAVEMIN != pol->xmin) {SPECEX_ERROR("requires same WAVEMIN for all traces and have" << WAVEMIN << " != " << pol->xmin);}
    if(WAVEMAX==0) WAVEMAX=pol->xmax;
    else if(WAVEMAX != pol->xmax) {SPECEX_ERROR("requires same WAVEMAX for all traces and have" << WAVEMAX << " != " << pol->xmax);} 
  }
  WAVEMIN=floor(WAVEMIN);
  WAVEMAX=floor(WAVEMAX);
  
  int FIBERMIN=100000;
  int FIBERMAX=0;
  for(std::map<int,specex::Trace>::const_iterator it=psf->FiberTraces.begin(); it!=psf->FiberTraces.end(); it++) {
    FIBERMIN=min(FIBERMIN,it->second.fiber);
    FIBERMAX=max(FIBERMAX,it->second.fiber);
    
  }
  int NFIBERS=(FIBERMAX-FIBERMIN+1);

  specex::image_data coeff2d(ncoeff,NFIBERS);
  
  for(std::map<int,specex::Trace>::const_iterator it=psf->FiberTraces.begin(); it!=psf->FiberTraces.end(); it++) {
    int fiber = it->second.fiber;
    const specex::Legendre1DPol *pol=0;
    if(is_x) pol = &(it->second.X_vs_W);
    else pol = &(it->second.Y_vs_W);

    if( (pol->xmin==WAVEMIN) && (pol->xmax==WAVEMAX) ) {
      for(int c=0;c<pol->coeff.size();c++)
	coeff2d(c,fiber-FIBERMIN)=pol->coeff(c);
    }else{ // need to refit
      if(pol->coeff.size()>0) {
	SPECEX_DEBUG("need to refit trace coeff.size=" << pol->coeff.size());
	try {
	  unhrp::vector_double wave(pol->coeff.size());
	  unhrp::vector_double x(pol->coeff.size());
	  for(int i=0;i<pol->coeff.size();i++) {
	    wave[i]=WAVEMIN+i*((WAVEMAX-WAVEMIN)/(pol->coeff.size()-1));
	    x[i]=pol->Value(wave[i]);
	  }
	  specex::Legendre1DPol npol(pol->coeff.size()-1,WAVEMIN,WAVEMAX);
	  npol.Fit(wave,x,0,false);
	  for(int c=0;c<npol.coeff.size();c++)
	    coeff2d(c,fiber-FIBERMIN)=npol.coeff(c);
	}catch(std::exception& e) {
	  SPECEX_ERROR("Fit failed " << e.what());
	}
      }
    }
  }

  SPECEX_DEBUG("done write trace in HDU ");

  // save key data 
  psf->pydata.FIBERMIN = FIBERMIN;
  psf->pydata.FIBERMAX = FIBERMAX;
  psf->pydata.trace_WAVEMIN = WAVEMIN;
  psf->pydata.trace_WAVEMAX = WAVEMAX;
    
  // save value of coeff2d 
  psf->pydata.SetCoeff2d(coeff2d,is_x);
  
}
  
void _load_xtrace(specex::PSF_p psf) {
  SPECEX_INFO("prepare XTRACE");
  return _load_trace(psf,true); // X 
}
void _load_ytrace(specex::PSF_p psf) {
  SPECEX_INFO("write YTRACE");
  return _load_trace(psf,false); // Y 
}

void _load_psf(specex::PSF_p psf) {
  
  SPECEX_INFO("load PSF");
  SPECEX_DEBUG("load_gauss_hermite_psf_version_3");
  
  ////////////////////////////
  string PSFVER = "3";
  ////////////////////////////
  
  int NPIX_X = psf->ccd_image_n_cols;
  int NPIX_Y = psf->ccd_image_n_rows;
  int GHDEGX = psf->Degree();
  int GHDEGY = psf->Degree();

  // FIBERMIN FIBERMAX WAVEMIN WAVEMAX are always defined by the trace set
  // to ensure consistency
  double WAVEMIN=1000000;
  double WAVEMAX=0;
  int    FIBERMIN=100000;
  int    FIBERMAX=0;
  for(std::map<int,specex::Trace>::const_iterator it=psf->FiberTraces.begin(); it!=psf->FiberTraces.end(); it++) {
    FIBERMIN=min(FIBERMIN,it->second.fiber);
    FIBERMAX=max(FIBERMAX,it->second.fiber);
    if(it->second.Y_vs_W.xmin>0)
      WAVEMIN=min(WAVEMIN,it->second.Y_vs_W.xmin);
    WAVEMAX=max(WAVEMAX,it->second.Y_vs_W.xmax);    
  }
  WAVEMIN = floor(WAVEMIN);
  WAVEMAX = floor(WAVEMAX)+1.;
  int NFIBERS=(FIBERMAX-FIBERMIN+1);
  
  int status = 0;
    
  // get largest legendre degree of all params in all bundles and check param numbers match !!
  
  int nparams=psf->LocalNAllPar();
  int nparams_all = nparams;
  nparams_all += 1; // GH trivial order zero
#ifdef CONTINUUM
  nparams_all += 1; // continuum
#endif
  
  int ncoeff=0;
  for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf->ParamsOfBundles.begin();
      bundle_it != psf->ParamsOfBundles.end(); ++bundle_it) {
    if( int(bundle_it->second.AllParPolXW.size()) != nparams ) SPECEX_ERROR("Fatal inconsistency in number of parameters in psf between bundles: AllParPolXW.size=" << bundle_it->second.AllParPolXW.size() << " psf->LocalNAllPar()=" << nparams);
    for(size_t p=0;p<bundle_it->second.AllParPolXW.size();p++) {
      ncoeff=max(ncoeff,bundle_it->second.AllParPolXW[p]->ydeg+1);
    }
  }
  SPECEX_DEBUG("ncoeff = " << ncoeff);
  psf->ncoeff = ncoeff;
  
  specex::FitsTable table;
  
  { // column description
    table.AddColumnDescription("PARAM","8A","","");
    char coeff_tform[100];
    sprintf(coeff_tform,"%dD",ncoeff*NFIBERS);
    vector <int> dim; dim.resize(2);
    dim[0]=ncoeff;
    dim[1]=NFIBERS;
    
    string sdim = table.encode_dimension(dim);
    table.AddColumnDescription("COEFF",coeff_tform,sdim,"");
    dim.resize(1); dim[0]=1;
    sdim = table.encode_dimension(dim);
    table.AddColumnDescription("LEGDEGX","J",sdim,"");
    table.AddColumnDescription("LEGDEGW","J",sdim,"");
  }

  unhrp::vector_double wave(ncoeff);
  {
    double wavestep = (WAVEMAX-WAVEMIN)/(ncoeff-1);
    for(int w=0;w<ncoeff;w++) {
      wave[w]   = WAVEMIN + wavestep*w;
    }
  }
  
  vector<string> keys;
  
  unhrp::vector_double coeff(ncoeff*NFIBERS);
  unhrp::vector_double values(ncoeff);
  
  bool need_to_add_first_gh = true;
  for(int p=0;p<nparams;p++) {  // loop on all PSF parameters      
        
    string pname = psf->ParamName(p);    
    if(need_to_add_first_gh && pname.find("GH-")<pname.npos) { 
      // insert now GH param 0 which is not a PSF parameter in C++ but needed but specter
      
      for(int fiber=0;fiber<NFIBERS;fiber++) 
	coeff(fiber*ncoeff)=1; 
      _AddRow2(table,"GH-0-0",coeff,0,0);
      need_to_add_first_gh = false;
    }
        
    int legdegx=0;
    int legdegw=0;
    // loop on bundles
    for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf->ParamsOfBundles.begin();
	bundle_it != psf->ParamsOfBundles.end(); ++bundle_it) {
      const specex::PSF_Params & params_of_bundle = bundle_it->second;      
      const specex::Pol_p pol2d = params_of_bundle.AllParPolXW[p]; // this is the 2D polynomiald of x_ccd and wave for this param and bundle
      legdegx=max(legdegx,pol2d->xdeg);
      legdegw=max(legdegw,pol2d->ydeg);

      for(int fiber=params_of_bundle.fiber_min; fiber<=params_of_bundle.fiber_max; fiber++) {
	
	const specex::Trace& trace = psf->FiberTraces.find(fiber)->second; // X_vs_W.Value();
	
	// build a Legendre1DPol out of the Legendre2DPol
	specex::Legendre1DPol pol1d(ncoeff-1,WAVEMIN,WAVEMAX);
	for(int w=0;w<ncoeff;w++) {
	  values[w] = pol2d->Value(trace.X_vs_W.Value(wave[w]),wave[w]);
	}
	pol1d.Fit(wave,values,0,false);
	
	// now copy parameters;	
	for(int w = 0; w < ncoeff ; w++) {
	  coeff((fiber-FIBERMIN)*ncoeff+w) = pol1d.coeff(w); // this is the definition of the ordering, (wave,fiber)
	}	
      } // end of loop on fibers of bundle      
    } // end of loop on bundles

    _AddRow2(table,pname,coeff,legdegx,legdegw); 
  } // end of loop on params

  { // add a parameter to link fibers and bundles in fit
    coeff.clear();
    for(int fiber=FIBERMIN;fiber<=FIBERMAX;fiber++)
      coeff((fiber-FIBERMIN)*ncoeff) = -1; // no bundle
    for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf->ParamsOfBundles.begin();
	bundle_it != psf->ParamsOfBundles.end(); ++bundle_it) {
      const specex::PSF_Params & params_of_bundle = bundle_it->second;
      for(int fiber=params_of_bundle.fiber_min; fiber<=params_of_bundle.fiber_max; fiber++) {
	coeff((fiber-FIBERMIN)*ncoeff) = bundle_it->first;
      }
    }
    _AddRow2(table,"BUNDLE",coeff,0,0); 
  }
  { // add a parameter with the fit status
    coeff.clear();
    for(int fiber=FIBERMIN;fiber<=FIBERMAX;fiber++)
      coeff((fiber-FIBERMIN)*ncoeff) = -1; // no bundle
    for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf->ParamsOfBundles.begin();
	bundle_it != psf->ParamsOfBundles.end(); ++bundle_it) {
      const specex::PSF_Params & params_of_bundle = bundle_it->second;
      for(int fiber=params_of_bundle.fiber_min; fiber<=params_of_bundle.fiber_max; fiber++) {
	coeff((fiber-FIBERMIN)*ncoeff) = params_of_bundle.fit_status;
      }
    }
    _AddRow2(table,"STATUS",coeff,0,0); 
  }
  
#ifdef CONTINUUM  
  { // add a parameter for continuum
    coeff.clear();
    int fiber_index=0;
    int legdegx=0;
    int legdegw=0;
    for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf->ParamsOfBundles.begin();
	bundle_it != psf->ParamsOfBundles.end(); ++bundle_it) {
      const specex::PSF_Params & params_of_bundle = bundle_it->second;
      legdegw=params_of_bundle.ContinuumPol.deg;
      specex::Legendre1DPol pol1d(ncoeff-1,WAVEMIN,WAVEMAX);
      for(int w=0;w<ncoeff;w++) {
	values[w]   = params_of_bundle.ContinuumPol.Value(wave[w]);
      }
      pol1d.Fit(wave,values,0,false);
      for(int fiber=params_of_bundle.fiber_min; fiber<=params_of_bundle.fiber_max; fiber++,fiber_index++) {
	for(int w = 0; w < ncoeff ; w++) {
	  coeff((fiber-FIBERMIN)*ncoeff+w)   =  pol1d.coeff(w);
	}    
      }
    }
    _AddRow2(table,"CONT",coeff,legdegx,legdegw);
  }
#endif    
    
  // save key data
  psf->pydata.table_WAVEMIN = WAVEMIN;
  psf->pydata.table_WAVEMAX = WAVEMAX;
  
  // save table in pydata object for export to python
  psf->pydata.table = table; 

} // end of routine

void specex::load_psf_work(specex::PSF_p psf){
    
  _load_xtrace(psf);
  _load_ytrace(psf);
  _load_psf(psf);
  SPECEX_INFO("loaded [x,y]trace and PSF");
  
}
