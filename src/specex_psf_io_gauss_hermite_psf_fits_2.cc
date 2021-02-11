#include <fstream>
#include <boost/algorithm/string.hpp>
#include <harp.hpp>
#include <specex_fits.h>
#include <specex_trace.h>
#include <specex_gauss_hermite_psf.h>

using namespace std ;

static void AddRow1(specex::FitsTable& table,const string& PARAM, double wavemin, double wavemax, harp::vector_double& coeff) {
  std::vector<specex::FitsTableEntry> row;
  {specex::FitsTableEntry entry; entry.string_val = PARAM; row.push_back(entry);}
  {specex::FitsTableEntry entry; entry.double_vals.resize(1); entry.double_vals[0] = wavemin; row.push_back(entry);}
  {specex::FitsTableEntry entry; entry.double_vals.resize(1); entry.double_vals[0] = wavemax; row.push_back(entry);}
  {specex::FitsTableEntry entry; entry.double_vals = coeff; row.push_back(entry);}
  table.data.push_back(row);
}

void write_gauss_hermite_psf_fits_version_2(const specex::GaussHermitePSF& psf, fitsfile* fp, int first_hdu, std::vector<specex::Spot_p> *spots) {
  
  SPECEX_DEBUG("write_gauss_hermite_psf_fits_version_2");
  
  ////////////////////////////
  string PSFVER = "2";
  ////////////////////////////
  
  int NPIX_X = psf.ccd_image_n_cols;
  int NPIX_Y = psf.ccd_image_n_rows;
  int BUNDLMIN = 0;
  int BUNDLMAX = 0;
  int FIBERMIN = 0;
  int FIBERMAX = 0; 
  int NFIBERS=0;
  for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf.ParamsOfBundles.begin();
      bundle_it != psf.ParamsOfBundles.end(); ++bundle_it) {
    
    if(bundle_it == psf.ParamsOfBundles.begin()) {
      BUNDLMIN = bundle_it->second.bundle_id;
      BUNDLMAX = bundle_it->second.bundle_id;
      FIBERMIN = bundle_it->second.fiber_min;
      FIBERMAX = bundle_it->second.fiber_max;
    }
    
    BUNDLMIN = min(BUNDLMIN,bundle_it->second.bundle_id);
    BUNDLMAX = max(BUNDLMAX,bundle_it->second.bundle_id);
    FIBERMIN = min(FIBERMIN,bundle_it->second.fiber_min);
    FIBERMAX = max(FIBERMAX,bundle_it->second.fiber_max);
  
    NFIBERS += (bundle_it->second.fiber_max-bundle_it->second.fiber_min+1);
  }
  
  SPECEX_DEBUG("BUNDLMIN=" << BUNDLMIN << " BUNDLMAX=" << BUNDLMAX << " FIBERMIN=" << FIBERMIN << " FIBERMAX=" << FIBERMAX << " NFIBERS=" << NFIBERS);
  
  // number of fibers per bundle from first bundle
  if(NFIBERS != (FIBERMAX+1)) {
    SPECEX_WARNING("will fill with zeros the coefficients of fibers that have not been fit, there are " << (FIBERMAX+1)-NFIBERS << " of those");
    NFIBERS=(FIBERMAX+1);
  }
  
  int GHDEGX = psf.Degree();
  int GHDEGY = psf.Degree();
  

  int status = 0;
  fits_movabs_hdu ( fp, first_hdu, NULL, &status ); harp::fits::check ( status );
  
  
  
  // get largest legendre degree of all params in all bundles and check param numbers match !!
  int ncoeff_param_fit=0;
  int nparams=psf.LocalNAllPar();
  int nparams_all = nparams;
  nparams_all += 2; // X and Y
  nparams_all += 1; // GH trivial order zero
#ifdef CONTINUUM
  nparams_all += 1; // continuum
#endif
  
  for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf.ParamsOfBundles.begin();
      bundle_it != psf.ParamsOfBundles.end(); ++bundle_it) {

    if( int(bundle_it->second.AllParPolXW.size()) != nparams ) SPECEX_ERROR("Fatal inconsistency in number of parameters in psf between bundles: AllParPolXW.size=" << bundle_it->second.AllParPolXW.size() << " psf.LocalNAllPar()=" << nparams);

    for(size_t p=0;p<bundle_it->second.AllParPolXW.size();p++) {
      ncoeff_param_fit=max(ncoeff_param_fit,bundle_it->second.AllParPolXW[p]->ydeg+1);
    }
  }

  
  
  int ncoeff_xtrace_fit=0;
  int ncoeff_ytrace_fit=0;
  
  // add traces
  for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf.ParamsOfBundles.begin();
      bundle_it != psf.ParamsOfBundles.end(); ++bundle_it) {
    const specex::PSF_Params & params_of_bundle = bundle_it->second;
    for(int fiber=params_of_bundle.fiber_min; fiber<=params_of_bundle.fiber_max; fiber++) {
      ncoeff_xtrace_fit=max(ncoeff_xtrace_fit,int(psf.FiberTraces.find(fiber)->second.X_vs_W.coeff.size()));
      ncoeff_ytrace_fit=max(ncoeff_ytrace_fit,int(psf.FiberTraces.find(fiber)->second.Y_vs_W.coeff.size()));
      
    }
  }
  int ncoeff=ncoeff_param_fit;
  ncoeff=max(ncoeff,ncoeff_xtrace_fit);
  ncoeff=max(ncoeff,ncoeff_ytrace_fit);

  int delta_deg = 1; // increase degree of legendre polynomials to minimize mapping errors
  ncoeff += delta_deg; 

  // no need to add degree for traces because same mapping
  
  SPECEX_DEBUG("ncoeff = " << ncoeff);
  

  specex::FitsTable table;

  { // column description
    table.AddColumnDescription("PARAM","8A","","");
    table.AddColumnDescription("WAVEMIN","D","","");
    table.AddColumnDescription("WAVEMAX","D","","");
    char coeff_tform[100];
    sprintf(coeff_tform,"%dD",ncoeff*NFIBERS);
    vector <int> dim; dim.resize(2);
    dim[0]=ncoeff;
    dim[1]=NFIBERS;
    
    string sdim = table.encode_dimension(dim);
    // debug
    vector <int> dim2 = table.decode_dimension(sdim);
    if (dim2.size()!=2 || dim2[0]!= dim[0]|| dim2[1]!=dim[1]) {
      SPECEX_ERROR("error encode/decode dimension");
    }

    table.AddColumnDescription("COEFF",coeff_tform,sdim,"");
  }
  
 
  int LEGWMIN=1000000;
  int LEGWMAX=0;
  vector<string> keys;
  
  

  { // data
    
    
    // now loop on real psf parameters

    harp::vector_double wave(ncoeff);
    harp::vector_double values(ncoeff);
    
    // get the max range of wavelength and convert to int 
    
    for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf.ParamsOfBundles.begin();
	  bundle_it != psf.ParamsOfBundles.end(); ++bundle_it) {
      const specex::PSF_Params & params_of_bundle = bundle_it->second;
      for(int p=0;p<nparams;p++) { 
	LEGWMIN=min(LEGWMIN,int(floor(params_of_bundle.AllParPolXW[p]->ymin)));
	LEGWMAX=max(LEGWMAX,int(floor(params_of_bundle.AllParPolXW[p]->ymax)));
      }
    }
    
    
    double wavemin = double(LEGWMIN);
    double wavemax = double(LEGWMAX);
    {
      double wavestep = (wavemax-wavemin)/(ncoeff-1);
      for(int w=0;w<ncoeff;w++) {
	wave[w]   = wavemin + wavestep*w;
      }
    }
    
    bool need_to_add_first_gh = true;

    harp::vector_double coeff(ncoeff*NFIBERS);

    // first deal with X and Y
    {
      harp::vector_double coeff_y(ncoeff*NFIBERS);
      harp::vector_double values_y(ncoeff);
      coeff.clear();
      coeff_y.clear();
      int fiber_index=0;
      
      for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf.ParamsOfBundles.begin();
	  bundle_it != psf.ParamsOfBundles.end(); ++bundle_it) {
	const specex::PSF_Params & params_of_bundle = bundle_it->second;
	for(int fiber=params_of_bundle.fiber_min; fiber<=params_of_bundle.fiber_max; fiber++,fiber_index++) {
	  
	  const specex::Trace &trace = psf.FiberTraces.find(fiber)->second;
	  specex::Legendre1DPol pol1d_x(ncoeff-1,wavemin,wavemax);
	  specex::Legendre1DPol pol1d_y(ncoeff-1,wavemin,wavemax);
	  for(int w=0;w<ncoeff;w++) {
	    values[w]   = trace.X_vs_W.Value(wave[w]);
	    values_y[w] = trace.Y_vs_W.Value(wave[w]);
	  }
	  pol1d_x.Fit(wave,values,0,false);
	  pol1d_y.Fit(wave,values_y,0,false);
	  for(int w = 0; w < ncoeff ; w++) {
	    coeff(fiber*ncoeff+w)   =  pol1d_x.coeff(w);
	    coeff_y(fiber*ncoeff+w) =  pol1d_y.coeff(w);
	  }
	  
	  
	} // end of loop on fiber
      } // end of loop on bundles
      
      AddRow1(table,"X",LEGWMIN,LEGWMAX,coeff);
      AddRow1(table,"Y",LEGWMIN,LEGWMAX,coeff_y);
      
    } // end of X Y context

    for(int p=0;p<nparams;p++) { 
      
      coeff.clear();
      
      string pname = psf.ParamName(p);

      if(need_to_add_first_gh && pname.find("GH-")<pname.npos) { // insert now GH param 0
	
	for(int fiber=0;fiber<NFIBERS;fiber++) 
	  coeff(fiber*ncoeff)=1; 
	AddRow1(table,"GH-0-0",LEGWMIN,LEGWMAX,coeff);
	need_to_add_first_gh = false;
      }

      
      int fiber_index=0;
      for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf.ParamsOfBundles.begin();
	  bundle_it != psf.ParamsOfBundles.end(); ++bundle_it) {
	const specex::PSF_Params & params_of_bundle = bundle_it->second;
	
	const specex::Pol_p pol2d = params_of_bundle.AllParPolXW[p]; // this is the 2D polynomiald of x_ccd and wave for this param and bundle
	

	
	for(int fiber=params_of_bundle.fiber_min; fiber<=params_of_bundle.fiber_max; fiber++,fiber_index++) {
	  
	  const specex::Trace& trace = psf.FiberTraces.find(fiber)->second; // X_vs_W.Value();

	  // build a Legendre1DPol out of the Legendre2DPol
	  specex::Legendre1DPol pol1d(ncoeff-1,wavemin,wavemax);
	  for(int w=0;w<ncoeff;w++) {
	    values[w] = pol2d->Value(trace.X_vs_W.Value(wave[w]),wave[w]);
	  }
	  pol1d.Fit(wave,values,0,false);

	  // now copy parameters;
	  
	  for(int w = 0; w < ncoeff ; w++) {
	    coeff(fiber*ncoeff+w) = pol1d.coeff(w); // this is the definition of the ordering, (wave,fiber)
	  }
	  
	} // end of loop on fibers of bundle
	      
      } // end of loop on bundles

      AddRow1(table,pname,LEGWMIN,LEGWMAX,coeff);
      
    } // end of loop on params

#ifdef CONTINUUM
    {
      coeff.clear();
      int fiber_index=0;
      for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf.ParamsOfBundles.begin();
	  bundle_it != psf.ParamsOfBundles.end(); ++bundle_it) {
	const specex::PSF_Params & params_of_bundle = bundle_it->second;
	specex::Legendre1DPol pol1d(ncoeff-1,wavemin,wavemax);
	for(int w=0;w<ncoeff;w++) {
	  values[w]   = params_of_bundle.ContinuumPol.Value(wave[w]);
	}
	pol1d.Fit(wave,values,0,false);
	for(int fiber=params_of_bundle.fiber_min; fiber<=params_of_bundle.fiber_max; fiber++,fiber_index++) {
	  for(int w = 0; w < ncoeff ; w++) {
	    coeff(fiber*ncoeff+w)   =  pol1d.coeff(w);
	  }    
	}
      }
      AddRow1(table,"CONT",LEGWMIN,LEGWMAX,coeff);
    }
#endif    


  }
  
  // write table
  table.Write(fp);
  
 
  // write keywords
  {
    
    fits_write_comment(fp,"------------------------------------------------------------------------",&status); harp::fits::check ( status );
    fits_write_comment(fp,"PSF generated by specex, https://github.com/julienguy/specex",&status); harp::fits::check ( status );
    
    {
      char date_comment[80];
      time_t t = time(0);   // get time now
      struct tm * now = localtime( & t );
      sprintf(date_comment,"PSF fit date %04d-%02d-%02d",(now->tm_year + 1900),(now->tm_mon + 1),now->tm_mday);
      fits_write_comment(fp,date_comment,&status); harp::fits::check ( status );
    }
    ///////////////////////xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx////////
    fits_write_comment(fp,"-  ",&status); harp::fits::check ( status );
    fits_write_comment(fp,"Each row of the table contains the data vector of one PSF parameter",&status); harp::fits::check ( status );
    fits_write_comment(fp,"The size of the vector is ((FIBERMAX-FIBERMIN+1)*(LEGDEG+1))",&status); harp::fits::check ( status );
    fits_write_comment(fp,"Description of  the NPARAMS parameters : ",&status); harp::fits::check ( status );
    fits_write_comment(fp,"X        : CCD column coordinate (as a function of fiber and wavelength)",&status); harp::fits::check ( status );
    fits_write_comment(fp,"Y        : CCD row coordinate (as a function of fiber and wavelength)",&status); harp::fits::check ( status ); 
    fits_write_comment(fp,"         (X,Y)=(0,0) means that PSF is centered on center of first pixel",&status); harp::fits::check ( status );
    fits_write_comment(fp,"GHSIGX   : Sigma of first Gaussian along CCD columns for PSF core",&status); harp::fits::check ( status );
    fits_write_comment(fp,"GHSIGY   : Sigma of first Gaussian along CCD rows for PSF core",&status); harp::fits::check ( status );
    fits_write_comment(fp,"GH-i-j   : Hermite pol. coefficents, i along columns, j along rows,",&status); harp::fits::check ( status );
    fits_write_comment(fp,"         i is integer from 0 to GHDEGX, j is integer from 0 to GHDEGY,",&status); harp::fits::check ( status );
    fits_write_comment(fp,"         there are (GHDEGX+1)*(GHDEGY+1) such coefficents.",&status); harp::fits::check ( status );
    fits_write_comment(fp,"TAILAMP  : Amplitude of PSF tail",&status); harp::fits::check ( status );
    fits_write_comment(fp,"TAILCORE : Size in pixels of PSF tail saturation in PSF core",&status); harp::fits::check ( status );
    fits_write_comment(fp,"TAILXSCA : Scaling apply to CCD coordinate along columns for PSF tail",&status); harp::fits::check ( status );
    fits_write_comment(fp,"TAILYSCA : Scaling apply to CCD coordinate along rows for PSF tail",&status); harp::fits::check ( status );
    fits_write_comment(fp,"TAILINDE : Asymptotic power law index of PSF tail",&status); harp::fits::check ( status );
    fits_write_comment(fp,"CONT     : Continuum flux in arc image (not part of PSF)",&status); harp::fits::check ( status );
    fits_write_comment(fp,"-  ",&status); harp::fits::check ( status );
    fits_write_comment(fp,"PSF_core(X,Y) = [ SUM_ij (GH-i-j)*HERM(i,X/GHSIGX)*HERM(j,Y/GHSIGX) ]",&status); harp::fits::check ( status );
    fits_write_comment(fp,"                                       *GAUS(X,GHSIGX)*GAUS(Y,GHSIGY)",&status); harp::fits::check ( status );
    fits_write_comment(fp,"-  ",&status); harp::fits::check ( status );
    fits_write_comment(fp,"PSF_tail(X,Y) = TAILAMP*R^2/(TAILCORE^2+R^2)^(1+TAILINDE/2)",&status); harp::fits::check ( status );
    fits_write_comment(fp,"                with R^2=(X/TAILXSCA)^2+(Y/TAILYSCA)^2",&status); harp::fits::check ( status );
    fits_write_comment(fp,"-  ",&status); harp::fits::check ( status );
    fits_write_comment(fp,"PSF_core is integrated in pixel",&status); harp::fits::check ( status );
    fits_write_comment(fp,"PSF_tail is not, it is evaluated at center of pixel",&status); harp::fits::check ( status );
    fits_write_comment(fp,"------------------------------------------------------------------------",&status); harp::fits::check ( status );   
    harp::fits::key_write(fp,"EXTNAME","PSF","");
    harp::fits::key_write(fp,"PSFTYPE","GAUSS-HERMITE","");
    harp::fits::key_write(fp,"PSFVER",PSFVER,"");
    
    harp::fits::key_write(fp,"MJD",(long long int)psf.mjd,"MJD of arc lamp exposure");
    harp::fits::key_write(fp,"PLATEID",psf.plate_id,"plate ID of arc lamp exposure");
    harp::fits::key_write(fp,"CAMERA",psf.camera_id,"camera ID");
    harp::fits::key_write(fp,"ARCEXP",psf.arc_exposure_id,"ID of arc lamp exposure used to fit PSF");
    
    harp::fits::key_write(fp,"NPIX_X",(long long int)NPIX_X,"number of columns in input CCD image");
    harp::fits::key_write(fp,"NPIX_Y",(long long int)NPIX_Y,"number of rows in input CCD image");
    harp::fits::key_write(fp,"HSIZEX",(long long int)psf.hSizeX,"Half size of PSF in fit, NX=2*HSIZEX+1");
    harp::fits::key_write(fp,"HSIZEY",(long long int)psf.hSizeY,"Half size of PSF in fit, NY=2*HSIZEY+1");
    harp::fits::key_write(fp,"BUNDLMIN",(long long int)BUNDLMIN,"first bundle of fibers (starting at 0)");
    harp::fits::key_write(fp,"BUNDLMAX",(long long int)BUNDLMAX,"last bundle of fibers (included)");
    harp::fits::key_write(fp,"FIBERMIN",(long long int)FIBERMIN,"first fiber (starting at 0)");
    harp::fits::key_write(fp,"FIBERMAX",(long long int)FIBERMAX,"last fiber (included)");
    harp::fits::key_write(fp,"NPARAMS",(long long int)nparams_all,"number of PSF parameters");
    harp::fits::key_write(fp,"LEGDEG",(long long int)(ncoeff-1),"degree of Legendre pol.(wave) for parameters");
    harp::fits::key_write(fp,"PALEGDEG",(long long int)(ncoeff_param_fit-1),"fitted degree of Legendre pol. for parameters");
    harp::fits::key_write(fp,"TXLEGDEG",(long long int)(ncoeff_xtrace_fit-1),"fitted degree of Legendre pol. for trace along x");
    harp::fits::key_write(fp,"TYLEGDEG",(long long int)(ncoeff_ytrace_fit-1),"fitted degree of Legendre pol. for trace along y");    
    harp::fits::key_write(fp,"GHDEGX",(long long int)GHDEGX,"degree of Hermite polynomial along CCD columns");
    harp::fits::key_write(fp,"GHDEGY",(long long int)GHDEGY,"degree of Hermite polynomial along CCD rows");
    
    // add chi2
    harp::fits::key_write(fp,"PSFERROR",psf.psf_error,"assumed PSF fractional error in chi2");
    harp::fits::key_write(fp,"READNOIS",psf.readout_noise,"assumed read out noise in chi2");
    harp::fits::key_write(fp,"GAIN",psf.gain,"assumed gain in chi2");
    
    for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf.ParamsOfBundles.begin();
	bundle_it != psf.ParamsOfBundles.end(); ++bundle_it) {
      
      const specex::PSF_Params & params_of_bundle = bundle_it->second;
      
      
      int ndf = params_of_bundle.ndata - params_of_bundle.nparams;
      double chi2pdf = 0;
      if(ndf>0) chi2pdf = params_of_bundle.chi2/ndf;
      
      char key[20];
      char comment[800];

      
      sprintf(key,"B%02dRCHI2",params_of_bundle.bundle_id);
      sprintf(comment,"best fit chi2/ndf for fiber bundle %d",params_of_bundle.bundle_id);
      harp::fits::key_write(fp,key,chi2pdf,comment);
      
      sprintf(key,"B%02dNDATA",params_of_bundle.bundle_id);
      sprintf(comment,"number of pixels in fit for fiber bundle %d",params_of_bundle.bundle_id);   
      harp::fits::key_write(fp,key,(long long int)params_of_bundle.ndata,comment);
      
      sprintf(key,"B%02dNPAR",params_of_bundle.bundle_id);
      sprintf(comment,"number of parameters in fit for fiber bundle %d",params_of_bundle.bundle_id);   
      harp::fits::key_write(fp,key,(long long int)params_of_bundle.nparams,comment);
      
    }

    // writing only few key words in the previous hdu
    fits_movrel_hdu ( fp, -1, NULL, &status ); harp::fits::check ( status );
    harp::fits::key_write(fp,"PSFTYPE","GAUSS-HERMITE","");
    harp::fits::key_write(fp,"PSFVER",PSFVER,"");
    fits_write_comment(fp,"PSF generated by specex, https://github.com/julienguy/specex",&status); harp::fits::check ( status );
    fits_write_comment(fp,"Data are in the form a binary table in the next hdu",&status); harp::fits::check ( status );   

    
    
  } // end of write key words
  
  if( spots != NULL ) { // write spots
    
    fits_movrel_hdu ( fp, 1, NULL, &status ); harp::fits::check ( status );

    specex::FitsTable table;
    { // column description
      table.AddColumnDescription("X","D","","");
      table.AddColumnDescription("Y","D","","");
      table.AddColumnDescription("FLUX","D","","");
      table.AddColumnDescription("EFLUX","D","","");
      table.AddColumnDescription("WAVE","D","","");
      table.AddColumnDescription("FIBER","D","","");
      table.AddColumnDescription("BUNDLE","D","","");
      table.AddColumnDescription("CHI2","D","","");
      table.AddColumnDescription("STATUS","D","","");
    }
    
    for(size_t s=0;s<spots->size();s++) {
      const specex::Spot_p spot=(*spots)[s];      
      std::vector<specex::FitsTableEntry> row;
      {specex::FitsTableEntry entry; entry.double_vals.resize(1); entry.double_vals[0] = spot->xc; row.push_back(entry);}
      {specex::FitsTableEntry entry; entry.double_vals.resize(1); entry.double_vals[0] = spot->yc; row.push_back(entry);}
      {specex::FitsTableEntry entry; entry.double_vals.resize(1); entry.double_vals[0] = spot->flux; row.push_back(entry);}
      {specex::FitsTableEntry entry; entry.double_vals.resize(1); entry.double_vals[0] = spot->eflux; row.push_back(entry);}
      {specex::FitsTableEntry entry; entry.double_vals.resize(1); entry.double_vals[0] = spot->wavelength; row.push_back(entry);}
      {specex::FitsTableEntry entry; entry.double_vals.resize(1); entry.double_vals[0] = spot->fiber; row.push_back(entry);}
      {specex::FitsTableEntry entry; entry.double_vals.resize(1); entry.double_vals[0] = spot->fiber_bundle; row.push_back(entry);}
      {specex::FitsTableEntry entry; entry.double_vals.resize(1); entry.double_vals[0] = spot->chi2; row.push_back(entry);}
      {specex::FitsTableEntry entry; entry.double_vals.resize(1); entry.double_vals[0] = spot->status; row.push_back(entry);}
      table.data.push_back(row);
    }
    table.Write(fp);
    
    harp::fits::key_write(fp,"EXTNAME","SPOTS","");
    fits_write_comment(fp,"Fitted spots properties",&status); harp::fits::check ( status );
    
  }
} // end of routine





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
  std::map<std::string,harp::vector_double > param_coeff;
  for(int i=0;i<table.data.size();i++) { 
    std::string pname=table.data[i][param_col].string_val;
    boost::trim(pname);
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
      param_coeff[pname] = ublas::project(param_coeff[pname],ublas::range(new_fiber_begin*ncoef_per_fiber,new_fiber_end*ncoef_per_fiber));
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
    harp::vector_double w(npts);
    harp::vector_double x(npts);
    harp::vector_double y(npts);
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
      harp::matrix_double A(npar,npar);
      harp::vector_double B(npar);
      A *= 0;
      B *= 0;
      for(int fiber=bundle_fibermin;fiber<=bundle_fibermax;fiber++) {
	const specex::Trace& trace = psf->FiberTraces[fiber];
	specex::Legendre1DPol fiberpol(legdeg,param_wavemin[pname],param_wavemax[pname]);
	for(int i=0;i<=legdeg;i++)
	  fiberpol.coeff[i]=param_coeff[pname][i+(fiber-fibermin)*(legdeg+1)];
	for(double wave=wavemin;wave<=wavemax;wave+=(wavemax-wavemin)/(legdeg_param_fit+2)) {
	  double x    = trace.X_vs_W.Value(wave);
	  double pval = fiberpol.Value(wave);
	  harp::vector_double der = pol->Monomials(x,wave);
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


