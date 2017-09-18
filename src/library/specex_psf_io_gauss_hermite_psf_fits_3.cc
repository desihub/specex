#include <fstream>
#include <boost/algorithm/string.hpp>
#include <harp.hpp>
#include <specex_fits.h>
#include <specex_trace.h>
#include <specex_gauss_hermite_psf.h>

using namespace std ;

static void AddRow2(specex::FitsTable& table,const string& PARAM, harp::vector_double& coeff, int legdegx, int legdegw) {
  std::vector<specex::FitsTableEntry> row;
  {specex::FitsTableEntry entry; entry.string_val = PARAM; row.push_back(entry);}
  {specex::FitsTableEntry entry; entry.double_vals = coeff; row.push_back(entry);}
  {specex::FitsTableEntry entry; entry.int_vals.resize(1); entry.int_vals[0] = legdegx; row.push_back(entry);}
  {specex::FitsTableEntry entry; entry.int_vals.resize(1); entry.int_vals[0] = legdegw; row.push_back(entry);}
  table.data.push_back(row);
}


void write_gauss_hermite_psf_fits_version_3(const specex::GaussHermitePSF& psf, fitsfile* fp) {
  
  SPECEX_INFO("write PSF");
  SPECEX_DEBUG("write_gauss_hermite_psf_fits_version_3");
  
  ////////////////////////////
  string PSFVER = "3";
  ////////////////////////////
  
  int NPIX_X = psf.ccd_image_n_cols;
  int NPIX_Y = psf.ccd_image_n_rows;
  int GHDEGX = psf.Degree();
  int GHDEGY = psf.Degree();

  // FIBERMIN FIBERMAX WAVEMIN WAVEMAX are always defined by the trace set
  // to ensure consistency
  double WAVEMIN=1000000;
  double WAVEMAX=0;
  int    FIBERMIN=100000;
  int    FIBERMAX=0;
  for(std::map<int,specex::Trace>::const_iterator it=psf.FiberTraces.begin(); it!=psf.FiberTraces.end(); it++) {
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
  
  int nparams=psf.LocalNAllPar();
  int nparams_all = nparams;
  nparams_all += 1; // GH trivial order zero
#ifdef CONTINUUM
  nparams_all += 1; // continuum
#endif
  
  int ncoeff=0;
  for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf.ParamsOfBundles.begin();
      bundle_it != psf.ParamsOfBundles.end(); ++bundle_it) {
    if( int(bundle_it->second.AllParPolXW.size()) != nparams ) SPECEX_ERROR("Fatal inconsistency in number of parameters in psf between bundles: AllParPolXW.size=" << bundle_it->second.AllParPolXW.size() << " psf.LocalNAllPar()=" << nparams);
    for(size_t p=0;p<bundle_it->second.AllParPolXW.size();p++) {
      ncoeff=max(ncoeff,bundle_it->second.AllParPolXW[p]->ydeg+1);
    }
  }
  SPECEX_DEBUG("ncoeff = " << ncoeff);

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

  
  
  harp::vector_double wave(ncoeff);
  {
    double wavestep = (WAVEMAX-WAVEMIN)/(ncoeff-1);
    for(int w=0;w<ncoeff;w++) {
      wave[w]   = WAVEMIN + wavestep*w;
    }
  }
  
  vector<string> keys;
  
  
  harp::vector_double coeff(ncoeff*NFIBERS);
  harp::vector_double values(ncoeff);
  
  bool need_to_add_first_gh = true;
  for(int p=0;p<nparams;p++) {  // loop on all PSF parameters      
        
    string pname = psf.ParamName(p);    
    if(need_to_add_first_gh && pname.find("GH-")<pname.npos) { 
      // insert now GH param 0 which is not a PSF parameter in C++ but needed but specter
      
      for(int fiber=0;fiber<NFIBERS;fiber++) 
	coeff(fiber*ncoeff)=1; 
      AddRow2(table,"GH-0-0",coeff,0,0);
      need_to_add_first_gh = false;
    }
        
    int legdegx=0;
    int legdegw=0;
    // loop on bundles
    for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf.ParamsOfBundles.begin();
	bundle_it != psf.ParamsOfBundles.end(); ++bundle_it) {
      const specex::PSF_Params & params_of_bundle = bundle_it->second;      
      const specex::Pol_p pol2d = params_of_bundle.AllParPolXW[p]; // this is the 2D polynomiald of x_ccd and wave for this param and bundle
      legdegx=max(legdegx,pol2d->xdeg);
      legdegw=max(legdegw,pol2d->ydeg);

      for(int fiber=params_of_bundle.fiber_min; fiber<=params_of_bundle.fiber_max; fiber++) {
	
	const specex::Trace& trace = psf.FiberTraces.find(fiber)->second; // X_vs_W.Value();
	
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

    AddRow2(table,pname,coeff,legdegx,legdegw); 
  } // end of loop on params

  { // add a parameter to link fibers and bundles in fit
    coeff.clear();
    for(int fiber=FIBERMIN;fiber<=FIBERMAX;fiber++)
      coeff((fiber-FIBERMIN)*ncoeff) = -1; // no bundle
    for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf.ParamsOfBundles.begin();
	bundle_it != psf.ParamsOfBundles.end(); ++bundle_it) {
      const specex::PSF_Params & params_of_bundle = bundle_it->second;
      for(int fiber=params_of_bundle.fiber_min; fiber<=params_of_bundle.fiber_max; fiber++) {
	coeff((fiber-FIBERMIN)*ncoeff) = bundle_it->first;
      }
    }
    AddRow2(table,"BUNDLE",coeff,0,0); 
  }
  { // add a parameter with the fit status
    coeff.clear();
    for(int fiber=FIBERMIN;fiber<=FIBERMAX;fiber++)
      coeff((fiber-FIBERMIN)*ncoeff) = -1; // no bundle
    for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf.ParamsOfBundles.begin();
	bundle_it != psf.ParamsOfBundles.end(); ++bundle_it) {
      const specex::PSF_Params & params_of_bundle = bundle_it->second;
      for(int fiber=params_of_bundle.fiber_min; fiber<=params_of_bundle.fiber_max; fiber++) {
	coeff((fiber-FIBERMIN)*ncoeff) = params_of_bundle.fit_status;
      }
    }
    AddRow2(table,"STATUS",coeff,0,0); 
  }
  
    
#ifdef CONTINUUM  
  { // add a parameter for continuum
    coeff.clear();
    int fiber_index=0;
    int legdegx=0;
    int legdegw=0;
    for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf.ParamsOfBundles.begin();
	bundle_it != psf.ParamsOfBundles.end(); ++bundle_it) {
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
    AddRow2(table,"CONT",coeff,legdegx,legdegw);
  }
#endif    
    
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
    harp::fits::key_write(fp,"FIBERMIN",(long long int)FIBERMIN,"first fiber (starting at 0)");
    harp::fits::key_write(fp,"FIBERMAX",(long long int)FIBERMAX,"last fiber (included)");
    harp::fits::key_write(fp,"NPARAMS",(long long int)nparams_all,"number of PSF parameters");
    harp::fits::key_write(fp,"LEGDEG",(long long int)(ncoeff-1),"degree of Legendre pol.(wave) for parameters");
    harp::fits::key_write(fp,"GHDEGX",(long long int)GHDEGX,"degree of Hermite polynomial along CCD columns");
    harp::fits::key_write(fp,"GHDEGY",(long long int)GHDEGY,"degree of Hermite polynomial along CCD rows");
    harp::fits::key_write(fp,"WAVEMIN",WAVEMIN,"minimum wavelength (A), used for the Legendre polynomials");
    harp::fits::key_write(fp,"WAVEMAX",WAVEMAX,"maximum wavelength (A), used for the Legendre polynomials");
    
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

  } // end of write key words
} // end of routine





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
  std::map<std::string,harp::vector_double > param_coeff;
  std::map<std::string,int > param_degx;
  std::map<std::string,int > param_degw;
  
  for(int i=0;i<table.data.size();i++) { 
    std::string pname=table.data[i][param_col].string_val;
    boost::trim(pname);
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
      harp::matrix_double A(npar,npar); A.clear();
      harp::vector_double B(npar); B.clear();
      
      int npoints=0;
      for(int fiber=bundle_fibermin;fiber<=bundle_fibermax;fiber++) {
	const specex::Trace& trace = psf->FiberTraces[fiber];
	specex::Legendre1DPol fiberpol(LEGDEG,WAVEMIN,WAVEMAX);
	for(int cj=0;cj<=LEGDEG;cj++)
	  fiberpol.coeff[cj]=param_coeff[pname][cj+(fiber-FIBERMIN)*(LEGDEG+1)];
	
	for(double wave=WAVEMIN;wave<WAVEMAX+0.01;wave+=(WAVEMAX-WAVEMIN)/(degw+1)) {
	  double x    = trace.X_vs_W.Value(wave);
	  double pval = fiberpol.Value(wave);
	  harp::vector_double der = pol->Monomials(x,wave);
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


