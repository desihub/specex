
void write_gauss_hermite_two_psf_fits_version_2(const specex::GaussHermite2PSF& psf, fitsfile* fp, int first_hdu) {
  SPECEX_DEBUG("write_gauss_hermite_two_psf_fits_version_2");
  
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
  
  


  int GHDEGX = psf.core_degree;
  int GHDEGY = psf.core_degree;
  int GHDEGX2 = psf.second_degree;
  int GHDEGY2 = psf.second_degree;
  

  int status = 0;
  fits_movabs_hdu ( fp, first_hdu, NULL, &status ); harp::fits::check ( status );
  
 
  
  // count parameters
  int nparams=psf.LocalNAllPar();
  int nparams_all = nparams;
  nparams_all += 2; // X and Y
  nparams_all += 1; // GH trivial order zero
#ifdef CONTINUUM
  nparams_all += 1; // continuum
#endif
  
  // get largest legendre degree of all params in all bundles and check param numbers match !!
  int ncoeff_max=0;
  for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf.ParamsOfBundles.begin();
      bundle_it != psf.ParamsOfBundles.end(); ++bundle_it) {

    if( int(bundle_it->second.AllParPolXW.size()) != nparams ) SPECEX_ERROR("Fatal inconsistency in number of parameters in psf between bundles: AllParPolXW.size=" << bundle_it->second.AllParPolXW.size() << " psf.LocalNAllPar()=" << nparams);

    for(size_t p=0;p<bundle_it->second.AllParPolXW.size();p++) {
      ncoeff_max=max(ncoeff_max,bundle_it->second.AllParPolXW[p]->ydeg+1);
    }
  }
  for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf.ParamsOfBundles.begin();
      bundle_it != psf.ParamsOfBundles.end(); ++bundle_it) {
    const specex::PSF_Params & params_of_bundle = bundle_it->second;
    for(int fiber=params_of_bundle.fiber_min; fiber<=params_of_bundle.fiber_max; fiber++) {
      ncoeff_max=max(ncoeff_max,int(psf.FiberTraces.find(fiber)->second.X_vs_W.coeff.size()));
      ncoeff_max=max(ncoeff_max,int(psf.FiberTraces.find(fiber)->second.Y_vs_W.coeff.size()));
    }
  }
  
  int delta_deg = 1; // increase degree of legendre polynomials to minimize mapping errors
  ncoeff_max += delta_deg; 
  SPECEX_DEBUG("ncoeff_max = " << ncoeff_max);
  
  
  

  specex::FitsTable table;

  { // column description
    table.AddColumnDescription("PARAM","8A","","");
    table.AddColumnDescription("WAVEMIN","D","","");
    table.AddColumnDescription("WAVEMAX","D","","");
    table.AddColumnDescription("NCOEFF","D","","");
    char coeff_tform[100];
    sprintf(coeff_tform,"%dD",ncoeff_max*NFIBERS);
    vector <int> dim; dim.resize(2);
    dim[0]=ncoeff_max;
    dim[1]=NFIBERS;
    
    string sdim = table.encode_dimension(dim);
    // debug
    vector <int> dim2 = table.decode_dimension(sdim);
    if (dim2.size()!=2 || dim2[0]!= dim[0]|| dim2[1]!=dim[1]) {
      SPECEX_ERROR("error encode/decode dimension");
    }
#ifdef CONTINUUM
    table.AddColumnDescription("COEFF",coeff_tform,sdim,"");
#endif
  }
  
 
  // get the max range of wavelength and convert to int 
  int LEGWMIN=1000000;
  int LEGWMAX=0;
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

  
  vector<string> keys;
  
  int number_of_non_zero_coeffs=0;
  harp::vector_double coeff;
  
  // first deal with X and Y
  coeff=coeffs_from_trace_x_vs_w(wavemin,wavemax,psf,ncoeff_max,NFIBERS,delta_deg,number_of_non_zero_coeffs);
  AddRow(table,"X",LEGWMIN,LEGWMAX,number_of_non_zero_coeffs,coeff);
  
  coeff=coeffs_from_trace_y_vs_w(wavemin,wavemax,psf,ncoeff_max,NFIBERS,delta_deg,number_of_non_zero_coeffs);
  AddRow(table,"Y",LEGWMIN,LEGWMAX,number_of_non_zero_coeffs,coeff);
  
  
  bool need_to_add_first_gh = true;

  // loop on all params, insert GH00
  for(int p=0;p<nparams;p++) { 
    string pname = psf.ParamName(p);
    
    if(need_to_add_first_gh && pname.find("GH-")<pname.npos) { // insert now GH param 0
      harp::vector_double  gh200coeff = coeffs_from_pold2d(psf.ParamIndex("GH2-0-0"),wavemin,wavemax,psf,ncoeff_max,NFIBERS,delta_deg,number_of_non_zero_coeffs);
      harp::vector_double  gh100coeff = -gh200coeff;
      for(int fiber=0;fiber<NFIBERS;fiber++)
	gh100coeff(fiber*ncoeff_max) += 1;
      AddRow(table,"GH-0-0",LEGWMIN,LEGWMAX,number_of_non_zero_coeffs,gh100coeff);
      need_to_add_first_gh=false;
    }
    coeff=coeffs_from_pold2d(p,wavemin,wavemax,psf,ncoeff_max,NFIBERS,delta_deg,number_of_non_zero_coeffs);
    AddRow(table,pname,LEGWMIN,LEGWMAX,number_of_non_zero_coeffs,coeff);
  }
  
#ifdef CONTINUUM
  coeff=coeffs_from_continuum(wavemin,wavemax,psf,ncoeff_max,NFIBERS,number_of_non_zero_coeffs);
  AddRow(table,"CONT",LEGWMIN,LEGWMAX,number_of_non_zero_coeffs,coeff);
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
    fits_write_comment(fp,"The size of the vector is ((FIBERMAX-FIBERMIN+1)*(NCOEFMAX))",&status); harp::fits::check ( status );
    fits_write_comment(fp,"Description of  the NPARAMS parameters : ",&status); harp::fits::check ( status );
    fits_write_comment(fp,"X        : CCD column coordinate (as a function of fiber and wavelength)",&status); harp::fits::check ( status );
    fits_write_comment(fp,"Y        : CCD row coordinate (as a function of fiber and wavelength)",&status); harp::fits::check ( status ); 
    fits_write_comment(fp,"         (X,Y)=(0,0) means that PSF is centered on center of first pixel",&status); harp::fits::check ( status );
    fits_write_comment(fp,"GHSIGX   : Sigma of first Gaussian along CCD columns for PSF core",&status); harp::fits::check ( status );
    fits_write_comment(fp,"GHSIGY   : Sigma of first Gaussian along CCD rows for PSF core",&status); harp::fits::check ( status );
    fits_write_comment(fp,"GHNSIG   : NxSigma cutoff for first Gaussian (based oncenter of pixel)",&status); harp::fits::check ( status );
    fits_write_comment(fp,"GHSIGX2  : Sigma of second Gaussian along CCD columns for PSF wings",&status); harp::fits::check ( status );
    fits_write_comment(fp,"GHSIGY2  : Sigma of second Gaussian along CCD rows for PSF wings",&status); harp::fits::check ( status );
    fits_write_comment(fp,"GH-i-j   : Hermite pol. coefficents, i along columns, j along rows,",&status); harp::fits::check ( status );
    fits_write_comment(fp,"         i is integer from 0 to GHDEGX, j is integer from 0 to GHDEGY,",&status); harp::fits::check ( status );
    fits_write_comment(fp,"         there are (GHDEGX+1)*(GHDEGY+1) such coefficents.",&status); harp::fits::check ( status );
    fits_write_comment(fp,"GH2-i-j  : Hermite pol. coefficents for sec. GH psf, i along columns, j along rows,",&status); harp::fits::check ( status );
    fits_write_comment(fp,"TAILAMP  : Amplitude of PSF tail",&status); harp::fits::check ( status );
    fits_write_comment(fp,"TAILCORE : Size in pixels of PSF tail saturation in PSF core",&status); harp::fits::check ( status );
    fits_write_comment(fp,"TAILXSCA : Scaling apply to CCD coordinate along columns for PSF tail",&status); harp::fits::check ( status );
    fits_write_comment(fp,"TAILYSCA : Scaling apply to CCD coordinate along rows for PSF tail",&status); harp::fits::check ( status );
    fits_write_comment(fp,"TAILINDE : Asymptotic power law index of PSF tail",&status); harp::fits::check ( status );
    fits_write_comment(fp,"CONT     : Continuum flux in arc image (not part of PSF)",&status); harp::fits::check ( status );
    fits_write_comment(fp,"-  ",&status); harp::fits::check ( status );
    fits_write_comment(fp,"PSF_core(X,Y) = SUM_ij (GH-i-j)*HERM(i,X/GHSIGX)*HERM(j,Y/GHSIGY)",&status); harp::fits::check ( status );
    fits_write_comment(fp,"               * GAUS(X,GHSIGX)*GAUS(Y,GHSIGY)",&status); harp::fits::check ( status );
    fits_write_comment(fp,"              +SUM_ij (GH2-i-j)*HERM(i,X/GHSIGX2)*HERM(j,Y/GHSIGY2)",&status); harp::fits::check ( status );
    fits_write_comment(fp,"                * GAUS(X,GHSIGX2)*GAUS(Y,GHSIGY2)",&status); harp::fits::check ( status );
    fits_write_comment(fp,"-  ",&status); harp::fits::check ( status );
    fits_write_comment(fp,"PSF_tail(X,Y) = TAILAMP*R^2/(TAILCORE^2+R^2)^(1+TAILINDE/2)",&status); harp::fits::check ( status );
    fits_write_comment(fp,"                with R^2=(X*TAILXSCA)^2+(Y*TAILYSCA)^2",&status); harp::fits::check ( status );
    fits_write_comment(fp,"-  ",&status); harp::fits::check ( status );
    fits_write_comment(fp,"PSF_core is integrated in pixel",&status); harp::fits::check ( status );
    fits_write_comment(fp,"PSF_tail is not, it is evaluated at center of pixel",&status); harp::fits::check ( status );
    fits_write_comment(fp,"------------------------------------------------------------------------",&status); harp::fits::check ( status );   
    harp::fits::key_write(fp,"EXTNAME","PSF","");
    harp::fits::key_write(fp,"PSFTYPE","GAUSS-HERMITE2","");
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
    harp::fits::key_write(fp,"NCOEFMAX",(long long int)(ncoeff_max),"largest value of NCOEF in the table");
    harp::fits::key_write(fp,"GHDEGX",(long long int)GHDEGX,"degree of Hermite polynomial along CCD columns");
    harp::fits::key_write(fp,"GHDEGY",(long long int)GHDEGY,"degree of Hermite polynomial along CCD rows");
    harp::fits::key_write(fp,"GHDEGX2",(long long int)GHDEGX2,"degree of Hermite polynomial along CCD columns (sec. term)");
    harp::fits::key_write(fp,"GHDEGY2",(long long int)GHDEGY2,"degree of Hermite polynomial along CCD rows (sec. term)");

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
      
      double chi2ndata_core = 0;
      if(params_of_bundle.ndata_in_core>0)
	chi2ndata_core = params_of_bundle.chi2_in_core/params_of_bundle.ndata_in_core;

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
      
      sprintf(key,"B%02dCCHI2",params_of_bundle.bundle_id);
      sprintf(comment,"chi2/ndata for fiber bundle %d in 5x5 pix. core",params_of_bundle.bundle_id);
      harp::fits::key_write(fp,key,chi2ndata_core,comment);
      
    }

 
    // writing only few key words in the previous hdu
    fits_movrel_hdu ( fp, -1, NULL, &status ); harp::fits::check ( status );
    
  {
    harp::fits::key_write(fp,"EXTNAME","PSF","");
    harp::fits::key_write(fp,"PSFTYPE","GAUSS-HERMITE2","");
    harp::fits::key_write(fp,"PSFVER",PSFVER,"");
    fits_write_comment(fp,"PSF generated by specex, https://github.com/julienguy/specex",&status); harp::fits::check ( status );
    fits_write_comment(fp,"Data are in the form a binary table in the next hdu",&status); harp::fits::check ( status );   
  }




  } // end of write keywords


  
} // end of routine

