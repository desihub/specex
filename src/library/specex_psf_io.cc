#include <fstream>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>

#include <harp.hpp>

#include <specex_psf.h>
#include <specex_psf_io.h>
#include <specex_message.h>
#include <specex_image_data.h>
#include <specex_fits.h>

#include <specex_gauss_hermite_psf.h>
#include <specex_gauss_hermite_two_psf.h>
#include <specex_hat_hermite_psf.h>

// also need #included <specex_serialisation.h>  before main() 

void specex::write_psf_fits_image(const specex::PSF_p psf, const string& filename, const int fiber, const double& wavelength, int bundle, int oversampling) {
  
  double x=psf->Xccd(fiber,wavelength);
  double y=psf->Yccd(fiber,wavelength);

  x = int(x);
  y = int(y);

  harp::vector_double P=psf->AllLocalParamsFW(fiber,wavelength,bundle);
  if(0) {
  for(int i=5;i<P.size();i++) P(i)=0; 
  P(0)=0.05;
  P(3)=0.0;
  SPECEX_WARNING("Hacking the parameters for test purpose!");
  }
  SPECEX_INFO("PSF Params " << P);
  SPECEX_INFO("PSF value at center = " << psf->PSFValueWithParamsXY(int(x),int(y),int(x),int(y),P,0,0));
  
  int nx = 2*psf->hSizeX*oversampling+1;
  int ny = 2*psf->hSizeY*oversampling+1;
  

  double sum=0;
  specex::image_data img(nx,ny);
  for(int j=0;j<ny;j++) {
    for(int i=0;i<nx;i++) {
      
      int ib = (i-nx/2)/oversampling;
      int jb = (j-ny/2)/oversampling;
      double dx = (i-nx/2)/double(oversampling)-ib;
      double dy = (j-ny/2)/double(oversampling)-jb;
      
      img(i,j)=psf->PSFValueWithParamsXY(x-dx,y-dy,ib+int(x),jb+int(y),P,0,0);
      sum += img(i,j);
    }
  }
  SPECEX_INFO("PSF integral in image = " << sum/(oversampling*oversampling));
 
  
  // get maximum of psf profile numerically
  {
    double maxval=0;
    int imax=0;
    int jmax=0;
   
    for(int j=int(y)-3;j<=int(y+3);j++)
      for(int i=int(x)-3;i<=int(x+3);i++)
	{
	  double val=psf->PSFValueWithParamsXY(x,y,i,j,P,0,0);
	  if(val>maxval) {maxval=val; imax=i; jmax=j;}
	}
    SPECEX_INFO("for x,y=" << x << "," << y << " max at i,j=" << imax << "," << jmax);
  }
  
  specex::write_new_fits_image(filename,img);
}



void specex::write_psf_xml(const specex::PSF_p psf, const std::string& filename) {
  
  std::ofstream os(filename.c_str());
  boost::archive::xml_oarchive xml_oa ( os );

  xml_oa << BOOST_SERIALIZATION_NVP(psf);
  
  os.close();
  
  SPECEX_INFO("wrote psf in " << filename);
}

void specex::read_psf_xml(specex::PSF_p& psf, const std::string& filename) {
  
  
  std::ifstream is(filename.c_str());
  
  boost::archive::xml_iarchive xml_ia ( is );

  xml_ia >> BOOST_SERIALIZATION_NVP(psf);
  
  is.close();

  SPECEX_INFO("read psf in " << filename);
  
}


void write_gauss_hermite_psf_fits_version_1(const specex::GaussHermitePSF& psf, fitsfile* fp, int first_hdu);
void write_gauss_hermite_psf_fits_version_2(const specex::GaussHermitePSF& psf, fitsfile* fp, int first_hdu);
void write_gauss_hermite_two_psf_fits_version_2(const specex::GaussHermite2PSF& psf, fitsfile* fp, int first_hdu);
void write_gauss_hermite_two_psf_fits_version_1(const specex::GaussHermite2PSF& psf, fitsfile* fp, int first_hdu);


void specex::write_psf_fits(const specex::PSF_p psf, const string& path) {
  fitsfile * fp;  
  harp::fits::create(fp,path);
  write_psf_fits(psf,fp,1);
  harp::fits::close(fp);
  SPECEX_INFO("wrote PSF in " << path);
}

void specex::write_psf_fits(const specex::PSF_p psf, fitsfile* fp, int first_hdu) {

  if(psf->Name()=="GaussHermitePSF") 
    return write_gauss_hermite_psf_fits_version_2((const specex::GaussHermitePSF&)*psf,fp,first_hdu);
  //return write_gauss_hermite_psf_fits_version_1((const specex::GaussHermitePSF&)*psf,fp,first_hdu);
  if(psf->Name()=="GaussHermite2PSF") 
    return write_gauss_hermite_two_psf_fits_version_2((const specex::GaussHermite2PSF&)*psf,fp,first_hdu);
  //return write_gauss_hermite_two_psf_fits_version_1((const specex::GaussHermite2PSF&)*psf,fp,first_hdu);
  
  SPECEX_ERROR("specex::write_psf_fits not implemented for PSF '" << psf->Name() << '"');
}

void specex::read_psf_fits(specex::PSF_p& psf, const string& filename) {
  SPECEX_ERROR("not implemented");
}

void specex::read_psf_fits(specex::PSF_p& psf, fitsfile* fp, int first_hdu) {
  SPECEX_ERROR("not implemented");
}


///////////////////////////////////


///////////////////////////////////
static void AddRow(specex::FitsTable& table,const string& PARAM, double wavemin, double wavemax, int ncoef, harp::vector_double& coeff) {
  std::vector<specex::FitsTableEntry> row;
  {specex::FitsTableEntry entry; entry.string_val = PARAM; row.push_back(entry);}
  {specex::FitsTableEntry entry; entry.double_vals.resize(1); entry.double_vals[0] = wavemin; row.push_back(entry);}
  {specex::FitsTableEntry entry; entry.double_vals.resize(1); entry.double_vals[0] = wavemax; row.push_back(entry);}
  {specex::FitsTableEntry entry; entry.double_vals.resize(1); entry.double_vals[0] = ncoef; row.push_back(entry);}
  {specex::FitsTableEntry entry; entry.double_vals = coeff; row.push_back(entry);}
  table.data.push_back(row);
}

static void AddRow1(specex::FitsTable& table,const string& PARAM, double wavemin, double wavemax, harp::vector_double& coeff) {
  std::vector<specex::FitsTableEntry> row;
  {specex::FitsTableEntry entry; entry.string_val = PARAM; row.push_back(entry);}
  {specex::FitsTableEntry entry; entry.double_vals.resize(1); entry.double_vals[0] = wavemin; row.push_back(entry);}
  {specex::FitsTableEntry entry; entry.double_vals.resize(1); entry.double_vals[0] = wavemax; row.push_back(entry);}
  {specex::FitsTableEntry entry; entry.double_vals = coeff; row.push_back(entry);}
  table.data.push_back(row);
}

harp::vector_double sample_wavelength(int ncoeff, const double& wavemin, const double& wavemax) {
  harp::vector_double wave(ncoeff);
  double wavestep = (wavemax-wavemin)/(ncoeff-1);
  for(int w=0;w<ncoeff;w++) {
    wave[w]   = wavemin + wavestep*w;
  }
  return wave;
}

harp::vector_double coeffs_from_trace_x_vs_w(const double& wavemin, const double& wavemax, const specex::PSF& psf, int ncoeff_max, int NFIBERS, int delta_deg, int& ncoeff) {
  
  // get ncoef for this param
  ncoeff=0;
  for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf.ParamsOfBundles.begin();
      bundle_it != psf.ParamsOfBundles.end(); ++bundle_it) {
    const specex::PSF_Params & params_of_bundle = bundle_it->second;
    for(int fiber=params_of_bundle.fiber_min; fiber<=params_of_bundle.fiber_max; fiber++) {
      ncoeff=max(ncoeff,int(psf.FiberTraces.find(fiber)->second.X_vs_W.coeff.size()));
    }
  }
  ncoeff += delta_deg; 
  
  harp::vector_double wave = sample_wavelength(ncoeff,wavemin,wavemax);
  harp::vector_double values(ncoeff); 
  harp::vector_double coeff(ncoeff_max*NFIBERS);
  coeff.clear();
  int fiber_index=0;
  for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf.ParamsOfBundles.begin();
      bundle_it != psf.ParamsOfBundles.end(); ++bundle_it) {
    const specex::PSF_Params & params_of_bundle = bundle_it->second;
    for(int fiber=params_of_bundle.fiber_min; fiber<=params_of_bundle.fiber_max; fiber++,fiber_index++) {
      const specex::Trace &trace = psf.FiberTraces.find(fiber)->second;
      specex::Legendre1DPol pol1d(ncoeff-1,wavemin,wavemax);
      for(int w=0;w<ncoeff;w++) {
	values[w] = trace.X_vs_W.Value(wave[w]);
      }
      pol1d.Fit(wave,values,0,false);
      for(int w = 0; w < ncoeff ; w++) {
	coeff(fiber_index*ncoeff_max+w) = pol1d.coeff(w);
      }
    }
  }
  return coeff;
}

harp::vector_double coeffs_from_trace_y_vs_w(const double& wavemin, const double& wavemax, const specex::PSF& psf, int ncoeff_max, int NFIBERS, int delta_deg, int& ncoeff) {
  
  // get ncoef for this param
  ncoeff=0;
  for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf.ParamsOfBundles.begin();
      bundle_it != psf.ParamsOfBundles.end(); ++bundle_it) {
    const specex::PSF_Params & params_of_bundle = bundle_it->second;
    for(int fiber=params_of_bundle.fiber_min; fiber<=params_of_bundle.fiber_max; fiber++) {
      ncoeff=max(ncoeff,int(psf.FiberTraces.find(fiber)->second.Y_vs_W.coeff.size()));
    }
  }
  ncoeff += delta_deg; 
  
  harp::vector_double wave = sample_wavelength(ncoeff,wavemin,wavemax);
  harp::vector_double values(ncoeff); 
  harp::vector_double coeff(ncoeff_max*NFIBERS);
  coeff.clear();
  int fiber_index=0;
  for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf.ParamsOfBundles.begin();
      bundle_it != psf.ParamsOfBundles.end(); ++bundle_it) {
    const specex::PSF_Params & params_of_bundle = bundle_it->second;
    for(int fiber=params_of_bundle.fiber_min; fiber<=params_of_bundle.fiber_max; fiber++,fiber_index++) {
      const specex::Trace &trace = psf.FiberTraces.find(fiber)->second;
      specex::Legendre1DPol pol1d(ncoeff-1,wavemin,wavemax);
      for(int w=0;w<ncoeff;w++) {
	values[w] = trace.Y_vs_W.Value(wave[w]);
      }
      pol1d.Fit(wave,values,0,false);
      for(int w = 0; w < ncoeff ; w++) {
	coeff(fiber_index*ncoeff_max+w) = pol1d.coeff(w);
      }
    }
  }
  return coeff;
}

#ifdef CONTINUUM
harp::vector_double coeffs_from_continuum(const double& wavemin, const double& wavemax, const specex::PSF& psf, int ncoeff_max, int NFIBERS, int& ncoeff) {
  
  
  ncoeff = 0;
  for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf.ParamsOfBundles.begin();
      bundle_it != psf.ParamsOfBundles.end(); ++bundle_it) {
    ncoeff = max(ncoeff,int(bundle_it->second.ContinuumPol.coeff.size()));
  }
  // direct mapping so no need to increase degree

  // recompute anyway because of change of wavemin,wavemax
  harp::vector_double wave = sample_wavelength(ncoeff,wavemin,wavemax);
  harp::vector_double values(ncoeff); 
  harp::vector_double coeff(ncoeff_max*NFIBERS);
  coeff.clear();
  int fiber_index=0;
  for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf.ParamsOfBundles.begin();
      bundle_it != psf.ParamsOfBundles.end(); ++bundle_it) {
    const specex::PSF_Params & params_of_bundle = bundle_it->second;
    for(int w=0;w<ncoeff;w++) {
      values[w] = params_of_bundle.ContinuumPol.Value(wave[w]);
    }
    specex::Legendre1DPol pol1d(ncoeff-1,wavemin,wavemax);
    pol1d.Fit(wave,values,0,false);
    
    // now copy parameters;
    for(int fiber=params_of_bundle.fiber_min; fiber<=params_of_bundle.fiber_max; fiber++,fiber_index++) {     
      for(int w = 0; w < ncoeff ; w++) {
	coeff(fiber_index*ncoeff_max+w) = pol1d.coeff(w); // this is the definition of the ordering, (wave,fiber)
      }
    }
  }
  return coeff;
}
#endif



harp::vector_double coeffs_from_pold2d(int param_index, const double& wavemin, const double& wavemax, const specex::PSF& psf, int ncoeff_max, int NFIBERS, int delta_deg, int& ncoeff) {

  
  // get ncoef for this param
  ncoeff=0;
  for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf.ParamsOfBundles.begin();
      bundle_it != psf.ParamsOfBundles.end(); ++bundle_it) {
    ncoeff=max(ncoeff,bundle_it->second.AllParPolXW[param_index]->ydeg+1);
  }
  ncoeff += delta_deg; 
  
  harp::vector_double wave = sample_wavelength(ncoeff,wavemin,wavemax);
  harp::vector_double values(ncoeff); 
  harp::vector_double coeff(ncoeff_max*NFIBERS);
  coeff.clear();
  int fiber_index=0;
  for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf.ParamsOfBundles.begin();
      bundle_it != psf.ParamsOfBundles.end(); ++bundle_it) {
    const specex::PSF_Params & params_of_bundle = bundle_it->second;
    const specex::Pol_p pol2d = params_of_bundle.AllParPolXW[param_index]; // this is the 2D polynomiald of x_ccd and wave for this param and bundle
    
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
	coeff(fiber_index*ncoeff_max+w) = pol1d.coeff(w); // this is the definition of the ordering, (wave,fiber)
      }
    }
  }
  
  return coeff;
}



void write_gauss_hermite_two_psf_fits_version_2(const specex::GaussHermite2PSF& psf, fitsfile* fp, int first_hdu) {
   
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

    if( int(bundle_it->second.AllParPolXW.size()) != nparams ) SPECEX_ERROR("Fatal inconsistency in number of parameters in psf between bundles");

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
  cout << "ncoeff_max = " << ncoeff_max << endl;
  
  

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
      LEGWMAX=max(LEGWMAX,int(floor(params_of_bundle.AllParPolXW[p]->ymax))+1);
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
  }
  
}


void write_gauss_hermite_two_psf_fits_version_1(const specex::GaussHermite2PSF& psf, fitsfile* fp, int first_hdu) {
   
  ////////////////////////////
  string PSFVER = "1";
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
  
  // get largest legendre degree of all params in all bundles and check param numbers match !!
  int ncoeff=0;
  int nparams=psf.LocalNAllPar();
  int nparams_all = nparams;
  nparams_all += 2; // X and Y
  nparams_all += 1; // GH trivial order zero
#ifdef CONTINUUM
  nparams_all += 1; // continuum
#endif
  
  for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf.ParamsOfBundles.begin();
      bundle_it != psf.ParamsOfBundles.end(); ++bundle_it) {

    if( int(bundle_it->second.AllParPolXW.size()) != nparams ) SPECEX_ERROR("Fatal inconsistency in number of parameters in psf between bundles");

    for(size_t p=0;p<bundle_it->second.AllParPolXW.size();p++) {
      ncoeff=max(ncoeff,bundle_it->second.AllParPolXW[p]->ydeg+1);
    }
  }

  ncoeff += 2; // increase degree of legendre polynomials to minimize mapping errors
  
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
	LEGWMAX=max(LEGWMAX,int(floor(params_of_bundle.AllParPolXW[p]->ymax))+1);
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
	    coeff(fiber_index*ncoeff+w)   =  pol1d_x.coeff(w);
	    coeff_y(fiber_index*ncoeff+w) =  pol1d_y.coeff(w);
	  }
	  
	  
	} // end of loop on fiber
      } // end of loop on bundles
      
      AddRow1(table,"X",LEGWMIN,LEGWMAX,coeff);
      AddRow1(table,"Y",LEGWMIN,LEGWMAX,coeff_y);
      
    } // end of X Y context

    // precompute coeffs of GH2-0-0
    harp::vector_double  gh200coeff;
    harp::vector_double  gh100coeff;
    {
      int gh200_index = psf.ParamIndex("GH2-0-0");
      int fiber_index=0;
      gh200coeff.resize(coeff.size());
      gh200coeff.clear();
      gh100coeff.resize(coeff.size());
      gh100coeff.clear();
      
      // duplicated code that need to be simplified
      for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf.ParamsOfBundles.begin();
	  bundle_it != psf.ParamsOfBundles.end(); ++bundle_it) {
	const specex::PSF_Params & params_of_bundle = bundle_it->second;
	const specex::Pol_p pol2d = params_of_bundle.AllParPolXW[gh200_index]; // this is the 2D polynomiald of x_ccd and wave for this param and bundle
      
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
	    gh200coeff(fiber_index*ncoeff+w) = pol1d.coeff(w); // this is the definition of the ordering, (wave,fiber)
	    gh100coeff(fiber_index*ncoeff+w) = -gh200coeff(fiber_index*ncoeff+w); // this is the definition of the ordering, (wave,fiber)
	  }
	  gh100coeff(fiber_index*ncoeff) += 1;
	}
      }
    }

    for(int p=0;p<nparams;p++) { 
      
      coeff.clear();
      
      string pname = psf.ParamName(p);

      if(need_to_add_first_gh && pname.find("GH-")<pname.npos) { // insert now GH param 0
	AddRow1(table,"GH-0-0",LEGWMIN,LEGWMAX,gh100coeff);
	need_to_add_first_gh = false;
	cout << "GH-0-0  : " << gh100coeff << endl;
	cout << "GH2-0-0 : " << gh200coeff << endl;
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
	    coeff(fiber_index*ncoeff+w) = pol1d.coeff(w); // this is the definition of the ordering, (wave,fiber)
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
	    coeff(fiber_index*ncoeff+w)   =  pol1d.coeff(w);
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
    fits_write_comment(fp,"GHNSIG   : NxSigma cutoff for first Gaussian",&status); harp::fits::check ( status );
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
    harp::fits::key_write(fp,"LEGDEG",(long long int)(ncoeff-1),"degree of Legendre pol.(wave) for parameters");
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
  }
  
}



void write_gauss_hermite_psf_fits_version_2(const specex::GaussHermitePSF& psf, fitsfile* fp, int first_hdu) {
   
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
  
  int GHDEGX = psf.Degree();
  int GHDEGY = psf.Degree();
  

  int status = 0;
  fits_movabs_hdu ( fp, first_hdu, NULL, &status ); harp::fits::check ( status );
  
  
  
  // get largest legendre degree of all params in all bundles and check param numbers match !!
  int ncoeff=0;
  int nparams=psf.LocalNAllPar();
  int nparams_all = nparams;
  nparams_all += 2; // X and Y
  nparams_all += 1; // GH trivial order zero
#ifdef CONTINUUM
  nparams_all += 1; // continuum
#endif
  
  for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf.ParamsOfBundles.begin();
      bundle_it != psf.ParamsOfBundles.end(); ++bundle_it) {

    if( int(bundle_it->second.AllParPolXW.size()) != nparams ) SPECEX_ERROR("Fatal inconsistency in number of parameters in psf between bundles");

    for(size_t p=0;p<bundle_it->second.AllParPolXW.size();p++) {
      ncoeff=max(ncoeff,bundle_it->second.AllParPolXW[p]->ydeg+1);
    }
  }
  

  

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
	LEGWMAX=max(LEGWMAX,int(floor(params_of_bundle.AllParPolXW[p]->ymax))+1);
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
	    coeff(fiber_index*ncoeff+w)   =  pol1d_x.coeff(w);
	    coeff_y(fiber_index*ncoeff+w) =  pol1d_y.coeff(w);
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
	
	for(int fiber_index=0;fiber_index<NFIBERS;fiber_index++) 
	  coeff(fiber_index*ncoeff)=1; // need check
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
	    coeff(fiber_index*ncoeff+w) = pol1d.coeff(w); // this is the definition of the ordering, (wave,fiber)
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
	    coeff(fiber_index*ncoeff+w)   =  pol1d.coeff(w);
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
    fits_write_comment(fp,"GHSIGX2  : Sigma of second Gaussian along CCD columns for PSF wings",&status); harp::fits::check ( status );
    fits_write_comment(fp,"GHSIGY2  : Sigma of second Gaussian along CCD rows for PSF wings",&status); harp::fits::check ( status );
    fits_write_comment(fp,"GHSCAL2  : Amplitude of second Gaussian; amp. of first is (1-GHSCAL2)",&status); harp::fits::check ( status );
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
    fits_write_comment(fp,"              * [ (1-GHSCAL2)*GAUS(X,GHSIGX)*GAUS(Y,GHSIGY)",&status); harp::fits::check ( status );
    fits_write_comment(fp,"                 + GHSCAL2*GAUS(X,GHSIGX2)*GAUS(Y,GHSIGY2) ]",&status); harp::fits::check ( status );
    fits_write_comment(fp,"-  ",&status); harp::fits::check ( status );
    fits_write_comment(fp,"PSF_tail(X,Y) = TAILAMP*R^2/(TAILCORE^2+R^2)^(1+TAILINDE/2)",&status); harp::fits::check ( status );
    fits_write_comment(fp,"                with R^2=(X/TAILXSCA)^2+(Y/TAILYSCA)^2",&status); harp::fits::check ( status );
    fits_write_comment(fp,"-  ",&status); harp::fits::check ( status );
    fits_write_comment(fp,"PSF_core is integrated in pixel",&status); harp::fits::check ( status );
    fits_write_comment(fp,"PSF_tail is not, it is evaluated at center of pixel",&status); harp::fits::check ( status );
    fits_write_comment(fp,"------------------------------------------------------------------------",&status); harp::fits::check ( status );   
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
  }
  
}












/// OLDER VERSIONS ///


















void write_gauss_hermite_psf_fits_version_1(const specex::GaussHermitePSF& psf, fitsfile* fp, int first_hdu) {
   
  ////////////////////////////
  string PSFVER = "1";
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
    FIBERMAX = min(FIBERMAX,bundle_it->second.fiber_max);
  
    NFIBERS += (bundle_it->second.fiber_max-bundle_it->second.fiber_min+1);
  }
  
  int status = 0;
  fits_movabs_hdu ( fp, first_hdu, NULL, &status ); harp::fits::check ( status );
  
  
  
  // get largest legendre degree of all params in all bundles and check param numbers match !!
  int ncoeff=0;
  int nparams=psf.LocalNAllPar();
  
  for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf.ParamsOfBundles.begin();
      bundle_it != psf.ParamsOfBundles.end(); ++bundle_it) {

    if( int(bundle_it->second.AllParPolXW.size()) != nparams ) SPECEX_ERROR("Fatal inconsistency in number of parameters in psf between bundles");

    for(size_t p=0;p<bundle_it->second.AllParPolXW.size();p++) {
      ncoeff=max(ncoeff,bundle_it->second.AllParPolXW[p]->ydeg+1);
    }
  }
  
  specex::image_data image;
  vector<string> keys;

  int LEGWMIN=1000000;
  int LEGWMAX=0;

  int nparams_all = nparams;
  nparams_all += 1; // GH trivial order zero
#ifdef EXTERNAL_TAIL
  nparams_all += 5; // tail params
#endif

  { // data
    image = specex::image_data(ncoeff*(nparams_all),NFIBERS);
    image.data.clear();
    
    // now loop on real psf parameters

    harp::vector_double wave(ncoeff);
    harp::vector_double values(ncoeff);
    
    // get the max range of wavelength and convert to int 
    
    for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf.ParamsOfBundles.begin();
	  bundle_it != psf.ParamsOfBundles.end(); ++bundle_it) {
      const specex::PSF_Params & params_of_bundle = bundle_it->second;
      for(int p=0;p<nparams;p++) { 
	LEGWMIN=min(LEGWMIN,int(floor(params_of_bundle.AllParPolXW[p]->ymin)));
	LEGWMAX=max(LEGWMAX,int(floor(params_of_bundle.AllParPolXW[p]->ymax))+1);
      }
    }
    
    
    double wavemin = double(LEGWMIN);
    double wavemax = double(LEGWMAX);

    bool need_to_add_first_gh = true;
    for(int p=0;p<nparams;p++) { 
      string pname = psf.ParamName(p);
      if(need_to_add_first_gh && pname.find("GH-")<pname.npos) { // insert now GH param 0
	// first gauss-hermite param
	
	keys.push_back("GH-0-0"); for(int f=0;f<NFIBERS;f++) image(0,f)=1;	
	need_to_add_first_gh = false;
      }

      keys.push_back(pname);
      
      int fiber_index=0;
      for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf.ParamsOfBundles.begin();
	  bundle_it != psf.ParamsOfBundles.end(); ++bundle_it) {
	const specex::PSF_Params & params_of_bundle = bundle_it->second;
	
	const specex::Pol_p pol2d = params_of_bundle.AllParPolXW[p]; // this is the 2D polynomiald of x_ccd and wave for this param and bundle
	

	
	for(int fiber=params_of_bundle.fiber_min; fiber<=params_of_bundle.fiber_max; fiber++,fiber_index++) {
	  
	  const specex::Trace& trace = psf.FiberTraces.find(fiber)->second; // X_vs_W.Value();

	  // build a Legendre1DPol out of the Legendre2DPol
	  specex::Legendre1DPol pol1d(ncoeff-1,wavemin,wavemax);
	  double wavestep = (wavemax-wavemin)/(ncoeff-1);
	  for(int w=0;w<ncoeff;w++) {
	    wave[w]   = wavemin + wavestep*w;
	    values[w] = pol2d->Value(trace.X_vs_W.Value(wave[w]),wave[w]);
	  }
	  pol1d.Fit(wave,values,0,false);

	  // now copy parameters;
	  
	  for(int w = 0; w < ncoeff ; w++) {
	    image((p+1)*ncoeff+w,fiber_index) = pol1d.coeff(w); // this is the definition of the ordering, (wave,fiber)
	  }
	  
	} // end of loop on fibers of bundle
	      
      } // end of loop on bundles
    } // end of loop on params
    
  }
  
  // write image
  {
    harp::fits::img_append < double > ( fp, image.n_rows(), image.n_cols() );
    harp::fits::img_write ( fp, image.data );
  }
  
  // write keywords
  {
    fits_write_comment(fp,"PSF generated by specex, https://github.com/julienguy/specex",&status); harp::fits::check ( status );
    
    {
      char date_comment[80];
      time_t t = time(0);   // get time now
      struct tm * now = localtime( & t );
      sprintf(date_comment,"PSF fit date %04d-%02d-%02d",(now->tm_year + 1900),(now->tm_mon + 1),now->tm_mday);
      fits_write_comment(fp,date_comment,&status); harp::fits::check ( status );
    }
    ///////////////////////xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx////////
    fits_write_comment(fp,"Each row of the image contains the PSF parameters of a fiber in the form",&status); harp::fits::check ( status );
    fits_write_comment(fp,"of Legendre coefficients. The coefficients of a given parameter are",&status); harp::fits::check ( status );
    fits_write_comment(fp,"contiguous in a row of the image. ",&status); harp::fits::check ( status );
    fits_write_comment(fp,"LEGDEG gives the degree of Legendre pol (LEGDEG+1 coeffs. per parameter)",&status); harp::fits::check ( status );
    fits_write_comment(fp,"LEGWMIN and LEGWMAX are needed to compute reduced wavelength.",&status); harp::fits::check ( status );
    fits_write_comment(fp,"PXXX(wave) = Legendre(2*(wave-LEGWMIN)/(LEGWMAX-LEGWMIN)-1,COEFFXXX)" ,&status); harp::fits::check ( status );
    fits_write_comment(fp,"The definition of the parameter XXX is given by the key word PXXX.",&status); harp::fits::check ( status );
    fits_write_comment(fp,"The number of parameters is given by the key word NPARAMS,",&status); harp::fits::check ( status );
    fits_write_comment(fp,"and the number of fibers = NAXIS2 = FIBERMAX-FIBERMIN+1",&status); harp::fits::check ( status );
    
    
    harp::fits::key_write(fp,"MJD",(long long int)psf.mjd,"MJD of arc lamp exposure");
    harp::fits::key_write(fp,"PLATEID",psf.plate_id,"plate ID of arc lamp exposure");
    harp::fits::key_write(fp,"CAMERA",psf.camera_id,"camera ID");
    
    harp::fits::key_write(fp,"PSFTYPE","GAUSS-HERMITE","");
    harp::fits::key_write(fp,"PSFVER",PSFVER,"");
    harp::fits::key_write(fp,"ARCEXP",psf.arc_exposure_id,"ID of arc lamp exposure used to fit PSF");
    
    harp::fits::key_write(fp,"NPIX_X",(long long int)NPIX_X,"number of columns in input CCD image");
    harp::fits::key_write(fp,"NPIX_Y",(long long int)NPIX_Y,"number of rows in input CCD image");
    harp::fits::key_write(fp,"BUNDLMIN",(long long int)BUNDLMIN,"first bundle of fibers (starting at 0)");
    harp::fits::key_write(fp,"BUNDLMAX",(long long int)BUNDLMAX,"last bundle of fibers (included)");
    harp::fits::key_write(fp,"FIBERMIN",(long long int)FIBERMIN,"first fiber (starting at 0)");
    harp::fits::key_write(fp,"FIBERMAX",(long long int)FIBERMAX,"last fiber (included)");
    harp::fits::key_write(fp,"NPARAMS",(long long int)keys.size(),"number of PSF parameters");
    harp::fits::key_write(fp,"LEGDEG",(long long int)(ncoeff-1),"degree of Legendre pol.(wave) for parameters");
    harp::fits::key_write(fp,"LEGWMIN",(long long int)LEGWMIN,"min. wave (A) for Legendre pol.");
    harp::fits::key_write(fp,"LEGWMAX",(long long int)LEGWMAX,"max. wave (A) for Legendre pol.");
    // harp::fits::key_write(fp,"NFIBERS",(long long int)NFIBERS,"number of fibers");
    
    // write first dummy GH param
    char keyname[8];
    char comment[80];
    for(size_t k=0;k<keys.size(); k++) {
      sprintf(keyname,"P%03d",int(k));
      sprintf(comment,"Param. Leg. coeff in cols %d-%d (start. at 0)",int(k*ncoeff),int((k+1)*ncoeff-1));
      harp::fits::key_write(fp,keyname,keys[k].c_str(),comment);
    }
  }
  

}
