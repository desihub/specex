#include <fstream>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/algorithm/string.hpp>

#include <harp.hpp>

#include <specex_psf.h>
#include <specex_spot.h>
#include <specex_psf_io.h>
#include <specex_message.h>
#include <specex_image_data.h>
#include <specex_fits.h>

#include <specex_trace.h>

#include <specex_gauss_hermite_psf.h>
#include <specex_gauss_hermite_two_psf.h>
#include <specex_hat_hermite_psf.h>

void specex::read_psf(specex::PSF_p& psf, const std::string& filename, int first_bundle, int last_bundle) {
  
  if (filename.find(".xml") != std::string::npos) {
    SPECEX_INFO("read xml file " << filename);
    specex::read_psf_xml(psf,filename);
  } else if (filename.find(".fits") != std::string::npos) {
    SPECEX_INFO("read fits file " << filename);
    specex::read_psf_fits(psf,filename,first_bundle,last_bundle);
  } else {
    SPECEX_ERROR("not sure how to read this file (expect xxx.fits or xxx.xml) " << filename);
  }
}
void specex::write_psf(const specex::PSF_p psf, const std::string& filename, std::vector<specex::Spot_p> *spots) {
  
  if (filename.find(".xml") != std::string::npos) {
    SPECEX_INFO("write xml file " << filename);
    specex::write_psf_xml(psf,filename);
  } else if (filename.find(".fits") != std::string::npos) {
    SPECEX_INFO("write fits file " << filename);
    specex::write_psf_fits(psf,filename,spots);
  } else {
    SPECEX_ERROR("not sure how to write this file (expect xxx.fits or xxx.xml) " << filename);
  }
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
harp::vector_double coeffs_from_trace_x_vs_w(const specex::PSF& psf, int ncoeff, int NFIBERS) {
  
 harp::vector_double coeff(ncoeff*NFIBERS);
 coeff.clear();
 int fiber_index=0;
 for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf.ParamsOfBundles.begin();
     bundle_it != psf.ParamsOfBundles.end(); ++bundle_it) {
   const specex::PSF_Params & params_of_bundle = bundle_it->second;
   for(int fiber=params_of_bundle.fiber_min; fiber<=params_of_bundle.fiber_max; fiber++,fiber_index++) {
     const specex::Trace &trace = psf.FiberTraces.find(fiber)->second;
     int n=trace.X_vs_W.coeff.size();
     if(n>ncoeff) {
       SPECEX_ERROR("ncoeff in trace X_vs_W=" << n << ">" << ncoeff);
     }
     for(int w = 0; w < n ; w++) {
       coeff(fiber_index*ncoeff+w) = trace.X_vs_W.coeff(w);
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

harp::vector_double coeffs_from_trace_y_vs_w(const specex::PSF& psf, int ncoeff, int NFIBERS) {
  
 harp::vector_double coeff(ncoeff*NFIBERS);
 coeff.clear();
 int fiber_index=0;
 for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = psf.ParamsOfBundles.begin();
     bundle_it != psf.ParamsOfBundles.end(); ++bundle_it) {
   const specex::PSF_Params & params_of_bundle = bundle_it->second;
   for(int fiber=params_of_bundle.fiber_min; fiber<=params_of_bundle.fiber_max; fiber++,fiber_index++) {
     const specex::Trace &trace = psf.FiberTraces.find(fiber)->second;
     int n=trace.Y_vs_W.coeff.size();
     if(n>ncoeff) {
       SPECEX_ERROR("ncoeff in trace Y_vs_W=" << n << ">" << ncoeff);
     }
     for(int w = 0; w < n ; w++) {
       coeff(fiber_index*ncoeff+w) = trace.Y_vs_W.coeff(w);
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

void _write_trace(const specex::PSF &psf, fitsfile *fp, int hdu, bool is_x) {
  int ncoeff=0;
  double wavemin=0;
  double wavemax=0;
  for(std::map<int,specex::Trace>::const_iterator it=psf.FiberTraces.begin(); it!=psf.FiberTraces.end(); it++) {
    const specex::Legendre1DPol *pol=0;
    if(is_x) pol = &(it->second.X_vs_W);
    else pol = &(it->second.Y_vs_W);
    
    ncoeff=max(ncoeff,int(pol->coeff.size()));
    if(wavemin==0) wavemin=pol->xmin;
    else if(wavemin != pol->xmin) {SPECEX_ERROR("requires same wavemin for all traces");}
    if(wavemax==0) wavemax=pol->xmax;
    else if(wavemax != pol->xmax) {SPECEX_ERROR("requires same wavemax for all traces");}
  }
  int nfibers=psf.FiberTraces.size();
    
  specex::image_data coeff2d(ncoeff,nfibers);
  int f=0;
  for(std::map<int,specex::Trace>::const_iterator it=psf.FiberTraces.begin(); it!=psf.FiberTraces.end(); it++) {
    const specex::Legendre1DPol *pol=0;
    if(is_x) pol = &(it->second.X_vs_W);
    else pol = &(it->second.Y_vs_W);
    for(int c=0;c<pol->coeff.size();c++)
      coeff2d(c,f)=pol->coeff(c);
    f++;
  }

  int status = 0;
  fits_movabs_hdu ( fp, hdu, NULL, &status ); harp::fits::check ( status );
  
  harp::fits::img_append < double > ( fp, coeff2d.n_rows(), coeff2d.n_cols() );
  harp::fits::img_write ( fp, coeff2d.data, false );
  if(is_x)
    harp::fits::key_write(fp,"EXTNAME","XTRACE","");
  else
    harp::fits::key_write(fp,"EXTNAME","YTRACE","");
  harp::fits::key_write(fp,"WAVEMIN",wavemin,"");
  harp::fits::key_write(fp,"WAVEMAX",wavemax,"");

}
  
void specex::write_xtrace_fits_hdu(const specex::PSF& psf, fitsfile *fp, int hdu) {
  SPECEX_INFO("write_xtrace in hdu "<< hdu);
  return _write_trace(psf,fp,hdu,true); // X 
}
void specex::write_ytrace_fits_hdu(const specex::PSF& psf, fitsfile *fp, int hdu) {
  SPECEX_INFO("write_ytrace in hdu "<< hdu);
  return _write_trace(psf,fp,hdu,false); // Y 
}

void _read_trace(specex::PSF_p psf, fitsfile *fp, int hdu, bool isx) {
  // read image of coefficients
  harp::fits::img_seek ( fp, hdu);
  size_t nrows,ncols;  
  harp::fits::img_dims ( fp, nrows, ncols );  
  specex::image_data coeff(ncols,nrows);
  harp::fits::img_read ( fp, coeff.data , false);
  int nfibers = coeff.n_rows();
  int ncoefs  = coeff.n_cols();
  double wavemin; harp::fits::key_read(fp,"WAVEMIN",wavemin);
  double wavemax; harp::fits::key_read(fp,"WAVEMAX",wavemax);
    
  for(int fiber=0;fiber<nfibers; fiber++) {
    // create or not the trace
    specex::Trace& trace = psf->FiberTraces[fiber];
    trace.fiber=fiber;
    specex::Legendre1DPol pol(ncoefs-1,wavemin,wavemax);
    for(int i=0;i<ncoefs;i++)
      pol.coeff[i]=coeff(i,fiber);
    if(isx) trace.X_vs_W=pol;
    else trace.Y_vs_W=pol;
    trace.synchronized = false; // will do synchro after
  }
}
  
void specex::read_xtrace_fits_hdu(specex::PSF_p psf, fitsfile *fp, int hdu) {
  SPECEX_INFO("read XTRACE in HDU "<< hdu);
  return _read_trace(psf,fp,hdu,true); // X 
}
void specex::read_ytrace_fits_hdu(specex::PSF_p psf, fitsfile *fp, int hdu) {
  SPECEX_INFO("read YTRACE in HDU "<< hdu);
  return _read_trace(psf,fp,hdu,false); // Y 
}

void specex::write_spots_xml(const std::vector<specex::Spot_p>& spots, const std::string& filename) {
  std::ofstream os(filename.c_str());
  boost::archive::xml_oarchive xml_oa ( os );
  xml_oa << BOOST_SERIALIZATION_NVP(spots);
  os.close();
  SPECEX_INFO("wrote spots in " << filename);
}

void specex::write_psf_fits_image(const specex::PSF_p psf, const string& filename, const int fiber, const double& wavelength, int oversampling) {
  
  double x=psf->Xccd(fiber,wavelength);
  double y=psf->Yccd(fiber,wavelength);
  SPECEX_INFO("PSF center X Y = " << x << " " << y);
  x = int(x);
  y = int(y);
  
  harp::vector_double P=psf->AllLocalParamsFW(fiber,wavelength);
  
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

  
#include <specex_psf_io_gauss_hermite_two_psf_fits_1.cpp>
#include <specex_psf_io_gauss_hermite_two_psf_fits_2.cpp>
#include <specex_psf_io_gauss_hermite_psf_fits_1.cpp>
#include <specex_psf_io_gauss_hermite_psf_fits_2.cpp>
#include <specex_psf_io_gauss_hermite_psf_fits_3.cpp>



void specex::write_psf_fits(const specex::PSF_p psf, const string& path, std::vector<specex::Spot_p> *spots) {
  fitsfile * fp;  
  harp::fits::create(fp,path);

  specex::write_xtrace_fits_hdu(*psf,fp,1);
  specex::write_ytrace_fits_hdu(*psf,fp,2);
  
  if(psf->Name()=="GaussHermitePSF")
    write_gauss_hermite_psf_fits_version_3((const specex::GaussHermitePSF&)*psf,fp,3,spots);
  else 
    SPECEX_ERROR("specex::write_psf_fits not implemented for PSF '" << psf->Name() << '"');
  
  harp::fits::close(fp);
  SPECEX_INFO("wrote PSF in " << path);
}


void specex::read_psf_fits(specex::PSF_p& psf, const string& filename, int first_bundle, int last_bundle) {
  fitsfile * fp;
  harp::fits::open_read(fp,filename);

  // first look at primary header to get psf type and I/O version
  int status=0;
  fits_movabs_hdu ( fp, 1, NULL, &status ); harp::fits::check ( status );  
  string psftype; harp::fits::key_read(fp,"PSFTYPE",psftype);
  int psfver; harp::fits::key_read (fp,"PSFVER",psfver);
  SPECEX_INFO("PSFTYPE=" << psftype);  
  SPECEX_INFO("PSFVER=" << psfver);

  // create PSF 
  if ( psftype=="GAUSS-HERMITE") psf = specex::PSF_p(new specex::GaussHermitePSF());
  else SPECEX_ERROR("read fits not implemented for psf " << psftype);
  
  if(psfver<2) { // older formats, psf is in hdu 2 and there is no x y trace hdus

    if(psftype=="GAUSS-HERMITE" && psfver==2) read_gauss_hermite_psf_fits_version_2(psf,fp,2,first_bundle,last_bundle);
    else SPECEX_ERROR("read fits not implemented for psf " << psftype << " and I/O version " << psfver);
    
  } else { // with newer format need to parse EXTNAME in HDUs to find XTRACE,YTRACE,PSF keyword
        
    read_xtrace_fits_hdu(psf,fp,find_hdu(fp,"XTRACE"));
    read_ytrace_fits_hdu(psf,fp,find_hdu(fp,"YTRACE"));

    if(psftype=="GAUSS-HERMITE" && psfver==2) read_gauss_hermite_psf_fits_version_2(psf,fp,find_hdu(fp,"PSF"),first_bundle,last_bundle);
    else SPECEX_ERROR("read fits not implemented for psf " << psftype << " and I/O version " << psfver);
    
  }
  harp::fits::close(fp);
  SPECEX_INFO("read PSF in " << filename);
}
