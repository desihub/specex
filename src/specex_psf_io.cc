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
//#include <specex_hat_hermite_psf.h>

#include "specex_psf_io_gauss_hermite_two_psf_fits_1.h"
#include "specex_psf_io_gauss_hermite_two_psf_fits_2.h"
#include "specex_psf_io_gauss_hermite_psf_fits_1.h"
#include "specex_psf_io_gauss_hermite_psf_fits_2.h"
#include "specex_psf_io_gauss_hermite_psf_fits_3.h"



void specex::read_psf_gen(specex::PSF_p& psf, const std::string& filename) {
  
  if (filename.find(".xml") != std::string::npos) {
    SPECEX_INFO("read xml file " << filename);
    specex::read_psf_xml(psf,filename);
  } else if (filename.find(".fits") != std::string::npos) {
    SPECEX_INFO("read fits file " << filename);
    specex::read_psf_fits(psf,filename);
  } else {
    SPECEX_ERROR("not sure how to read this file (expect xxx.fits or xxx.xml) " << filename);
  }
}

void specex::write_psf_xml(const specex::PSF_p psf, const std::string& filename) {
  // this function is no longer functional as xml serialization is removed
  std::ofstream os(filename.c_str());
  os.close();
  SPECEX_INFO("wrote psf in " << filename);
}

void specex::read_psf_xml(specex::PSF_p& psf, const std::string& filename) {
  
  // this function is no longer functional as xml serialization is removed
  std::ifstream is(filename.c_str());
  is.close();
  SPECEX_INFO("read psf in " << filename);
  
}

void _read_trace(specex::PSF_p psf, fitsfile *fp, int hdu, bool isx, int requested_deg) {
  // read image of coefficients
  harp::fits::img_seek ( fp, hdu);
  size_t nrows,ncols;  
  harp::fits::img_dims ( fp, nrows, ncols );  
  specex::image_data coeff(ncols,nrows);
  harp::fits::img_read ( fp, coeff.data , false);
  int nfibers = coeff.n_rows();
  int ncoefs  = coeff.n_cols();
  double WAVEMIN; harp::fits::key_read(fp,"WAVEMIN",WAVEMIN);
  double WAVEMAX; harp::fits::key_read(fp,"WAVEMAX",WAVEMAX);
  
  int requested_ncoefs=ncoefs;
  if(requested_deg>0)  requested_ncoefs=requested_deg+1;
  
  for(int fiber=0;fiber<nfibers; fiber++) {
    // create or not the trace
    specex::Trace& trace = psf->FiberTraces[fiber];
    trace.fiber=fiber;
    trace.mask=0;
    trace.synchronized = false; // will do synchro after
    specex::Legendre1DPol pol(requested_ncoefs-1,WAVEMIN,WAVEMAX);
    for(int i=0;i<min(ncoefs,requested_ncoefs);i++)
      pol.coeff[i]=coeff(i,fiber);
    if(isx) trace.X_vs_W=pol;
    else trace.Y_vs_W=pol;
  }
}
  
void specex::read_xtrace_fits_hdu(specex::PSF_p psf, fitsfile *fp, int hdu, int requested_deg) {
  SPECEX_INFO("read XTRACE in HDU "<< hdu);
  return _read_trace(psf,fp,hdu,true,requested_deg); // X 
}
void specex::read_ytrace_fits_hdu(specex::PSF_p psf, fitsfile *fp, int hdu, int requested_deg) {
  SPECEX_INFO("read YTRACE in HDU "<< hdu);
  return _read_trace(psf,fp,hdu,false,requested_deg); // Y 
}

void specex::synchronize_traces(specex::PSF_p psf) {
  SPECEX_INFO("synchronizing traces");
  for(std::map<int,specex::Trace>::iterator it=psf->FiberTraces.begin(); it!=psf->FiberTraces.end(); it++) {
    specex::Trace &trace = it->second;
    // X_vs_Y
    int deg     = trace.Y_vs_W.deg;
    double wmin = trace.Y_vs_W.xmin;
    double wmax = trace.Y_vs_W.xmax;
    
    int ddeg = 1; // add one degree for inversion
    harp::vector_double ty(deg+ddeg+1);
    harp::vector_double tx(deg+ddeg+1);
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
  
void specex::write_spots_xml(const std::vector<specex::Spot_p>& spots, const std::string& filename) {

  // this function is no longer functional as xml serialization is removed
  std::ofstream os(filename.c_str());
  os.close();
  SPECEX_INFO("wrote spots in " << filename);
}

void specex::write_psf_fits_image(const specex::PSF_p psf, const string& filename, const int fiber, const double& wavelength, int oversampling) {
  
  double x=psf->Xccd(fiber,wavelength);
  double y=psf->Yccd(fiber,wavelength);
  SPECEX_INFO("PSF center X Y = " << x << " " << y);
  bool force_center=false;
  if (force_center) {
    x = int(x);
    y = int(y);
  }
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
      //if(j==0) SPECEX_INFO("x[" << i << "]=" << ib+dx << " y[" << j << "]=" << jb+dy);
      
      img(i,j)=psf->PSFValueWithParamsXY(x-dx,y-dy,ib+int(x+0.5),jb+int(y+0.5),P,0,0);
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

void specex::write_psf_fits_dummy(const std::string& path){
  fitsfile * fp;
  int status = 0;
  long naxes[2];
  double dummy;
  naxes[0] = 10;
  naxes[1] = 15;
  harp::fits::create(fp,path);  
  fits_create_img(fp,-64, 2, naxes, &status);
  harp::fits::close(fp);
}

void specex::read_traceset_fits(specex::PSF_p psf, fitsfile * fp, int degx, int degy) {
  read_xtrace_fits_hdu(psf,fp,find_hdu(fp,"XTRACE","XCOEFF"),degx);
  read_ytrace_fits_hdu(psf,fp,find_hdu(fp,"YTRACE","YCOEFF"),degy);
  synchronize_traces(psf);
}

void specex::read_traceset_fits(specex::PSF_p psf, const string& filename, int degx, int degy) {
  fitsfile * fp;
  harp::fits::open_read(fp,filename);
  read_traceset_fits(psf,fp,degx,degy);
  harp::fits::close(fp);
  SPECEX_INFO("read trace set in " << filename);
}

void specex::read_psf_fits(specex::PSF_p& psf, const string& filename) {
  fitsfile * fp;
  harp::fits::open_read(fp,filename);

  // first look at primary header to get psf type and I/O version
  string psftype; harp::fits::key_read(fp,"PSFTYPE",psftype);
  int psfver; harp::fits::key_read (fp,"PSFVER",psfver);
  SPECEX_INFO("PSFTYPE=" << psftype);  
  SPECEX_INFO("PSFVER=" << psfver);

  // create PSF 
  if ( psftype=="GAUSS-HERMITE") psf = specex::PSF_p(new specex::GaussHermitePSF());
  else SPECEX_ERROR("read fits not implemented for psf " << psftype);

  if(psfver<2) { // older formats, psf is in hdu 2 and there is no x y trace hdus

    if(psftype=="GAUSS-HERMITE" && psfver==2) read_gauss_hermite_psf_fits_version_2(psf,fp,2);
    else SPECEX_ERROR("read fits not implemented for psf " << psftype << " and I/O version " << psfver);
    
  } else { // with newer format need to parse EXTNAME in HDUs to find XTRACE,YTRACE,PSF keyword
    
    read_traceset_fits(psf,fp);    
    
    if(psftype=="GAUSS-HERMITE" && psfver==2) read_gauss_hermite_psf_fits_version_2(psf,fp,find_hdu(fp,"PSF"));
    else if(psftype=="GAUSS-HERMITE" && psfver==3) read_gauss_hermite_psf_fits_version_3(psf,fp,find_hdu(fp,"PSF"));
    else SPECEX_ERROR("read fits not implemented for psf " << psftype << " and I/O version " << psfver);
    
  }

  harp::fits::close(fp);
  SPECEX_INFO("read PSF in " << filename);
}
