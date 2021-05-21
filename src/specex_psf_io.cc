#include <fstream>
#include <boost/algorithm/string.hpp>

#include <unhrp.h>

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
