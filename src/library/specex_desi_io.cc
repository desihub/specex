#include <iostream>

#include <harp.hpp>

#include "specex_psf.h"
#include "specex_trace.h"
#include "specex_desi_io.h"
#include "specex_message.h"
#include "specex_fits.h"


using namespace std;


void specex::read_DESI_traceset_in_fits(
					specex::TraceSet& traceset,
					const std::string& x_vs_wave_filename, 
					int x_vs_wave_hdu_number, 
					const std::string& y_vs_wave_filename, 
					int y_vs_wave_hdu_number,
					int required_x_vs_wave_degree ,
					int required_y_vs_wave_degree
					) {
  size_t nrows,ncols;
  specex::image_data x_vs_wave_coefs;
  specex::image_data y_vs_wave_coefs;
  double x_vs_wave_wavemin=0;
  double y_vs_wave_wavemin=0;
  double x_vs_wave_wavemax=0;
  double y_vs_wave_wavemax=0;
  
  
  // reading
  
  {
    fitsfile * fp= 0;
    harp::fits::open_read (fp,x_vs_wave_filename);
    
    if(x_vs_wave_hdu_number<0) x_vs_wave_hdu_number = find_hdu(fp,"XTRACE");
    if(x_vs_wave_hdu_number<0) { // backward compatibility
      x_vs_wave_hdu_number = find_hdu(fp,"XCOEFF");
      if(x_vs_wave_hdu_number<0) SPECEX_ERROR("didn't find extension XTRACE of XCOEFF");
      SPECEX_WARNING("Using deprecated extension XCOEFF for XTRACE");
    }
    
    harp::fits::img_seek ( fp, x_vs_wave_hdu_number);
    harp::fits::img_dims ( fp, nrows, ncols );
    x_vs_wave_coefs.resize(ncols,nrows); // note my ordering in images, first is x=col, second is y=row
    harp::fits::img_read ( fp, x_vs_wave_coefs.data, false );
    harp::fits::key_read (fp,"WAVEMIN",x_vs_wave_wavemin);
    harp::fits::key_read (fp,"WAVEMAX",x_vs_wave_wavemax);
    harp::fits::close ( fp );
  }
  {
    fitsfile * fp= 0;
    harp::fits::open_read (fp,y_vs_wave_filename);
    
    if(y_vs_wave_hdu_number<0) y_vs_wave_hdu_number = find_hdu(fp,"YTRACE");
    if(y_vs_wave_hdu_number<0) { // backward compatibility
      y_vs_wave_hdu_number = find_hdu(fp,"YCOEFF");
      if(y_vs_wave_hdu_number<0) SPECEX_ERROR("didn't find extension YTRACE of YCOEFF");
      SPECEX_WARNING("Using deprecated extension YCOEFF for YTRACE");
    }
    
    harp::fits::img_seek ( fp, y_vs_wave_hdu_number);
    harp::fits::img_dims ( fp, nrows, ncols );
    y_vs_wave_coefs.resize(ncols,nrows); // note my ordering in images, first is x=col, second is y=row
    harp::fits::img_read ( fp, y_vs_wave_coefs.data, false );
    harp::fits::key_read (fp,"WAVEMIN",y_vs_wave_wavemin);
    harp::fits::key_read (fp,"WAVEMAX",y_vs_wave_wavemax);
    harp::fits::close ( fp );
  }
  
  SPECEX_INFO("read traces in '" << x_vs_wave_filename << "' HDU "<< x_vs_wave_hdu_number << " and '" << y_vs_wave_filename << "' HDU " <<  y_vs_wave_hdu_number);    
  SPECEX_INFO("X WAVE MIN MAX, SIZE " << x_vs_wave_wavemin << " " << x_vs_wave_wavemax << " " << x_vs_wave_coefs.n_cols() << "x" << x_vs_wave_coefs.n_rows());
  SPECEX_INFO("Y WAVE MIN MAX, SIZE " << y_vs_wave_wavemin << " " << y_vs_wave_wavemax << " " << y_vs_wave_coefs.n_cols() << "x" << y_vs_wave_coefs.n_rows());
  
  int nfibers=x_vs_wave_coefs.n_rows();
  if(nfibers != y_vs_wave_coefs.n_rows()) {
    SPECEX_ERROR("not same number of fibers for x and y " << x_vs_wave_coefs.n_rows() << " " << y_vs_wave_coefs.n_rows() );
  }

  traceset.resize(nfibers);
  
  for(int fiber=0;fiber<nfibers;fiber++) {
    specex::Trace& trace = traceset[fiber];
    trace.fiber = fiber;
    trace.mask  = 0; // need to do something here?
    trace.X_vs_W.xmin = x_vs_wave_wavemin;
    trace.X_vs_W.xmax = x_vs_wave_wavemax;
    if(required_x_vs_wave_degree>=0) {
      trace.X_vs_W.deg  = required_x_vs_wave_degree;
      if(required_x_vs_wave_degree<(x_vs_wave_coefs.n_cols()-1) && fiber==0)
	SPECEX_WARNING("Will use lower x trace degree that input");
    } else {
      trace.X_vs_W.deg  = x_vs_wave_coefs.n_cols()-1;
    }
    trace.X_vs_W.coeff.resize(trace.X_vs_W.deg+1);
    for(int i=0;i<min(int(trace.X_vs_W.coeff.size()),int(x_vs_wave_coefs.n_cols()));i++)
      trace.X_vs_W.coeff[i]= x_vs_wave_coefs(i,fiber);
    trace.Y_vs_W.xmin = y_vs_wave_wavemin;
    trace.Y_vs_W.xmax = y_vs_wave_wavemax;
    if(required_y_vs_wave_degree>=0) {
      trace.Y_vs_W.deg  = required_y_vs_wave_degree;
      if(required_y_vs_wave_degree<(y_vs_wave_coefs.n_cols()-1) && fiber==0)
	SPECEX_WARNING("Will use lower y trace degree that input");
    } else {
      trace.Y_vs_W.deg  = y_vs_wave_coefs.n_cols()-1;
    }
    trace.Y_vs_W.coeff.resize(trace.Y_vs_W.deg+1);
    for(int i=0;i<min(int(trace.Y_vs_W.coeff.size()),int(y_vs_wave_coefs.n_cols()));i++)
      trace.Y_vs_W.coeff[i]= y_vs_wave_coefs(i,fiber);

    trace.X_vs_Y.deg  = trace.Y_vs_W.deg + 1; // add one for inversion
    trace.W_vs_Y.deg  = trace.Y_vs_W.deg + 1; // add one for inversion
    int npts=100;
    harp::vector_double w(npts);
    harp::vector_double x(npts);
    harp::vector_double y(npts);
    for(int i=0;i<npts;i++) {
      w[i]=y_vs_wave_wavemin+((y_vs_wave_wavemax-y_vs_wave_wavemin)/(npts-1))*i;
      x[i]=trace.X_vs_W.Value(w[i]);
      y[i]=trace.Y_vs_W.Value(w[i]);      
    }
    trace.X_vs_Y.Fit(y,x);
    trace.W_vs_Y.Fit(y,w);
    trace.synchronized = true;
  }
  
  SPECEX_DEBUG("done reading traceset");
}

typedef std::map<std::string,std::string> Header;
Header read_header(const std::string& filename, int hdu) {
  fitsfile * fp= 0;
  harp::fits::open_read (fp,filename);
  harp::fits::img_seek ( fp, hdu);
  
  int status = 0;
  int nkeys;
  int morekeys;
  int ret = fits_get_hdrspace ( fp, &nkeys, &morekeys, &status );
  harp::fits::check ( status );

  // get each key in order and append to the ptree
  
  
  char keyname[200];
  char keyval[200];
  char keycom[200];
  Header head;
  for ( int i = 0; i < nkeys; ++i ) {
    ret = fits_read_keyn ( fp, i+1, keyname, keyval, keycom, &status );
    harp::fits::check ( status );
    head[keyname]=keyval;
    //SPECEX_INFO("READ IN HEADER " << keyname << " " << keyval); 
  }
  return head;
}


void specex::read_DESI_preprocessed_image(const std::string& filename, image_data &image, image_data &weight, image_data &mask, image_data &rdnoise, std::map<std::string,std::string>& header) {
  
  SPECEX_INFO("Reading DESI preprocessed image in " << filename);
  
  // first read list of headers and extensions to see if this is correct
  fitsfile * fp= 0;
  harp::fits::open_read (fp,filename);
  int nhdus = harp::fits::nhdus(fp);
  std::vector<std::string> extname;
  std::vector<int> naxis;
  std::vector<int> naxis1;
  std::vector<int> naxis2;
  for(int i=0;i<nhdus;i++) {
    int hdu=i+1; // starts at 1
    harp::fits::img_seek ( fp, hdu);
    extname.push_back("");
    naxis.push_back(0);
    naxis1.push_back(0);
    naxis2.push_back(0);
    try {
      harp::fits::key_read (fp,"EXTNAME",extname[i]);
    }catch(...) { 
      if(i!=0) SPECEX_WARNING("could not read EXTNAME in hdu " << hdu << " of " << filename);
    }
    try {
      harp::fits::key_read (fp,"NAXIS",naxis[i]);
    }catch(...) { 
      SPECEX_WARNING("could not read NAXIS in hdu " << hdu << " of " << filename);
    }
    try {
      harp::fits::key_read (fp,"NAXIS1",naxis1[i]);
    }catch(...) { 
      SPECEX_WARNING("could not read NAXIS1 in hdu " << hdu << " of " << filename);
    }
    try {
      harp::fits::key_read (fp,"NAXIS2",naxis2[i]);
    }catch(...) { 
      SPECEX_WARNING("could not read NAXIS2 in hdu " << hdu << " of " << filename);
    }
  }
  
  int flux_hdu = -1;
  int ivar_hdu = -1;
  int mask_hdu = -1;
  int rdnoise_hdu = -1;
  for(int i=0;i<nhdus;i++) {
    int hdu=i+1; // starts at 1
    if(extname[i]=="IMAGE")
      flux_hdu=hdu;
    else if(extname[i]=="IVAR") 
      ivar_hdu=hdu;
    else if(extname[i]=="MASK") 
      mask_hdu=hdu;
    else if(extname[i]=="READNOISE") 
      rdnoise_hdu=hdu;
    else if (naxis[i]==2) {
      if(flux_hdu == -1)
	flux_hdu = hdu; // the first image is supposed to be flux
      else if (ivar_hdu == -1)
	ivar_hdu = hdu; // next is supposed to be ivar
      else if (mask_hdu == -1)
	mask_hdu = hdu; // next is supposed to be mask
      else if (rdnoise_hdu == -1)
	rdnoise_hdu = hdu; // next is supposed to be rdnoise
    }
    SPECEX_DEBUG("HDU=" << hdu << " EXTNAME='" << extname[i] << "' NAXIS=" << naxis[i] << " NAXIS1=" << naxis1[i] << " NAXIS2=" << naxis2[i]);
  }
  if(flux_hdu==-1) 
    SPECEX_ERROR("Could not figure out which HDU has flux");
  if(flux_hdu>-1)  
    SPECEX_DEBUG("Will read flux in HDU " << flux_hdu);
  if(ivar_hdu>-1)  
    SPECEX_DEBUG("Will read ivar in HDU " << ivar_hdu);
  if(mask_hdu>-1)  
    SPECEX_DEBUG("Will read mask in HDU " << mask_hdu);
  if(rdnoise_hdu>-1)  
    SPECEX_DEBUG("Will read rdnoise in HDU " << rdnoise_hdu);
  
  if(flux_hdu>-1) 
    read_fits_image(filename,flux_hdu,image);
  if(ivar_hdu>-1) 
    read_fits_image(filename,ivar_hdu,weight);
  if(mask_hdu>-1) 
    read_fits_image(filename,mask_hdu,mask);
  if(rdnoise_hdu>-1) 
    read_fits_image(filename,rdnoise_hdu,rdnoise);

  // read header
  header = read_header(filename,flux_hdu);
  
  if(ivar_hdu==-1) {
    SPECEX_WARNING("Didn't find ivar in file " << filename << ", set ivar=1");
    weight.resize(image.n_cols(),image.n_rows());
    for(int j=0;j<weight.Ny();j++) {
      for(int i=0;i<weight.Nx();i++) {
	weight(i,j)=1;
      }
    }    
  }
  if(mask_hdu==-1) {
    SPECEX_WARNING("Didn't find mask in file " << filename << ", set mask=0");
    mask.resize(image.n_cols(),image.n_rows());
    for(int j=0;j<mask.Ny();j++) {
      for(int i=0;i<mask.Nx();i++) {
	mask(i,j)=0;
      }
    } 
  }else{
    SPECEX_INFO("Set ivar=0 to pixels with mask!=0");
    for(int j=0;j<weight.Ny();j++) {
      for(int i=0;i<weight.Nx();i++) {
	if(mask(i,j)!=0)
	  weight(i,j)=0;
      }
    }
  }
  
  if(rdnoise_hdu==-1) {
    SPECEX_WARNING("Didn't find rdnoise image in file " << filename << ", look for header keyword");
    rdnoise.resize(image.n_cols(),image.n_rows());
    
    double default_rdnoise_value = 0.001; // dummy default to avoid trouble in the code;
    double rdnoise_value = default_rdnoise_value;
    
    for(Header::const_iterator it=header.begin(); it != header.end(); ++it) {
      const string& key   = it->first;
      const string& value = it->second;
      if (key.find("NOISE")!=std::string::npos && ( key.find("RD")!=std::string::npos || key.find("READ")!=std::string::npos )) {
      	SPECEX_INFO("Possible read noise header key '" << key << "'=" << value);
	if (rdnoise_value==default_rdnoise_value) {
	  SPECEX_INFO("Use this one :" << value);
	  rdnoise_value=atof(value.c_str());
	}
      }
    }
    
    for(int j=0;j<rdnoise.Ny();j++) {
      for(int i=0;i<rdnoise.Nx();i++) {
	rdnoise(i,j)=rdnoise_value;
      }
    } 
  }
}
