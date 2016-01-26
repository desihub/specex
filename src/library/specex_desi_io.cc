#include <iostream>

#include <harp.hpp>

#include "specex_psf.h"
#include "specex_trace.h"
#include "specex_desi_io.h"
#include "specex_message.h"



using namespace std;


/*
  input format of fits file
  HDU 0 : blank
  HDU 1 ELECTRONS : 2D image of noisy electrons
  - RDNOISE is in header
  HDU 2 IVAR : Inverse variance [1/electrons^2]
  - ivar = 1/(pix.clip(0) + rdnoise**2)
  HDU 3 MASK : 0=good.  Currently all 0.
  HDU 4 XCOEFF : Legendre coefficients for mapping wavelength -> x
  - WAVEMIN, WAVEMAX : domain for mapping to [-1,1] for Legendre polynomials
  - image is coefficients for each fiber
  HDU 5 YCOEFF : Legendre coefficients for mapping wavelength -> y
  HDU 5 TRUE_ELECTRONS : original noiseless image in electrons
  ELECTRONS = poisson(TRUE_ELECTRONS) + gaussian(rdnoise)
 */

void specex::read_DESI_traceset_in_fits(
					specex::TraceSet& traceset,
					const std::string& x_vs_wave_filename, 
					int x_vs_wave_hdu_number, 
					const std::string& y_vs_wave_filename, 
					int y_vs_wave_hdu_number
					) {
  size_t nrows,ncols;
  specex::image_data x_vs_wave_coefs;
  specex::image_data y_vs_wave_coefs;
  double x_vs_wave_wavemin=0;
  double y_vs_wave_wavemin=0;
  double x_vs_wave_wavemax=0;
  double y_vs_wave_wavemax=0;
  
  
  // reading
  
  try{
    fitsfile * fp= 0;
    harp::fits::open_read (fp,x_vs_wave_filename);
    harp::fits::img_seek ( fp, x_vs_wave_hdu_number);
    harp::fits::img_dims ( fp, nrows, ncols );
    x_vs_wave_coefs.resize(ncols,nrows); // note my ordering in images, first is x=col, second is y=row
    harp::fits::img_read ( fp, x_vs_wave_coefs.data, false );
    harp::fits::key_read (fp,"WAVEMIN",x_vs_wave_wavemin);
    harp::fits::key_read (fp,"WAVEMAX",x_vs_wave_wavemax);
    harp::fits::close ( fp ); 
  }catch(...) { SPECEX_ERROR("could not read properly x traces in " << x_vs_wave_filename << " HDU " << x_vs_wave_hdu_number);}
  try{
    fitsfile * fp= 0;
    harp::fits::open_read (fp,y_vs_wave_filename);
    harp::fits::img_seek ( fp, y_vs_wave_hdu_number);
    harp::fits::img_dims ( fp, nrows, ncols );
    y_vs_wave_coefs.resize(ncols,nrows); // note my ordering in images, first is x=col, second is y=row
    harp::fits::img_read ( fp, y_vs_wave_coefs.data, false );
    harp::fits::key_read (fp,"WAVEMIN",y_vs_wave_wavemin);
    harp::fits::key_read (fp,"WAVEMAX",y_vs_wave_wavemax);
    
    harp::fits::close ( fp ); 
  }catch(...) { SPECEX_ERROR("could not read properly y traces in " << y_vs_wave_filename << " HDU " << y_vs_wave_hdu_number);}
    
  SPECEX_INFO("X WAVEMIN WAVEMAX " << x_vs_wave_wavemin << " " << x_vs_wave_wavemax);
  SPECEX_INFO("X WAVE COEF SIZE  " << x_vs_wave_coefs.n_cols() << " " << x_vs_wave_coefs.n_rows());
  SPECEX_INFO("Y WAVEMIN WAVEMAX " << y_vs_wave_wavemin << " " << y_vs_wave_wavemax);
  SPECEX_INFO("Y WAVE COEF SIZE  " << y_vs_wave_coefs.n_cols() << " " << y_vs_wave_coefs.n_rows());
  SPECEX_INFO("read traces in '" << x_vs_wave_filename << "' HDU "<< x_vs_wave_hdu_number << " and '" << y_vs_wave_filename << "' HDU " <<  y_vs_wave_hdu_number);
  
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
    trace.X_vs_W.deg  = x_vs_wave_coefs.n_cols()-1; // add 1 to minimise errors
    trace.X_vs_W.coeff.resize(trace.X_vs_W.deg+1);
    for(int i=0;i<=trace.X_vs_W.deg;i++)
      trace.X_vs_W.coeff[i]= x_vs_wave_coefs(i,fiber);
    trace.Y_vs_W.xmin = y_vs_wave_wavemin;
    trace.Y_vs_W.xmax = y_vs_wave_wavemax;
    trace.Y_vs_W.deg  = y_vs_wave_coefs.n_cols()-1; // add 1 to minimise errors
    trace.Y_vs_W.coeff.resize(trace.Y_vs_W.deg+1);
    for(int i=0;i<=trace.Y_vs_W.deg;i++)
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
  
  SPECEX_INFO("done reading traceset");
}

void specex::read_DESI_keywords(const std::string& arc_image_filename, std::map<std::string,std::string>& infos, int hdu) {
  
  fitsfile * fp= 0;
  harp::fits::open_read (fp,arc_image_filename);
  harp::fits::img_seek ( fp, hdu);
  infos.clear();
  string keys[]={"RDNOISE"};
  for(int k=0;k<1;k++) {
    harp::fits::key_read (fp,keys[k],infos[keys[k]]);
  }

  harp::fits::close(fp);
}
