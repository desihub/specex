
#include <cmath>
#include <assert.h>
#include <ctime>



#include "harp.hpp"

#include "specex_hermite.h"
#include "specex_gauss_hermite_psf.h"
#include "specex_linalg.h"
#include "specex_message.h"

using namespace std;



specex::GaussHermitePSF::GaussHermitePSF(int ideg) {
  name = "GaussHermitePSF";
  SetDegree(ideg);
  
}

#define TWO_GAUSSIANS

void specex::GaussHermitePSF::SetDegree(const int ideg) {
  SPECEX_INFO("Gauss-Hermite PSF set degree " << ideg);
  degree = ideg;
  
  paramNames.clear();
  paramNames.push_back("GHSIGX");
  paramNames.push_back("GHSIGY");

#ifdef TWO_GAUSSIANS
  paramNames.push_back("GHSIGX2");
  paramNames.push_back("GHSIGY2");
  paramNames.push_back("GHSCAL2");
#endif 
  
  char n[10];
  for(int j=0;j<degree+1;j++) {
    for(int i=0;i<degree+1;i++) {
      if(i==0 && j==0) continue;
      //if(i==1 && j==0) continue;
      //if(i==0 && j==1) continue;
      //if(i==1 && j==1) continue;
      
      sprintf(n,"GH-%d-%d",i,j);
      paramNames.push_back(n);
    }
  }
}

int specex::GaussHermitePSF::LocalNAllPar() const {
    
  int npar = 2; // sigma_x and sigma_y
  
#ifdef TWO_GAUSSIANS
  npar += 3;
#endif

  // npar += (degree+1)*(degree+1)-4;// skip (0,0)(1,0)(0,1)(1,1) : -1 because normalized, -3 because centered 
  npar += (degree+1)*(degree+1)-1;// skip (0,0) : -1 because normalized
  

  return npar;
}

double specex::GaussHermitePSF::Profile(const double &input_X, const double &input_Y,
			  const harp::vector_double &Params,
			  harp::vector_double *PosDer,
			  harp::vector_double *ParamDer) const
{

  double sigma_x_inv = 1./Params(0);
  double sigma_y_inv = 1./Params(1);
  double x = input_X*sigma_x_inv;
  double y = input_Y*sigma_y_inv;
#ifdef TWO_GAUSSIANS
  double sigma_x_2_inv = 1./Params(2);
  double sigma_y_2_inv = 1./Params(3);
  double scale_2       = Params(4);
  double scale_1       = 1-scale_2;
  double x_2 = input_X*sigma_x_2_inv;
  double y_2 = input_Y*sigma_y_2_inv;
#endif
  

 
  int nx=(degree+1);
  int ny=(degree+1);
  //int nc=nx*ny-4; // skip (0,0)(1,0)(0,1)(1,1)
  int nc=nx*ny-1; // skip (0,0)
  
  // precompute to go faster
  harp::vector_double Monomials;
  harp::vector_double Monomials_dx;
  harp::vector_double Monomials_dy;
  double prefactor=1;
#ifndef TWO_GAUSSIANS
  int first_hermite_param_index = 2; // first 2 params are sigmas
  double expfact=1./(2*M_PI)*sigma_x_inv*sigma_y_inv*exp(-0.5*(x*x+y*y));
#else
  int first_hermite_param_index = 5; // first 5 params are sigmas
  double expfact_1_part=1/(2*M_PI)*sigma_x_inv*sigma_y_inv*exp(-0.5*(x*x+y*y));
  double expfact_2_part=1/(2*M_PI)*sigma_x_2_inv*sigma_y_2_inv*exp(-0.5*(x_2*x_2+y_2*y_2));
  double expfact_1=scale_1*expfact_1_part;
  double expfact_2=scale_2*expfact_2_part;
  double expfact=expfact_1+expfact_2;
#endif


  if(PosDer==0 && ParamDer==0) {
    int param_index=first_hermite_param_index;
    for(int j=0;j<ny;j++) {
      double Hyj=HermitePol(j,y);
      //int imin=0; if(j<2) {imin=2;} // skip (0,0)(1,0)(0,1)(1,1)
      int imin=0; if(j==0) imin=1; // skip (0,0)
      for(int i=imin;i<nx;i++,param_index++) {
	prefactor+=Params(param_index)*Hyj*HermitePol(i,x);
      }
    }
    return expfact*prefactor;
    
  } else if(PosDer==0 && ParamDer) {
    Monomials.resize(nc);
    int index=0;
    for(int j=0;j<ny;j++) {
      double Hyj=HermitePol(j,y);
      //int imin=0; if(j<2) {imin=2;} // skip (0,0)(1,0)(0,1)(1,1)
      int imin=0; if(j==0) imin=1; // skip (0,0)
      for(int i=imin;i<nx;i++,index++) {
	Monomials[index]=Hyj*HermitePol(i,x);
      }
    }
    prefactor += specex::dot(Monomials,ublas::project(Params,ublas::range(first_hermite_param_index,first_hermite_param_index+nc)));
  
  } else if(PosDer && ParamDer==0) {
    Monomials_dx.resize(nc);
    Monomials_dy.resize(nc);
    int param_index=first_hermite_param_index;
    int index=0;
    for(int j=0;j<ny;j++) {
      double Hyj=HermitePol(j,y);
      double dHyj=HermitePolDerivative(j,y);
      
      //int imin=0; if(j<2) {imin=2;} // skip (0,0)(1,0)(0,1)(1,1)
      int imin=0; if(j==0) imin=1; // skip (0,0)
      for(int i=imin;i<nx;i++,index++,param_index++) {
	
	double Hxi=HermitePol(i,x);
	prefactor+=Params(param_index)*Hxi*Hyj;
	Monomials_dx[index]=Hyj*HermitePolDerivative(i,x);
	Monomials_dy[index]=dHyj*Hxi;
      }
    }
  } else if(PosDer && ParamDer) {
    Monomials.resize(nc);
    Monomials_dx.resize(nc);
    Monomials_dy.resize(nc);
    int index=0;
    for(int j=0;j<ny;j++) {
      double Hyj=HermitePol(j,y);
      double dHyj=HermitePolDerivative(j,y);
      
      //int imin=0; if(j<2) {imin=2;} // skip (0,0)(1,0)(0,1)(1,1)
      int imin=0; if(j==0) imin=1; // skip (0,0)
      for(int i=imin;i<nx;i++,index++) {
	double Hxi=HermitePol(i,x);
	Monomials[index]=Hyj*Hxi;
	Monomials_dx[index]=Hyj*HermitePolDerivative(i,x);
	Monomials_dy[index]=dHyj*Hxi;
      }
    }
    prefactor += specex::dot(Monomials,ublas::project(Params,ublas::range(first_hermite_param_index,first_hermite_param_index+nc)));
  }
  
  
  double psf_val = expfact*prefactor;
  
  if(ParamDer) {
#ifndef TWO_GAUSSIANS
    (*ParamDer)[0] = (x*x-1)*sigma_x_inv*psf_val; // exact ONLY if all gauss-hermite terms except zeroth order  = 0
    (*ParamDer)[1] = (y*y-1)*sigma_y_inv*psf_val; // exact ONLY if all gauss-hermite terms except zeroth order  = 0
#else
    double psf_val_1 = expfact_1*prefactor;
    double psf_val_2 = expfact_2*prefactor;
    
    (*ParamDer)[0] = (x*x-1)*sigma_x_inv*psf_val_1;
    (*ParamDer)[1] = (y*y-1)*sigma_y_inv*psf_val_1;
    (*ParamDer)[2] = (x_2*x_2-1)*sigma_x_2_inv*psf_val_2;
    (*ParamDer)[3] = (y_2*y_2-1)*sigma_y_2_inv*psf_val_2;
    (*ParamDer)[4] = (-expfact_1_part+expfact_2_part)*prefactor;
#endif
    ublas::project(*ParamDer,ublas::range(first_hermite_param_index,first_hermite_param_index+nc)) = expfact*Monomials;
  }
  
  if(PosDer) { // wrong sign on purpose (derivatives w.r.t -X)
    double& dvdx=(*PosDer)(0);
    double& dvdy=(*PosDer)(1);  
   
#ifndef TWO_GAUSSIANS 
    dvdx=x*sigma_x_inv*psf_val;
    dvdy=y*sigma_y_inv*psf_val;
    
    double d_poly_dx = specex::dot(Monomials_dx,ublas::project(Params,ublas::range(first_hermite_param_index,first_hermite_param_index+nc)));
    double d_poly_dy = specex::dot(Monomials_dy,ublas::project(Params,ublas::range(first_hermite_param_index,first_hermite_param_index+nc)));
    dvdx -= d_poly_dx*expfact*sigma_x_inv; // minus sign cause derivative wrt -x
    dvdy -= d_poly_dy*expfact*sigma_y_inv;
#else 
    
    double psf_val_1 = expfact_1*prefactor;
    double psf_val_2 = expfact_2*prefactor;
    
    dvdx=x*sigma_x_inv*psf_val_1;
    dvdy=y*sigma_y_inv*psf_val_1;
    dvdx += x_2*sigma_x_2_inv*psf_val_2;
    dvdy += y_2*sigma_y_2_inv*psf_val_2;
    
    double d_poly_dx = specex::dot(Monomials_dx,ublas::project(Params,ublas::range(first_hermite_param_index,first_hermite_param_index+nc)));
    double d_poly_dy = specex::dot(Monomials_dy,ublas::project(Params,ublas::range(first_hermite_param_index,first_hermite_param_index+nc)));
    dvdx -= d_poly_dx*expfact*sigma_x_inv; // minus sign cause derivative wrt -x
    dvdy -= d_poly_dy*expfact*sigma_y_inv;



#endif
}
  
  return psf_val;
}

  
harp::vector_double specex::GaussHermitePSF::DefaultParams() const 
{
  
  harp::vector_double Params(LocalNAllPar());
  Params.clear(); // all = zero at beginning = a pure gaussian
  Params(0) = 1.1; // this is sigma_x
  Params(1) = 1.1; // this is sigma_y
#ifdef TWO_GAUSSIANS
  Params(2) = 3.; // this is sigma_x_2
  Params(3) = 3.; // this is sigma_y_2
  Params(4) = 0.1; // this is the amplitude of the second gaussian
#endif
  return Params;
}



/////////////////////////////////////////////
//////////////  Fits format IO //////////////
/////////////////////////////////////////////




#include <specex_image_data.h>

void specex::GaussHermitePSF::WriteFits(fitsfile* fp, int first_hdu) const {
  WriteFits_v0(fp,first_hdu);
}

void specex::GaussHermitePSF::ReadFits(fitsfile* fp, int first_hdu) {
  
  int status = 0;  
  fits_movabs_hdu ( fp, first_hdu, NULL, &status ); harp::fits::check ( status );
  
  {
    string psf_type;
    harp::fits::key_read (fp,"PSFTYPE", psf_type);
    if(psf_type != "GAUSS-HERMITE") SPECEX_ERROR("PSF in file is not 'GAUSS-HERMITE' but '" << psf_type << "'");
  }
  
  string psf_version;
  harp::fits::key_read (fp,"PSFVER", psf_version);
  
  if(psf_version=="0")
    ReadFits_v0(fp,first_hdu);
  else 
    SPECEX_ERROR("Cannot read GAUSS-HERMITE PSF format version '" << psf_version << "' only version '0' is implemented so far");
  
}



/////////////////////////////////////////////
///////////// IO VERSION 0 //////////////////
/////////////////////////////////////////////


void specex::GaussHermitePSF::WriteFits_v0(fitsfile* fp, int first_hdu) const {
  
  ////////////////////////////
  string PSFVER = "0";
  ////////////////////////////
  
  int NPIX_X = ccd_image_n_cols;
  int NPIX_Y = ccd_image_n_rows;
  
  int bundle_min = 1000;
  int bundle_max = 0;
  for(std::map<int,PSF_Params>::const_iterator bundle_it = ParamsOfBundles.begin();
      bundle_it != ParamsOfBundles.end(); ++bundle_it) {
    if(bundle_it->first > bundle_max) bundle_max = bundle_it->first;
    if(bundle_it->first < bundle_min) bundle_min = bundle_it->first;
  }
  
  int NFLUX = NPIX_Y; // this has to be a parameter
  int NSPEC = FiberTraces.size(); // this also has to be a parameter
  double ystep = ccd_image_n_rows/NFLUX; // this has to be a parameter
  
  
  SPECEX_INFO("NFLUX=" << NFLUX << " NSPEC=" << NSPEC);
	      
  specex::image_data x_data(NFLUX,NSPEC);
  specex::image_data y_data(NFLUX,NSPEC);
  specex::image_data wave_data(NFLUX,NSPEC);
  
  int ispec=0;
  for(std::map<int,Trace>::const_iterator fiber_it = FiberTraces.begin(); fiber_it != FiberTraces.end(); ++fiber_it, ++ispec) {
    
    const specex::Trace& trace = fiber_it->second;
    
    for(int i=0;i<NFLUX;i++) {
      double y = ystep*i;

      y_data(i,ispec)=y;
      x_data(i,ispec)=trace.X_vs_Y.Value(y);
      wave_data(i,ispec)=trace.W_vs_Y.Value(y);
      //cout << i << " " << x_data(i,ispec) << " " << y << " " << wave_data(i,ispec) << endl;
    }
  }
  
  string default_comment = "";

  
  int status = 0;
  
  fits_movabs_hdu ( fp, first_hdu, NULL, &status ); harp::fits::check ( status );
  
  harp::fits::img_append < double > ( fp, NSPEC, NFLUX );
  harp::fits::img_write (fp, x_data.data);
  
  fits_write_comment(fp,"PSF generated by specex, https://github.com/julienguy/specex",&status); harp::fits::check ( status );
  
  {
    char date_comment[80];
    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );
    sprintf(date_comment,"PSF fit date %04d-%02d-%02d",(now->tm_year + 1900),(now->tm_mon + 1),now->tm_mday);
    fits_write_comment(fp,date_comment,&status); harp::fits::check ( status );
  }


  harp::fits::key_write(fp,"MJD",(long long int)mjd,"MJD of arc lamp exposure");
  harp::fits::key_write(fp,"PLATEID",plate_id,"plate ID of arc lamp exposure");
  harp::fits::key_write(fp,"CAMERA",camera_id,"camera ID");
  
  harp::fits::key_write(fp,"PSFTYPE","GAUSS-HERMITE",default_comment);
  harp::fits::key_write(fp,"PSFVER",PSFVER,default_comment);
  harp::fits::key_write(fp,"ARCEXP",arc_exposure_id,"ID of arc lamp exposure used to fit PSF");
  
  harp::fits::key_write(fp,"NPIX_X",(long long int)NPIX_X,"number of columns in input CCD image");
  harp::fits::key_write(fp,"NPIX_Y",(long long int)NPIX_Y,"number of rows in input CCD image");
  harp::fits::key_write(fp,"NFLUX",(long long int)NFLUX,"number of flux elements");
  harp::fits::key_write(fp,"NSPEC",(long long int)NSPEC,"number of fibers/spectra");
  harp::fits::key_write(fp,"BUNDLMIN",(long long int)bundle_min,"first bundle of fibers");
  harp::fits::key_write(fp,"BUNDLMAX",(long long int)bundle_max,"last bundle of fibers");
  
  harp::fits::key_write(fp,"PSFPARAM","X",default_comment);
  
     
  harp::fits::img_append < double > ( fp, NSPEC, NFLUX );
  harp::fits::img_write (fp, y_data.data);
  harp::fits::key_write(fp,"PSFPARAM","Y",default_comment);
  harp::fits::key_write(fp,"EXTNAME","Y",default_comment);
  
  
  harp::fits::img_append < double > ( fp, NSPEC, NFLUX );
  harp::fits::img_write (fp, wave_data.data);
  harp::fits::key_write(fp,"PSFPARAM","WAVELENGTH",default_comment);
  harp::fits::key_write(fp,"EXTNAME","WAVELENGTH",default_comment);
  

  
  // now create image of coefficients for each bundle of fibers
  // ============================================================================================== 
  
  for(std::map<int,PSF_Params>::const_iterator bundle_it = ParamsOfBundles.begin();
      bundle_it != ParamsOfBundles.end(); ++bundle_it) {
    
    const PSF_Params& params_of_bundle = bundle_it->second;
    
    int BUNDLEID = bundle_it->first;
    int FIBERMIN = params_of_bundle.fiber_min;
    int FIBERMAX = params_of_bundle.fiber_max;
    int LDEGX = params_of_bundle.AllParPolXW[0]->xdeg;
    int LDEGY = params_of_bundle.AllParPolXW[0]->ydeg;
    double LXMIN = params_of_bundle.AllParPolXW[0]->xmin;
    double LXMAX = params_of_bundle.AllParPolXW[0]->xmax;
    double LYMIN = params_of_bundle.AllParPolXW[0]->ymin;
    double LYMAX = params_of_bundle.AllParPolXW[0]->ymax;
    
    // check all params_of_bundle have same degree !!
    for(size_t i=1;i<params_of_bundle.AllParPolXW.size();i++) {
      if(LDEGX != params_of_bundle.AllParPolXW[i]->xdeg) SPECEX_ERROR("Need same degree for all Legendre2DPol for complying with fits format");
      if(LDEGY != params_of_bundle.AllParPolXW[i]->ydeg) SPECEX_ERROR("Need same degree for all Legendre2DPol for complying with fits format");
      if(LXMIN != params_of_bundle.AllParPolXW[i]->xmin) SPECEX_ERROR("Need same xy range for all Legendre2DPol for complying with fits format");
      if(LXMAX != params_of_bundle.AllParPolXW[i]->xmax) SPECEX_ERROR("Need same xy range for all Legendre2DPol for complying with fits format");
      if(LYMIN != params_of_bundle.AllParPolXW[i]->ymin) SPECEX_ERROR("Need same xy range for all Legendre2DPol for complying with fits format");
      if(LYMAX != params_of_bundle.AllParPolXW[i]->ymax) SPECEX_ERROR("Need same xy range for all Legendre2DPol for complying with fits format"); 
    }
    
    int GHDEGX  = degree; 
    int GHDEGY  = degree;
   
    size_t ncoefs_gauss_hermite = (GHDEGX+1)*(GHDEGY+1);
    size_t ncoefs_legendre = (LDEGX+1)*(LDEGY+1);
    
    long ncoefs  = ncoefs_gauss_hermite*ncoefs_legendre;
    
  
    // GHX = AXIS1 ?
    // GHY = AXIS2 ?
    // LX  = AXIS3 ?
    // LY  = AXIS4 ?
    // index = A1 + NAXIS1*(A2 + NAXIS2*(A3 + NAXIS3*(A4))) ?
    
    
    double* buffer = new double[ncoefs];
    for(size_t i=0;i<ncoefs;i++) buffer[i]=0; // set to 0
    
    // loop on gauss-hermite parameters
    buffer[0]=1; // order 0 legendre pol for GH00 is = 1
    
    
    // WILL BE first 4 gauss hermite coefficients are by default GH00=1, GH01=0, GH10=0 , GH11=0    
    // IS first 3 gauss hermite coefficients are by default GH00=1, GH01=0, GH10=0  
    
    int gh_index=0;
    int buffer_index=0;
    for(int j_gh=0;j_gh<=GHDEGY;j_gh++) { // loop on j=y=rows first
      for(int i_gh=0;i_gh<=GHDEGX;i_gh++) {
	// skip some coeffs
	if(i_gh == 0 && j_gh == 0) {buffer_index += ncoefs_legendre; continue;}
	if(i_gh == 1 && j_gh == 0) {buffer_index += ncoefs_legendre; continue;}
	if(i_gh == 0 && j_gh == 1) {buffer_index += ncoefs_legendre; continue;}
	
	const harp::vector_double& legendre_coefficents = params_of_bundle.AllParPolXW[ gh_index ]->coeff;
	gh_index++;

	// now loop on legendre coefficients m(i+j*(xdeg+1))
	for(size_t i=0;i<legendre_coefficents.size();i++,buffer_index++) {
	  buffer[buffer_index]=legendre_coefficents(i);
	  cout << buffer_index << " " << buffer[buffer_index] << endl;
	    
	}
      }
    }

    
    
    // and now create a 4 dimension fits image
    long naxes[4];
    naxes[0] = LDEGX+1;
    naxes[1] = LDEGY+1;
    naxes[2] = GHDEGX+1;
    naxes[3] = GHDEGY+1;
    
    long fpixel[4];
    fpixel[0] = 1;
    fpixel[1] = 1;
    fpixel[2] = 1;
    fpixel[3] = 1;
    
    fits_create_img ( fp, harp::fits::ftype< double >::bitpix(), 4, naxes, &status ); harp::fits::check ( status );
    fits_write_pix ( fp, harp::fits::ftype< double >::datatype(), fpixel, ncoefs, buffer , &status ); harp::fits::check ( status );
    harp::fits::key_write(fp,"FIBERMIN",(long long int)FIBERMIN,"First fiber for which this PSF is valid");
    harp::fits::key_write(fp,"FIBERMAX",(long long int)FIBERMAX,"Last fiber for which this PSF is valid");
    
    harp::fits::key_write(fp,"HSIZEX",(long long int)hSizeX,"Half size of PSF in fit, NX=2*HSIZEX+1");
    harp::fits::key_write(fp,"HSIZEY",(long long int)hSizeY,"Half size of PSF in fit, NY=2*HSIZEY+1");
    
    
    harp::fits::key_write(fp,"GHDEGX",(long long int)GHDEGX,"Gauss-Hermite pol. degree along x");
    harp::fits::key_write(fp,"GHDEGY",(long long int)GHDEGY,"Gauss-Hermite pol. degree along y");
    harp::fits::key_write(fp,"LDEGX",(long long int)LDEGX,"Legendre pol. degree along x");
    harp::fits::key_write(fp,"LDEGY",(long long int)LDEGY,"Legendre pol. degree along y");
    harp::fits::key_write(fp,"LXMIN",(long long int)LXMIN,"Legendre pol. xmin (for x=xmin red.x=-1)");
    harp::fits::key_write(fp,"LXMAX",(long long int)LXMAX,"Legendre pol. xmax (for x=xmax red.x=+1)");
    harp::fits::key_write(fp,"LYMIN",(long long int)LYMIN,"Legendre pol. ymin (for y=ymin red.y=-1)");
    harp::fits::key_write(fp,"LYMAX",(long long int)LYMAX,"Legendre pol. ymax (for y=ymax red.y=+1)");
    
    {
      char extname[80];
    
      sprintf(extname,"BUNDLE%02d",BUNDLEID);
      harp::fits::key_write(fp,"EXTNAME",extname,"Fiber bundle number in CCD");
    }
    
    delete [] buffer;
    
  } // end of loop on fiber bundles
}

void specex::GaussHermitePSF::ReadFits_v0(fitsfile* fp, int first_hdu) {

  SPECEX_INFO("starting reading GaussHermitePSF v0"); 
  SPECEX_WARNING("this routine specex::GaussHermitePSF::ReadFits_v0 has a bug somewhere"); 
    
  harp::fits::key_read(fp,"MJD",mjd);
  harp::fits::key_read(fp,"PLATEID",plate_id);
  harp::fits::key_read(fp,"CAMERA",camera_id); 
  harp::fits::key_read(fp,"ARCEXP",arc_exposure_id);
  
  long long int NPIX_X;
  long long int NPIX_Y;
  long long int NFLUX;
  long long int NSPEC;
  long long int bundle_min;
  long long int bundle_max;
  
  harp::fits::key_read(fp,"NPIX_X",NPIX_X);
  harp::fits::key_read(fp,"NPIX_Y",NPIX_Y);
  harp::fits::key_read(fp,"NFLUX",NFLUX);
  harp::fits::key_read(fp,"NSPEC",NSPEC);
  harp::fits::key_read(fp,"BUNDLMIN",bundle_min);
  harp::fits::key_read(fp,"BUNDLMAX",bundle_max);

  // check param
  {
    string param_name;
    harp::fits::key_read(fp,"PSFPARAM",param_name);
    if(param_name != "X") SPECEX_ERROR("Expect PSFPARAM='X' in first HDU");
  }
  
  // read X image here
  size_t nrows,ncols; 
  harp::fits::img_dims ( fp, nrows, ncols );
  if(ncols != NFLUX || nrows != NSPEC) SPECEX_ERROR("X image does not have correct dimension NFLUXxNSPEC");
  specex::image_data img_x(ncols,nrows);
  harp::fits::img_read ( fp, img_x.data );
  
  
  // change HDU
  {
    int status = 0;  
    fits_movrel_hdu ( fp, 1, NULL, &status ); harp::fits::check ( status );
  }
  // check param
  {
    string param_name;
    harp::fits::key_read(fp,"PSFPARAM",param_name);
    if(param_name != "Y") SPECEX_ERROR("Expect PSFPARAM='Y' in this HDU");
  }

  // read Y image here
  harp::fits::img_dims ( fp, nrows, ncols );
  if(ncols != NFLUX || nrows != NSPEC) SPECEX_ERROR("Y image does not have correct dimension NFLUXxNSPEC"); 
  specex::image_data img_y(ncols,nrows);
  harp::fits::img_read ( fp, img_y.data );
  
  
  // change HDU
  {
    int status = 0;
    fits_movrel_hdu ( fp, 1, NULL, &status ); harp::fits::check ( status );
  }
  // check param
  {
    string param_name;
    harp::fits::key_read(fp,"PSFPARAM",param_name);
    if(param_name != "WAVELENGTH") SPECEX_ERROR("Expect PSFPARAM='WAVELENGTH' in this HDU");
  }
  
  // read WAVELENGTH image here
  harp::fits::img_dims ( fp, nrows, ncols );
  if(ncols != NFLUX || nrows != NSPEC) SPECEX_ERROR("WAVELENGTH image does not have correct dimension NFLUXxNSPEC"); 
  specex::image_data img_wave(ncols,nrows);
  harp::fits::img_read ( fp, img_wave.data );
  

  // clear fiber traces, and psf params
  FiberTraces.clear();
  ParamsOfBundles.clear();
  
  int first_fiber_in_file = -1;

  // loop on bundles
  for(int bundle_id = bundle_min; bundle_id <= bundle_max; bundle_id ++) {
    
    SPECEX_INFO("Reading bundle " << bundle_id);

    // change HDU
    {
      int status = 0;
      fits_movrel_hdu ( fp, 1, NULL, &status ); harp::fits::check ( status );
    }
    
    {
      int naxis; 
      int status;
      fits_get_img_dim ( fp, &naxis, &status ); 
      cout << "fits status = " << status << endl; 
      harp::fits::check ( status ); 
      cout << "naxis = " << naxis << endl;
    }
    
    long long int FIBERMIN;
    long long int FIBERMAX; 
    long long int HSIZEX;
    long long int HSIZEY;
    long long int GHDEGX;
    long long int GHDEGY;
    long long int LDEGX;
    long long int LDEGY;
    long long int LXMIN;
    long long int LXMAX;
    long long int LYMIN;
    long long int LYMAX;
    //string extname;

    harp::fits::key_read(fp,"FIBERMIN",FIBERMIN);
    harp::fits::key_read(fp,"FIBERMAX",FIBERMAX);   
    harp::fits::key_read(fp,"HSIZEX",HSIZEX);
    harp::fits::key_read(fp,"HSIZEY",HSIZEY);
    harp::fits::key_read(fp,"GHDEGX",GHDEGX);
    harp::fits::key_read(fp,"GHDEGY",GHDEGY);
    harp::fits::key_read(fp,"LDEGX",LDEGX);
    harp::fits::key_read(fp,"LDEGY",LDEGY);
    harp::fits::key_read(fp,"LXMIN",LXMIN);
    harp::fits::key_read(fp,"LXMAX",LXMAX);
    harp::fits::key_read(fp,"LYMIN",LYMIN);
    harp::fits::key_read(fp,"LYMAX",LYMAX);    
    //harp::fits::key_read(fp,"EXTNAME",extname);
    
    

    if(bundle_id == bundle_min ) {
      // set global PSF properties
      hSizeX = HSIZEX;
      hSizeY = HSIZEY;
      

      first_fiber_in_file = FIBERMIN;    
      
    }else{ 
      // check they are the same or die
      if(HSIZEX != hSizeX) SPECEX_ERROR("currenlty same PSF size for all bundles described in file");
      if(HSIZEY != hSizeY) SPECEX_ERROR("currenlty same PSF size for all bundles described in file");
      if(GHDEGX != degree) SPECEX_ERROR("currenlty same PSF Gauss-Hermite degree for all bundles described in file");
      if(GHDEGY != degree) SPECEX_ERROR("currenlty same PSF Gauss-Hermite degree for all bundles described in file");
     
    }
    
    // and now read a 4 dimension fits image
    /*
    long naxes[4];
    naxes[0] = LDEGX+1;
    naxes[1] = LDEGY+1;
    naxes[2] = GHDEGX+1;
    naxes[3] = GHDEGY+1;
    */
    long fpixel[4];
    fpixel[0] = 1;
    fpixel[1] = 1;
    fpixel[2] = 1;
    fpixel[3] = 1;
    
    
    size_t ncoefs_gauss_hermite = (GHDEGX+1)*(GHDEGY+1);
    size_t ncoefs_legendre = (LDEGX+1)*(LDEGY+1);
    long ncoefs  = ncoefs_gauss_hermite*ncoefs_legendre;
    
    //fits_create_img ( fp, harp::fits::ftype< double >::bitpix(), 4, naxes, &status ); harp::fits::check ( status );
    //fits_write_pix ( fp, harp::fits::ftype< double >::datatype(), fpixel, ncoefs, buffer , &status ); harp::fits::check ( status );
    // int ffgpv(fitsfile *fptr, int  datatype, LONGLONG firstelem, LONGLONG nelem,
    //   void *nulval, void *array, int *anynul, int  *status);
    double nullval = 0;
    int anynull;
    int status;
    double* buffer = new double[ncoefs];
    for(size_t i=0;i<ncoefs;i++) buffer[i]=0; // set to 0
    
    cout << "here ncoefs = " << ncoefs << endl;
    //if(naxis !=4 ) SPECEX_ERROR("Expect a 4D image and get naxis=" << naxis);
    
    
    
    
     //fits_read_img(fp, harp::fits::ftype< double >::datatype(), 1, ncoefs, &nullval, buffer, &anynull, &status);
    fits_read_pix(fp, harp::fits::ftype< double >::datatype(), fpixel, ncoefs, &nullval, buffer, &anynull, &status);
    harp::fits::check ( status );

    // fill psf here
    
    // allocate psf parameters
    ParamsOfBundles[bundle_id] = PSF_Params();
    PSF_Params& params_of_bundle = ParamsOfBundles[bundle_id];
    params_of_bundle.bundle_id = bundle_id;
    params_of_bundle.fiber_min = FIBERMIN;
    params_of_bundle.fiber_max = FIBERMAX;
    
    //  load traces
    for(int fiber_id = params_of_bundle.fiber_min ; fiber_id <= params_of_bundle.fiber_max ; fiber_id ++ ) {
      
      SPECEX_INFO("loading trace for fiber " << fiber_id);

      FiberTraces[fiber_id] = specex::Trace();
      specex::Trace& trace = FiberTraces[fiber_id];
      trace.fiber=fiber_id;
      #warning NEED TO IMPLEMENT JUMPS
      trace.yjumplo = -1; 
      trace.yjumphi = -1;
      trace.yjumpval = -1;
      
      
      // load x,y,lw data
      harp::vector_double x(NFLUX);
      harp::vector_double y(NFLUX);
      harp::vector_double w(NFLUX);
      for(int i=0;i<NFLUX;i++) {
	x(i) = img_x(i,fiber_id-first_fiber_in_file);
	y(i) = img_y(i,fiber_id-first_fiber_in_file);
	w(i) = img_wave(i,fiber_id-first_fiber_in_file);
      }
      
      // define W vs Y legendre polynomial
      trace.X_vs_W.deg  = SPECEX_TRACE_DEFAULT_LEGENDRE_POL_DEGREE; trace.X_vs_W.coeff.resize(trace.X_vs_W.deg+1);
      trace.Y_vs_W.deg  = SPECEX_TRACE_DEFAULT_LEGENDRE_POL_DEGREE; trace.Y_vs_W.coeff.resize(trace.Y_vs_W.deg+1);
      trace.W_vs_Y.deg  = SPECEX_TRACE_DEFAULT_LEGENDRE_POL_DEGREE; trace.W_vs_Y.coeff.resize(trace.W_vs_Y.deg+1);
      trace.X_vs_Y.deg   = SPECEX_TRACE_DEFAULT_LEGENDRE_POL_DEGREE; trace.X_vs_Y.coeff.resize(trace.X_vs_Y.deg+1);
      
      // fit traces using these data
      trace.X_vs_W.Fit(w,x,NULL,true); // auto sets range
      trace.Y_vs_W.Fit(w,y,NULL,true); // auto sets range
      trace.W_vs_Y.Fit(y,w,NULL,true); // auto sets range
      trace.X_vs_Y.Fit(y,x,NULL,true); // auto sets range
      
      // ok
    }

    // load legendre psf params from buffer
    int buffer_index=0;
    for(int j_gh=0;j_gh<=GHDEGY;j_gh++) { // loop on j=y=rows first
      for(int i_gh=0;i_gh<=GHDEGX;i_gh++) {
	// skip some coeffs
	if(i_gh == 0 && j_gh == 0) {buffer_index += ncoefs_legendre; continue;}
	if(i_gh == 1 && j_gh == 0) {buffer_index += ncoefs_legendre; continue;}
	if(i_gh == 0 && j_gh == 1) {buffer_index += ncoefs_legendre; continue;}
	
	Legendre2DPol_p pol(new Legendre2DPol());
	pol->xdeg = LDEGX;
	pol->ydeg = LDEGY;
	pol->xmin = LXMIN;
	pol->xmax = LXMAX;
	pol->ymin = LYMIN;
	pol->ymax = LYMAX;
	pol->coeff.resize((pol->xdeg+1)*(pol->ydeg+1));
	
	// now loop on legendre coefficients m(i+j*(xdeg+1))
	for(size_t i=0;i<pol->coeff.size();i++,buffer_index++) {
	  pol->coeff(i) = buffer[buffer_index];
	}
	params_of_bundle.AllParPolXW.push_back(pol);
      }
    }
    
    
    delete [] buffer;
    
  
    SPECEX_INFO("read successfully fiber bundle " << bundle_id);
  }
  
}
 
