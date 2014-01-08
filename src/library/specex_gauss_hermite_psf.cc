
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
  sigma = 1;
  need_to_resize_buffer = true;
  SetDegree(ideg);

}

void specex::GaussHermitePSF::ResizeBuffer() {
  if(!need_to_resize_buffer) return;
  Hx.resize(degree+1);
  Hy.resize(degree+1);
  dHx.resize(degree+1);
  dHy.resize(degree+1);
  need_to_resize_buffer = false;
}

void specex::GaussHermitePSF::SetDegree(const int ideg) {
  SPECEX_INFO("Gauss-Hermite PSF set degree " << ideg);
  degree = ideg;
  
  need_to_resize_buffer = true;
    
  paramNames.clear();
  char n[10];
  for(int j=0;j<degree+1;j++) {
    for(int i=0;i<degree+1;i++) {
      if(i==0 && j==0) continue;
      if(i==1 && j==0) continue;
      if(i==0 && j==1) continue;
      sprintf(n,"P%d.%d",i,j);
      paramNames.push_back(n);
    }
  }
  
  tail_norm_index=-1;
#ifndef LORENTZIAN_TAILS
  tail_power_index=-1;
  tail_x_scale_plus_index=-1;
  tail_x_scale_minus_index=-1;
  tail_y_scale_minus_index=-1;
#endif
  y_tail_norm_index=-1;

#ifdef ADD_Y_TAILS_TO_GAUSS_HERMITE   
  
  paramNames.push_back("YTailNorm");  
  y_tail_norm_index = ParamIndex("YTailNorm");
    
#endif
  
#ifdef ADD_2D_TAILS_TO_GAUSS_HERMITE   
  
  paramNames.push_back("TailNorm"); tail_norm_index = ParamIndex("TailNorm");

#ifndef LORENTZIAN_TAILS
    paramNames.push_back("TailPower"); tail_power_index = ParamIndex("TailPower");
    paramNames.push_back("TailXposScale"); tail_x_scale_plus_index  = ParamIndex("TailXposScale");
    paramNames.push_back("TailXnegScale"); tail_x_scale_minus_index = ParamIndex("TailXnegScale");
    paramNames.push_back("TailYnegScale"); tail_y_scale_minus_index = ParamIndex("TailYnegScale");
#endif    
    
#endif
    
    /*
      for(int i=0;i<int(paramNames.size());i++) cout << i << " " <<  paramNames[i] << endl;
      cout << "Npar = " << NPar() << endl;
      exit(12);
    */
  }
size_t specex::GaussHermitePSF::NPar() const {
    
    int npar = (degree+1)*(degree+1)-3;// -1 because normalized, -2 because centered 
    
#ifdef ADD_Y_TAILS_TO_GAUSS_HERMITE
    npar += 1; // normalization
#endif
    
#ifdef ADD_2D_TAILS_TO_GAUSS_HERMITE
    npar += 1; //
#ifndef LORENTZIAN_TAILS
    npar += 4; //
#endif
#endif
    
  return npar;
  }

double specex::GaussHermitePSF::Profile(const double &input_X, const double &input_Y,
			  const harp::vector_double &Params,
			  harp::vector_double *PosDer,
			  harp::vector_double *ParamDer) const
{

  if(need_to_resize_buffer) const_cast<specex::GaussHermitePSF*>(this)->ResizeBuffer();

  // direct pointers to go faster
  //const double* params = Params.Data();
  //double* posder = 0;
  //if(PosDer) {posder = PosDer->NonConstData(); PosDer->Zero();}
  //double* paramder = 0;
  //if(ParamDer) {paramder = ParamDer->NonConstData(); ParamDer->Zero();}

  assert(NPar()<=Params.size());

  double x = input_X/sigma;
  double y = input_Y/sigma;
 
  double psf_val = 0;
  
  // if(x*x+y*y<25) { // sharp cut at 5 sigma BAD IDEA BECAUSE DERIVATIVES !!
  {
    { 
      harp::vector_double& vHx = const_cast<harp::vector_double&>(Hx); 
      for(int i=0;i<=degree;i++) {
	vHx[i]=HermitePol(i,x);
    }
      harp::vector_double& vHy = const_cast<harp::vector_double&>(Hy); 
      for(int i=0;i<=degree;i++) {
	vHy[i]=HermitePol(i,y);
      }
    }
    
    
    double expfact=1./(2*M_PI*sigma*sigma)*exp(-0.5*(x*x+y*y));
    
    double prefactor = 1; // first constant term is not a free parameter, it is = 1 because PSF integral = 1, all higher order gauss hermite terms have integral = 0 
    // we skip constant monomial (0,0), but also (1,0) and (0,1) because degenerate with position
    {
      int param_index=0;
      for(int j=0;j<=degree;j++) {
	int imin=0; 
	if(j==0) {imin=2;} // skip two first (0,0) and (1,0)
	else if(j==1) {imin=1;} // skip  first (0,1)
	
	for(int i=imin;i<=degree;i++,param_index++) {
	  prefactor+=Params[param_index]*Hx[i]*Hy[j];
	} 
      }
    }
    psf_val = prefactor*expfact;
    
    
    if(ParamDer) {
      int param_index=0;
      for(int j=0;j<=degree;j++) {
	int imin=0; 
	if(j==0) {imin=2;} // skip two first (0,0) and (1,0)
	else if(j==1) {imin=1;} // skip  first (0,1)
	
	for(int i=imin;i<=degree;i++,param_index++) {
	  (*ParamDer)[param_index] = expfact*Hx[i]*Hy[j];
	} 
      }
    }
    if(PosDer) { // wrong sign on purpose (derivatives w.r.t -X)
      double& dvdx=(*PosDer)[0];
      double& dvdy=(*PosDer)[1];
      
      dvdx=x/sigma*psf_val;
      dvdy=y/sigma*psf_val;
      
      // should take care of other terms here
      { 
	harp::vector_double& vdHx = const_cast<harp::vector_double&>(dHx); 
	for(int i=0;i<=degree;i++) {
	  vdHx[i]=HermitePolDerivative(i,x);
	}
	harp::vector_double& vdHy = const_cast<harp::vector_double&>(dHy); 
	for(int i=0;i<=degree;i++) {
	  vdHy[i]=HermitePolDerivative(i,y);
	}
      }
      
      double d_poly_dx = 0;
      double d_poly_dy = 0;
      
      /*
      const double* p=params;  
      const double *hy=Hy.Data();
      const double *dhy=dHy.Data();
      for(int j=0;j<=degree;j++,hy++,dhy++) {
	const double *hx=Hx.Data();
	const double *dhx=dHx.Data();
	
	int imin=0; 
	if(j==0) {imin=2; hx+=2; dhx +=2;} // skip two first (0,0) and (1,0)
	else if(j==1) {imin=1; hx+=1; dhx +=1;} // skip  first (0,1)
	
	for(int i=imin;i<=degree;i++,p++,hx++,dhx++) {
	  d_poly_dx+=(*p)*(*dhx)*(*hy);
	  d_poly_dy+=(*p)*(*hx)*(*dhy);
	} 
      }
      */
      
      int param_index=0;
      for(int j=0;j<=degree;j++) {
	int imin=0; 
	if(j==0) {imin=2;} // skip two first (0,0) and (1,0)
	else if(j==1) {imin=1;} // skip  first (0,1)
	for(int i=imin;i<=degree;i++,param_index++) {
	  d_poly_dx += Params[param_index]*dHx[i]*Hy[j];
	  d_poly_dy += Params[param_index]*Hx[i]*dHy[j];
	} 
      }
      
      dvdx -= d_poly_dx*expfact/sigma; // minus sign cause derivative wrt -x
      dvdy -= d_poly_dy*expfact/sigma;
      
    }

  }// end of cut at N sigma
  


#ifdef ADD_Y_TAILS_TO_GAUSS_HERMITE

  if(x*x<25) { // sharp cut at 5 sigma
    
    
    double cross_dispersion_profile = 1./sqrt(2*M_PI)/sigma*exp(-0.5*x*x);
    
    double tail_scale = 0.8; 
    double core_power = 2;
    double tail_power = 1.2;
    
    // define norm at a distance of 20 pixels from peak
    double tail_norm  = 0.0005*pow(tail_scale/20,-tail_power);
    
    
    double ys2   = pow(fabs(input_Y/tail_scale),core_power);
    double ys2p1 = 1+ys2;
    double tail_prof = cross_dispersion_profile*tail_norm*pow(ys2p1,-tail_power/core_power);
    
    double tail_val = Params(y_tail_norm_index)*tail_prof;
    
    
    
    if(PosDer) { // wrong sign on purpose (derivatives w.r.t -X)
      double& dvdx=(*PosDer)[0];
      double& dvdy=(*PosDer)[1];
      dvdx += x/sigma*tail_val;
      int signe=1;
      if(input_Y<0) signe = -1;
      
      dvdy += signe*tail_power/ys2p1/tail_scale*pow(fabs(input_Y/tail_scale),core_power-1)*tail_val;
      
    }
    
    if(ParamDer) {
      (*ParamDer)[y_tail_norm_index]=tail_prof;
    }
    
    psf_val += tail_val;
  }
#endif // end of ADD_Y_TAILS_TO_GAUSS_HERMITE


#ifdef ADD_2D_TAILS_TO_GAUSS_HERMITE
  
  

  
  // numbers tuned for ccd r1
  double core_size = 0.8; // pixel
  double y_scale_plus   = 1; // by construction , so that y_core_size = core_size


  double tail_norm_param = Params[tail_norm_index];
#ifndef LORENTZIAN_TAILS
  double y_scale_minus   = Params[tail_y_scale_minus_index]; // 0.8; // the other allows to change amplitude of tail lower scale means higher tail flux by scale**-power
  double x_scale_minus   = Params[tail_x_scale_minus_index]; //1
  double x_scale_plus    = Params[tail_x_scale_plus_index]; //0.6
  double tail_power      = Params[tail_power_index];
#endif
  


#ifdef LORENTZIAN_TAILS 
  double y_scale_minus  = 1.0;
  double x_scale_minus  = 1.0;
  double x_scale_plus   = 1.0;
#endif

  double x_scale = x_scale_plus;
  if(input_X<0) x_scale = -x_scale_minus;
  double y_scale = y_scale_plus;
  if(input_Y<0) y_scale = -y_scale_minus;
  
#ifndef LORENTZIAN_TAILS 
  double core_power = 2.;
#endif
  
  

  // define norm at a distance of 20 pixels from peak along y
  double ref_dist_inv = 1./20;
  
  core_size    *= ref_dist_inv;
  x_scale      *= ref_dist_inv;
  y_scale      *= ref_dist_inv;

#ifdef LORENTZIAN_TAILS 
  double rs2p1 = square(core_size)+square(input_X*x_scale)+square(input_Y*y_scale);
  double tail_val_ref = 0.00001/rs2p1;
#else 
  double rs2p1 = pow(core_size,core_power)+pow(input_X*x_scale,core_power)+pow(input_Y*y_scale,core_power);
  double tail_val_ref = 0.00001*pow(rs2p1,-tail_power/core_power);
#endif
  
  double tail_val = tail_norm_param*tail_val_ref;
  
  if(PosDer) {

    double& dvdx=(*PosDer)[0];
    double& dvdy=(*PosDer)[1];
    
    
#ifdef LORENTZIAN_TAILS 
    dvdx += 2/rs2p1*x_scale*(input_X*x_scale)*tail_val;
    dvdy += 2/rs2p1*y_scale*(input_Y*y_scale)*tail_val;
#else   
    dvdx += tail_power/rs2p1*x_scale*pow(input_X*x_scale,core_power-1)*tail_val;
    dvdy += tail_power/rs2p1*y_scale*pow(input_Y*y_scale,core_power-1)*tail_val;
#endif
  }
  
  
  if(ParamDer) {
    
    (*ParamDer)[tail_norm_index]=tail_val_ref; // norm of tail
    
#ifndef LORENTZIAN_TAILS     
    (*ParamDer)[tail_power_index]=-1./core_power*log(rs2p1)*tail_val; // slope of tail
    
    if(input_Y<0) {
      (*ParamDer)[tail_y_scale_minus_index] = tail_power/rs2p1*(ref_dist_inv*input_Y)*pow(input_Y*y_scale,core_power-1)*tail_val;
    }else{
      (*ParamDer)[tail_y_scale_minus_index] = 0;
    }
    if(input_X<0) { 
      (*ParamDer)[tail_x_scale_minus_index] = tail_power/rs2p1*(ref_dist_inv*input_X)*pow(input_X*x_scale,core_power-1)*tail_val;
      (*ParamDer)[tail_x_scale_plus_index] = 0;
    }else{
      (*ParamDer)[tail_x_scale_plus_index] = -tail_power/rs2p1*(ref_dist_inv*input_X)*pow(input_X*x_scale,core_power-1)*tail_val;
      (*ParamDer)[tail_x_scale_minus_index] = 0;
    }
#endif


  }
  
  //if(tail_val!=0) cout << "tail_val psf_val = " << tail_val << " " << psf_val  << endl;

  psf_val += tail_val;

#endif

  return psf_val;
}

  
harp::vector_double specex::GaussHermitePSF::DefaultParams() const 
{
  
  harp::vector_double Params(NPar());

  // all = zero at beginning = a pure gaussian
#ifdef ADD_Y_TAILS_TO_GAUSS_HERMITE
  Params[tail_norm_index] = 1;
#endif
#ifdef ADD_2D_TAILS_TO_GAUSS_HERMITE
  Params[tail_norm_index]  = 0; // norm

#ifndef LORENTZIAN_TAILS
  Params[tail_power_index] = 2; // slope
  Params[tail_x_scale_minus_index] = 1;
  Params[tail_x_scale_plus_index] = 1;
  Params[tail_y_scale_minus_index] = 1;
#endif
  
#endif

  return Params;
}
/*
  Fits format
  


*/

#include <specex_image_data.h>

void specex::GaussHermitePSF::WriteFits(fitsfile* fp, int first_hdu) const {
  
  ////////////////////////////
  string PSFVER = "2";
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
    int fiber_id = fiber_it->first;
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
    int LDEGX = params_of_bundle.Polynomials[0].xdeg;
    int LDEGY = params_of_bundle.Polynomials[0].ydeg;
    double LXMIN = params_of_bundle.Polynomials[0].xmin;
    double LXMAX = params_of_bundle.Polynomials[0].xmax;
    double LYMIN = params_of_bundle.Polynomials[0].ymin;
    double LYMAX = params_of_bundle.Polynomials[0].ymax;
    
    // check all params_of_bundle have same degree !!
    for(size_t i=1;i<params_of_bundle.Polynomials.size();i++) {
      if(LDEGX != params_of_bundle.Polynomials[i].xdeg) SPECEX_ERROR("Need same degree for all Legendre2DPol for complying with fits format");
      if(LDEGY != params_of_bundle.Polynomials[i].ydeg) SPECEX_ERROR("Need same degree for all Legendre2DPol for complying with fits format");
      if(LXMIN != params_of_bundle.Polynomials[i].xmin) SPECEX_ERROR("Need same xy range for all Legendre2DPol for complying with fits format");
      if(LXMAX != params_of_bundle.Polynomials[i].xmax) SPECEX_ERROR("Need same xy range for all Legendre2DPol for complying with fits format");
      if(LYMIN != params_of_bundle.Polynomials[i].ymin) SPECEX_ERROR("Need same xy range for all Legendre2DPol for complying with fits format");
      if(LYMAX != params_of_bundle.Polynomials[i].ymax) SPECEX_ERROR("Need same xy range for all Legendre2DPol for complying with fits format"); 
    }
    
    int GHDEGX  = degree; 
    int GHDEGY  = degree;
    double GHSIGX = sigma;
    double GHSIGY = sigma;
    
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
	
	const harp::vector_double& legendre_coefficents = params_of_bundle.Polynomials[ gh_index ].coeff;
	gh_index++;

	// now loop on legendre coefficients m(i+j*(xdeg+1))
	for(size_t i=0;i<legendre_coefficents.size();i++,buffer_index++) {
	  buffer[buffer_index]=legendre_coefficents[i];
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
    
    long fpixels[4];
    fpixels[0] = 1;
    fpixels[1] = 1;
    fpixels[2] = 1;
    fpixels[3] = 1;
    
    fits_create_img ( fp, harp::fits::ftype< double >::bitpix(), 4, naxes, &status ); harp::fits::check ( status );
    fits_write_pix ( fp, harp::fits::ftype< double >::datatype(), fpixels, ncoefs, buffer , &status ); harp::fits::check ( status );
    harp::fits::key_write(fp,"FIBERMIN",(long long int)FIBERMIN,"First fiber for which this PSF is valid");
    harp::fits::key_write(fp,"FIBERMAX",(long long int)FIBERMAX,"Last fiber for which this PSF is valid");
    
    harp::fits::key_write(fp,"HSIZEX",(long long int)hSizeX,"Half size of PSF in fit, NX=2*HSIZEX+1");
    harp::fits::key_write(fp,"HSIZEY",(long long int)hSizeY,"Half size of PSF in fit, NY=2*HSIZEY+1");
    
    harp::fits::key_write(fp,"GHSIGX",(long long int)GHSIGX,"Gauss-Hermite sigma along x");
    harp::fits::key_write(fp,"GHSIGY",(long long int)GHSIGY,"Gauss-Hermite sigma along y");
    
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

void specex::GaussHermitePSF::ReadFits(fitsfile* fp, int first_hdu) {

  SPECEX_INFO("starting reading GaussHermitePSF"); 

  int status = 0;  
  fits_movabs_hdu ( fp, first_hdu, NULL, &status ); harp::fits::check ( status );
  
  // read 
  { 
    string psf_type;
    harp::fits::key_read (fp,"PSFTYPE", psf_type);
    if(psf_type != "GAUSS-HERMITE") SPECEX_ERROR("PSF in file is not 'GAUSS-HERMITE' but '" << psf_type << "'");
  }
  { 
    string psf_version;
    harp::fits::key_read (fp,"PSFVER", psf_version);
    if(psf_version != "2") SPECEX_ERROR("Cannot read GAUSS-HERMITE PSF format version '" << psf_version << "' only version '2' is implemented so far");
  }
  
  
  
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

    // change HDU
    {
      int status = 0;
      fits_movrel_hdu ( fp, 1, NULL, &status ); harp::fits::check ( status );
    }

    long long int FIBERMIN;
    long long int FIBERMAX; 
    long long int HSIZEX;
    long long int HSIZEY;
    long long int GHSIGX;
    long long int GHSIGY;
    long long int GHDEGX;
    long long int GHDEGY;
    long long int LDEGX;
    long long int LDEGY;
    long long int LXMIN;
    long long int LXMAX;
    long long int LYMIN;
    long long int LYMAX;
    string extname;

    harp::fits::key_read(fp,"FIBERMIN",FIBERMIN);
    harp::fits::key_read(fp,"FIBERMAX",FIBERMAX);   
    harp::fits::key_read(fp,"HSIZEX",HSIZEX);
    harp::fits::key_read(fp,"HSIZEY",HSIZEY);
    harp::fits::key_read(fp,"GHSIGX",GHSIGX);
    harp::fits::key_read(fp,"GHSIGY",GHSIGY);
    harp::fits::key_read(fp,"GHDEGX",GHDEGX);
    harp::fits::key_read(fp,"GHDEGY",GHDEGY);
    harp::fits::key_read(fp,"LDEGX",LDEGX);
    harp::fits::key_read(fp,"LDEGY",LDEGY);
    harp::fits::key_read(fp,"LXMIN",LXMIN);
    harp::fits::key_read(fp,"LXMAX",LXMAX);
    harp::fits::key_read(fp,"LYMIN",LYMIN);
    harp::fits::key_read(fp,"LYMAX",LYMAX);    
    harp::fits::key_read(fp,"EXTNAME",extname);
    
    if(bundle_id == bundle_min ) {
      // set global PSF properties
      hSizeX = HSIZEX;
      hSizeY = HSIZEY;
      degree = GHDEGX;
      if(GHDEGY != GHDEGX) {
	SPECEX_ERROR("currently only same Gauss-Hermite degree in X and Y is implemented");
      }
      sigma = GHSIGX;
      if(GHSIGY != GHSIGX) {
	SPECEX_ERROR("currently only same Gauss-Hermite sigma in X and Y is implemented");
      }
      
      first_fiber_in_file = FIBERMIN;    
      
    }else{ 
      // check they are the same or die
      if(HSIZEX != hSizeX) SPECEX_ERROR("currenlty same PSF size for all bundles described in file");
      if(HSIZEY != hSizeY) SPECEX_ERROR("currenlty same PSF size for all bundles described in file");
      if(GHDEGX != degree) SPECEX_ERROR("currenlty same PSF Gauss-Hermite degree for all bundles described in file");
      if(GHDEGY != degree) SPECEX_ERROR("currenlty same PSF Gauss-Hermite degree for all bundles described in file");
      if(GHSIGX != degree) SPECEX_ERROR("currenlty same PSF Gauss-Hermite sigma for all bundles described in file");
      if(GHSIGY != degree) SPECEX_ERROR("currenlty same PSF Gauss-Hermite sigma for all bundles described in file");
    }
    
    // and now read a 4 dimension fits image
    long naxes[4];
    naxes[0] = LDEGX+1;
    naxes[1] = LDEGY+1;
    naxes[2] = GHDEGX+1;
    naxes[3] = GHDEGY+1;
    
    long fpixels[4];
    fpixels[0] = 1;
    fpixels[1] = 1;
    fpixels[2] = 1;
    fpixels[3] = 1;
    
    size_t ncoefs_gauss_hermite = (GHDEGX+1)*(GHDEGY+1);
    size_t ncoefs_legendre = (LDEGX+1)*(LDEGY+1);
    long ncoefs  = ncoefs_gauss_hermite*ncoefs_legendre;
    
    //fits_create_img ( fp, harp::fits::ftype< double >::bitpix(), 4, naxes, &status ); harp::fits::check ( status );
    //fits_write_pix ( fp, harp::fits::ftype< double >::datatype(), fpixels, ncoefs, buffer , &status ); harp::fits::check ( status );
    // int ffgpv(fitsfile *fptr, int  datatype, LONGLONG firstelem, LONGLONG nelem,
    //   void *nulval, void *array, int *anynul, int  *status);
    double nullval = 0;
    int anynull;
    int status;
    double* buffer = new double[ncoefs];
    for(size_t i=0;i<ncoefs;i++) buffer[i]=0; // set to 0
    fits_read_img(fp, harp::fits::ftype< double >::datatype(), 1, ncoefs, &nullval, buffer, &anynull, &status);
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
	x[i] = img_x(i,fiber_id-first_fiber_in_file);
	y[i] = img_y(i,fiber_id-first_fiber_in_file);
	w[i] = img_wave(i,fiber_id-first_fiber_in_file);
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
    int gh_index=0;
    int buffer_index=0;
    for(int j_gh=0;j_gh<=GHDEGY;j_gh++) { // loop on j=y=rows first
      for(int i_gh=0;i_gh<=GHDEGX;i_gh++) {
	// skip some coeffs
	if(i_gh == 0 && j_gh == 0) {buffer_index += ncoefs_legendre; continue;}
	if(i_gh == 1 && j_gh == 0) {buffer_index += ncoefs_legendre; continue;}
	if(i_gh == 0 && j_gh == 1) {buffer_index += ncoefs_legendre; continue;}
	
	Legendre2DPol pol;
	pol.xdeg = LDEGX;
	pol.ydeg = LDEGY;
	pol.xmin = LXMIN;
	pol.xmax = LXMAX;
	pol.ymin = LYMIN;
	pol.ymax = LYMAX;
	pol.coeff.resize((pol.xdeg+1)*(pol.ydeg+1));
	
	// now loop on legendre coefficients m(i+j*(xdeg+1))
	for(size_t i=0;i<pol.coeff.size();i++,buffer_index++) {
	  pol.coeff[i] = buffer[buffer_index];
	}
	params_of_bundle.Polynomials.push_back(pol);
      }
    }
    
    
    delete [] buffer;
    
  
    SPECEX_INFO("read successfully fiber bundle " << bundle_id);
  }
  
}
 
