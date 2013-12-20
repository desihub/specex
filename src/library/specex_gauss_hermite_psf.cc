
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
  
  if(x*x+y*y<25) { // sharp cut at 5 sigma

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
  
  
void specex::GaussHermitePSF::InitParams(const double &i_sigma, harp::vector_double &Params) 
{
  sigma = i_sigma;
  Params.resize(NPar());
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
}
/*
  Fits format
  


*/

#include <specex_image_data.h>

void specex::GaussHermitePSF::WriteFits(fitsfile* fp, int first_hdu) const {
  
  string PSFVER = "2";
  
  
  int LDEGX = Params[0].xdeg;
  int LDEGY = Params[0].ydeg;
  double LXMIN = Params[0].xmin;
  double LXMAX = Params[0].xmax;
  double LYMIN = Params[0].ymin;
  double LYMAX = Params[0].ymax;
  
  
  
  for(size_t i=1;i<Params.size();i++) {
    if(LDEGX != Params[i].xdeg) SPECEX_ERROR("Need same degree for all Legendre2DPol for complying with fits format");
    if(LDEGY != Params[i].ydeg) SPECEX_ERROR("Need same degree for all Legendre2DPol for complying with fits format");
    if(LXMIN != Params[i].xmin) SPECEX_ERROR("Need same xy range for all Legendre2DPol for complying with fits format");
    if(LXMAX != Params[i].xmax) SPECEX_ERROR("Need same xy range for all Legendre2DPol for complying with fits format");
    if(LYMIN != Params[i].ymin) SPECEX_ERROR("Need same xy range for all Legendre2DPol for complying with fits format");
    if(LYMAX != Params[i].ymax) SPECEX_ERROR("Need same xy range for all Legendre2DPol for complying with fits format");
    
  }
  
  int NPIX_X = ccd_image_n_cols;
  int NPIX_Y = ccd_image_n_rows;
  
  
  int NFLUX = NPIX_Y; // this has to be a parameter
  int NSPEC = FiberTraces.size(); // this also has to be a parameter
  double ystep = (LYMAX-LYMIN)/NFLUX; // this has to be a parameter
  
  
  SPECEX_INFO("NFLUX=" << NFLUX << " NSPEC=" << NSPEC);
	      
  specex::image_data x_data(NFLUX,NSPEC);
  specex::image_data y_data(NFLUX,NSPEC);
  specex::image_data wave_data(NFLUX,NSPEC);
  
  

  
  int ispec=0;
  for(std::map<int,Trace>::const_iterator fiber_it = FiberTraces.begin(); fiber_it != FiberTraces.end(); ++fiber_it, ++ispec) {
    int fiber_id = fiber_it->first;
    const specex::Trace& trace = fiber_it->second;
    
    for(int i=0;i<NFLUX;i++) {
      double y = LYMIN+ystep*i;

      y_data(i,ispec)=y;
      x_data(i,ispec)=trace.X_vs_Y.Value(y);
      wave_data(i,ispec)=pow(10.,trace.lW_vs_Y.Value(y));
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
  
  harp::fits::key_write(fp,"PSFPARAM","X",default_comment);
  
     
  harp::fits::img_append < double > ( fp, NSPEC, NFLUX );
  harp::fits::img_write (fp, y_data.data);
  harp::fits::key_write(fp,"PSFPARAM","Y",default_comment);
  
  
  harp::fits::img_append < double > ( fp, NSPEC, NFLUX );
  harp::fits::img_write (fp, wave_data.data);
  harp::fits::key_write(fp,"EXTNAME","WAVELENGTH",default_comment);
  
  // now create image of coefficients
  // first check all polynomials have same degree
  
  int GHDEGX  = degree; 
  int GHDEGY  = degree;
  double GHSIGX = sigma;
  double GHSIGY = sigma;
  
  int FIBERMIN = FiberTraces.begin()->first;
  int FIBERMAX = (--FiberTraces.end())->first;
  

  size_t ncoefs_gauss_hermite = (GHDEGX+1)*(GHDEGY+1);
  size_t ncoefs_legendre = (LDEGX+1)*(LDEGY+1);
  
  long ncoefs  = ncoefs_gauss_hermite*ncoefs_legendre;
 
  // GHX = AXIS1 ?
  // GHY = AXIS2 ?
  // LX  = AXIS3 ?
  // LY  = AXIS4 ?
  // index = A1 + NAXIS1*(A2 + NAXIS2*(A3 + NAXIS3*(A4))) ?

  
  double* buffer = new double[ncoefs];

  // loop on gauss-hermite parameters
  int index=0;

  // first 3 gauss hermite coefficients are by default GH00=1, GH01=0, GH10=0
  buffer[index++]=1; // order 0 legendre pol for GH00 is = 1
  for(;index<3*ncoefs_legendre; index++) {
    buffer[index] = 0;
  }
  
  for(size_t gh_param_index = 0; gh_param_index < Params.size(); ++ gh_param_index) {
    const harp::vector_double& legendre_coefficents = Params[ gh_param_index ].coeff;
    for(size_t i=0;i<legendre_coefficents.size();i++)
      buffer[index++]=legendre_coefficents[i];
  }
  
  // and now create a 4 dimension fits image
  long naxes[4];
  naxes[0] = GHDEGX+1;
  naxes[1] = GHDEGY+1;
  naxes[2] = LDEGX+1;
  naxes[3] = LDEGY+1;
  
  long fpixels[4];
  fpixels[0] = 1;
  fpixels[1] = 1;
  fpixels[2] = 1;
  fpixels[3] = 1;
  
  fits_create_img ( fp, harp::fits::ftype< double >::bitpix(), 4, naxes, &status ); harp::fits::check ( status );
  fits_write_pix ( fp, harp::fits::ftype< double >::datatype(), fpixels, ncoefs, buffer , &status ); harp::fits::check ( status );
  harp::fits::key_write(fp,"FIBERMIN",(long long int)FIBERMIN,"First fiber for which this PSF is valid");
  harp::fits::key_write(fp,"FIBERMAX",(long long int)FIBERMAX,"Last fiber for which this PSF is valid");
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
  
  delete [] buffer;
  
}

void specex::GaussHermitePSF::ReadFits(fitsfile* fp, int first_hu) {
  
}
 
