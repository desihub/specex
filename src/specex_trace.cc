#include <cmath>
#include <iostream>
#include <fstream>

#include <harp.hpp>

#include <specex_spot.h>
#include <specex_trace.h>
#include <specex_linalg.h>
#include <specex_message.h>
#include <specex_spot_array.h>
#include <specex_vector_utils.h>

specex::Trace::Trace(int i_fiber) :
  fiber(i_fiber)
{
  synchronized=false;
}

void specex::Trace::resize(int ncoeff) {
  harp::vector_double coeff;

  coeff=X_vs_W.coeff;
  X_vs_W.deg = ncoeff-1;
  X_vs_W.coeff.resize(ncoeff);
  X_vs_W.coeff = ublas::project(coeff,ublas::range(0,min(ncoeff,int(coeff.size()))));
  
  coeff=Y_vs_W.coeff;
  Y_vs_W.deg = ncoeff-1;
  Y_vs_W.coeff.resize(ncoeff);
  Y_vs_W.coeff = ublas::project(coeff,ublas::range(0,min(ncoeff,int(coeff.size()))));
  
  coeff=W_vs_Y.coeff;
  W_vs_Y.deg = ncoeff-1;
  W_vs_Y.coeff.resize(ncoeff);
  W_vs_Y.coeff = ublas::project(coeff,ublas::range(0,min(ncoeff,int(coeff.size()))));
  
  coeff=X_vs_Y.coeff;
  X_vs_Y.deg = ncoeff-1;
  X_vs_Y.coeff.resize(ncoeff);
  X_vs_Y.coeff = ublas::project(coeff,ublas::range(0,min(ncoeff,int(coeff.size()))));
  
  
  
}

bool specex::Trace::Fit(std::vector<specex::Spot_p> spots, bool set_xy_range) {

  if(fiber==-1) {
    SPECEX_ERROR("specex::Trace::Fit need to set fiber id first");
  }

  SPECEX_INFO("specex::Trace::Fit starting fit of trace of fiber " << fiber);

  
  int nspots=0;
  for(size_t s=0;s<spots.size();s++) {
    const specex::Spot &spot = *(spots[s]);
    if(spot.fiber == fiber) 
      nspots++;
  }
  
   
  if(set_xy_range) {
    X_vs_W.xmin = 1e20;
    X_vs_W.xmax = -1e20;
    Y_vs_W.xmin = 1e20;
    Y_vs_W.xmax = -1e20;
    for(size_t s=0;s<spots.size();s++) {
      const specex::Spot &spot = *(spots[s]);
      if(spot.fiber != fiber) 
	continue;
      if(spot.wavelength<X_vs_W.xmin) X_vs_W.xmin=spot.wavelength;
      if(spot.wavelength>X_vs_W.xmax) X_vs_W.xmax=spot.wavelength;
      if(spot.wavelength<Y_vs_W.xmin) Y_vs_W.xmin=spot.wavelength;
      if(spot.wavelength>Y_vs_W.xmax) Y_vs_W.xmax=spot.wavelength;
    }
  }
  if(X_vs_W.xmax == X_vs_W.xmin) X_vs_W.xmax = X_vs_W.xmin + 0.1;
  if(Y_vs_W.xmax == Y_vs_W.xmin) Y_vs_W.xmax = Y_vs_W.xmin + 0.1;
  

  SPECEX_INFO("number of spots for fiber " << fiber << " = " << nspots);
  
  if(nspots==0) SPECEX_ERROR("specex::Trace::Fit : no spots");
  
  X_vs_W.deg = min(SPECEX_TRACE_DEFAULT_LEGENDRE_POL_DEGREE,nspots-1);
  Y_vs_W.deg = min(SPECEX_TRACE_DEFAULT_LEGENDRE_POL_DEGREE,nspots-1);

  //#warning only deg2 trace fit for debug 
  //Y_vs_W.deg = min(2,nspots-1);
   
  

  if(nspots<max(X_vs_W.deg+1,Y_vs_W.deg+1)) {
    SPECEX_ERROR("specex::Trace::Fit only " << nspots << " spots to fit " << X_vs_W.deg 
		 << " and " << Y_vs_W.deg << " degree polynomials");
  }
  

  X_vs_W.coeff.resize(X_vs_W.deg+1);
  Y_vs_W.coeff.resize(Y_vs_W.deg+1);
  
  
  // fit x
  {
    int npar = X_vs_W.coeff.size();
    
    harp::matrix_double A(npar,npar); A.clear();
    harp::vector_double B(npar); B.clear();
    
    for(size_t s=0;s<spots.size();s++) {
      const specex::Spot &spot = *(spots[s]);
      if(spot.fiber != fiber) 
	continue;
      double w=1;
      double res=spot.xc;
      harp::vector_double h=X_vs_W.Monomials(spot.wavelength);
      
      
      specex::syr(w,h,A);  // A += w*Mat(h)*h.transposed();
      specex::axpy(w*res,h,B); // B += (w*res)*h;
    }
    int status = cholesky_solve(A,B);
    if(status != 0) {
      SPECEX_ERROR("failed to fit X vs wavelength");
    }
    X_vs_W.coeff=B;
  }

  // fit y
  {
    int npar = Y_vs_W.coeff.size();
    harp::matrix_double A(npar,npar); A.clear();
    harp::vector_double B(npar); B.clear();
    for(size_t s=0;s<spots.size();s++) {
      const specex::Spot &spot = *(spots[s]);
      if(spot.fiber != fiber) 
	continue;
      double w=1;
      double res=spot.yc;
      harp::vector_double h=Y_vs_W.Monomials(spot.wavelength);
      
      specex::syr(w,h,A);  // A += w*Mat(h)*h.transposed();
      specex::axpy(w*res,h,B); // B += (w*res)*h;
    }
    int status = cholesky_solve(A,B);
    if(status != 0) {
      SPECEX_ERROR("failed to fit Y vs wavelength");
    }
    Y_vs_W.coeff=B;
  }

  // monitoring results
  double x_sumw = 0;
  double x_sumwx  = 0;
  double x_sumwx2 = 0;
  double y_sumw = 0;
  double y_sumwy  = 0;
  double y_sumwy2 = 0;
  for(size_t s=0;s<spots.size();s++) {
    const specex::Spot &spot = *(spots[s]);
    if(spot.fiber != fiber) 
      continue;
    double xw=1;
    double yw=1;

    double xres=spot.xc-X_vs_W.Value(spot.wavelength);
    double yres=spot.yc-Y_vs_W.Value(spot.wavelength);
    
    
    x_sumw += xw;
    x_sumwx += xw*xres;
    x_sumwx2 += xw*xres*xres;
    y_sumw += yw;
    y_sumwy += yw*yres;
    y_sumwy2 += yw*yres*yres;
  }
  double x_mean = x_sumwx/x_sumw;
  double x_rms  = sqrt(x_sumwx2/x_sumw);
  double y_mean = y_sumwy/y_sumw;
  double y_rms  = sqrt(y_sumwy2/y_sumw);
  
  SPECEX_INFO(
    "specex::Trace::Fit fiber trace #" << fiber
    << " nspots=" << nspots << " xdeg=" << X_vs_W.deg << " ydeg=" << Y_vs_W.deg
    << " dx=" << x_mean << " xrms=" << x_rms
    << " dy=" << y_mean << " yrms=" << y_rms
    );
  if(x_rms>2 || y_rms>2) {
    
    
    for(size_t s=0;s<spots.size();s++) {
      specex::Spot& spot = *(spots[s]);
      if(spot.fiber != fiber) 
	continue;
      spot.xc=X_vs_W.Value(spot.wavelength);
      spot.yc=Y_vs_W.Value(spot.wavelength);
      
    }
    
    write_spots_list(spots,"debug_spots.list");
    SPECEX_ERROR("specex::Trace::Fit rms are too large");
  }
  synchronized = true;
  return true;
}
//#warning check meaning of fibermask
bool specex::Trace::Off() const { return mask==3;} // I am guessing here?

int specex::eval_bundle_size(const specex::TraceSet& traceset) {
  
  SPECEX_DEBUG("Guess number of bundles from traces");
  int nfibers = traceset.size();
  
  if(nfibers<30) return nfibers;
  
  double* central_waves = new double[nfibers];
  for(int f=0;f<nfibers;f++)
    central_waves[f]=(traceset.find(f)->second.X_vs_W.xmin+traceset.find(f)->second.X_vs_W.xmax)/2.;
  double central_wave=DConstArrayMedian(central_waves,nfibers);
  SPECEX_DEBUG("Central wavelength             = " << central_wave);
  delete [] central_waves;
  
  double* spacing = new double[nfibers-1];
  for(int f=0;f<nfibers-1;f++) {
    spacing[f] = traceset.find(f+1)->second.X_vs_W.Value(central_wave)-traceset.find(f)->second.X_vs_W.Value(central_wave);
    //SPECEX_DEBUG("Spacing=" << spacing[f]);
  }
  double median_spacing = DConstArrayMedian(spacing,nfibers-1);
  SPECEX_INFO("Median distance between fibers = " << median_spacing);
  
  int number_of_bundles=0;
  int bundle_size=0;
  int first_fiber=0;
    
  for(int f=0;f<nfibers-1;f++) {
    if(spacing[f]>median_spacing*1.5) {
      // we have a bundle
      int current_bundle_size = f-first_fiber+1;
      number_of_bundles += 1;
      SPECEX_DEBUG("Bundle of size " << current_bundle_size);
      
      if(bundle_size==0) {
	bundle_size = current_bundle_size;
      }else{
	if(current_bundle_size != bundle_size) {
	  SPECEX_ERROR("cannot deal with varying bundle size");
	}
      }
      first_fiber=f+1;
    }
  }
  number_of_bundles += 1;
  if(number_of_bundles==1) { // there is only one 
    bundle_size = nfibers;
  }
  
  SPECEX_INFO("number of fibers per bundle    = " << bundle_size);
  
  delete[] spacing;
  return bundle_size;
}
