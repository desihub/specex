#include <cmath>
#include <iostream>
#include <fstream>

#include <harp.hpp>

#include <specex_spot.h>
#include <specex_trace.h>
#include <specex_linalg.h>
#include <specex_message.h>
#include <specex_spot_array.h>

specex::Trace::Trace(int i_fiber) :
  fiber(i_fiber)
{
  
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
    
    harp::matrix_double A(npar,npar); specex::zero(A);
    harp::vector_double B(npar); specex::zero(B);
    
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
    harp::matrix_double A(npar,npar); specex::zero(A);
    harp::vector_double B(npar); specex::zero(B);
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
  double x_chi2 = 0;
  double x_sumw = 0;
  double x_sumwx  = 0;
  double x_sumwx2 = 0;
  double y_chi2 = 0;
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
  return true;
}


/* 
void specex::Trace::write(ostream &os) const {
  os << "Beginspecex::Trace" << endl;
  os << fiber << endl;
  X_vs_W.write(os);
  Y_vs_W.write(os);
  os << "Endspecex::Trace" << endl;
}



#define ERROR_MESSAGE {cout << "ERROR specex::Trace::read failed" << endl; abort();}

bool specex::Trace::read(istream& is) {
  string label;
  if(! (is >> label)) ERROR_MESSAGE;
  if(! (label=="Beginspecex::Trace")) ERROR_MESSAGE;

  if(! (is >> fiber)) ERROR_MESSAGE;
  if(! X_vs_W.read(is)) ERROR_MESSAGE;
  if(! Y_vs_W.read(is)) ERROR_MESSAGE;
  
  if(! (is >> label)) ERROR_MESSAGE;
  if(! (label=="Endspecex::Trace")) ERROR_MESSAGE;
  //cout << "INFO specex::Trace::read ok" << endl;
  return true;
}

#undef ERROR_MESSAGE

void specex::Trace::write(const std::string &FileName) const {
  ofstream os(FileName.c_str());
  write(os);
  os.close();
}

bool specex::Trace::read(const std::string &FileName) {
  ifstream is(FileName.c_str());
  bool ok = read(is);
  is.close();
  return ok;
}
*/ 

/* 
#include <fitstable.h>
bool specex::TraceSet::ReadSDSS_SingleSet_Fits(const string& filename, int hdu_number, TraceSetType ttype) {
  FitsTable table(filename,hdu_number);

  // test formatting
  if(table.data.empty()) {cout << "ERROR specex::Trace::read_sdss_fits data empty" << endl; return false;}

  double yjumplo=-12;
  double yjumphi=-12;
  double yjumpval=-12;

  std::vector<FitsTableEntry>& data = table.data[0];
  int data_size = data.size();
  if(! (data_size == 4 || data_size == 7)) {cout << "ERROR specex::Trace::read_sdss_fits data has " << data.size() << " rows and I expect 4 or 7" << endl;  return false;}
  if(data[0].string_val != "legendre") {cout << "ERROR specex::Trace::read_sdss_fits func is '" << data[0].string_val << "' and only legendre is implemented" << endl; return false;}
  if(data[1].double_vals.size()  != 1 ) {cout << "ERROR specex::Trace::read_sdss_fits col. 2 has wrong size (not 1)" << endl; return false;}
  if(data[2].double_vals.size()  != 1 ) {cout << "ERROR specex::Trace::read_sdss_fits col. 3 has wrong size (not 1)" << endl; return false;}
  if(data_size==7) {
    if(data[4].double_vals.size()  != 1 ) {cout << "ERROR specex::Trace::read_sdss_fits col. 5 has wrong size (not 1)" << endl; return false;}
    if(data[5].double_vals.size()  != 1 ) {cout << "ERROR specex::Trace::read_sdss_fits col. 6 has wrong size (not 1)" << endl; return false;}
    if(data[6].double_vals.size()  != 1 ) {cout << "ERROR specex::Trace::read_sdss_fits col. 7 has wrong size (not 1)" << endl; return false;}
  
    
    
    for(int i=4;i<7;i++) {
      const string& colname = table.columns[i].name;
      FitsTableEntry& entry = data[i];
      if(colname == "XJUMPLO") {
	yjumplo = entry.double_vals(0);
      }else if(colname == "XJUMPHI") {
	yjumphi = entry.double_vals(0); 
      }else if(colname == "XJUMPVAL") {
	yjumpval = entry.double_vals(0);
      }else{
	cout << "ERROR specex::Trace::read_sdss_fits dont' know what to do with column named '" << colname << "'" << endl; return false;
      }
    }
    cout << "INFO specex::Trace::read_sdss_fits, has y-jump : lo,hi,val = " << yjumplo << "," << yjumphi << "," << yjumpval << endl;
    if(yjumplo==-12 || yjumphi==-12 || yjumpval==-12) {
      cout << "ERROR specex::Trace::read_sdss_fits, missing something" << endl;
      return false;
    }
  }


  int ntotcoeff = data[3].double_vals.size();
  int ncoeff_per_fiber = ntotcoeff/NUMBER_OF_FIBERS_PER_CCD;
  if(ntotcoeff%NUMBER_OF_FIBERS_PER_CCD!=0) {
    cout << "ERROR specex::Trace::read_sdss_fits col. 3 has bizarre size " << ntotcoeff << endl; 
    return false;
  }
  
  
  if(size()!=NUMBER_OF_FIBERS_PER_CCD) resize(NUMBER_OF_FIBERS_PER_CCD);
  
  const double& xmin = data[1].double_vals(0);
  const double& xmax = data[2].double_vals(0);
  const Vect & coefficients = data[3].double_vals;
  
  int index=0;
  for(int fiber=0;fiber<NUMBER_OF_FIBERS_PER_CCD; fiber++) {

    specex::Trace& trace = (*this)[fiber];
    trace.fiber=fiber;
    
    trace.yjumplo = yjumplo;
    trace.yjumphi = yjumphi;
    trace.yjumpval = yjumpval;
    
    Legendre1DPol* pol = 0;
    switch(ttype) {
    case WY : pol = &trace.W_vs_Y; break;
    case XY : pol = &trace.X_vs_Y; break;
    case YW : pol = &trace.Y_vs_W; break;
    case XW : pol = &trace.X_vs_W; break;
    default : {cout << "unknown TraceSetType " << ttype << endl; return false;};
    }
    
    pol->xmin = xmin;
    pol->xmax = xmax;
    pol->deg = ncoeff_per_fiber-1;
    pol->coeff.resize(ncoeff_per_fiber);
      
    for(int i=0;i<ncoeff_per_fiber;i++,index++) 
      pol->coeff(i) = coefficients(index);
  }

  cout <<"DEBUG specex::TraceSet::ReadSDSS_SingleSet_Fits successful" << endl;
  return true;
}


bool specex::TraceSet::ReadSDSS_FullSet_Fits(const string& wave_vs_y_filename, int wave_vs_y_hdu_number,
					   const string& x_vs_y_filename, int x_vs_y_hdu_number) {
  
  
  cout << "INFO specex::TraceSet::ReadSDSS_FullSet_Fits starting" << endl;
  cout << "INFO Reading " << wave_vs_y_filename << "[" << wave_vs_y_hdu_number << "]" << endl;
  cout << "INFO Reading " << x_vs_y_filename << "[" << x_vs_y_hdu_number << "]" << endl;
  
  
  bool ok;
  ok=ReadSDSS_SingleSet_Fits(wave_vs_y_filename,wave_vs_y_hdu_number,WY);
  if(!ok) return false;
  
  ok=ReadSDSS_SingleSet_Fits(x_vs_y_filename,x_vs_y_hdu_number,XY);
  if(!ok) return false;
  
  for(int fiber=0;fiber<size();fiber++) {
    
    specex::Trace& trace = (*this)[fiber];
    // first need to invert W_vs_Y
    trace.Y_vs_W = trace.W_vs_Y.Invert(2); 
    
    // now compose
    trace.X_vs_W = composed_pol(trace.X_vs_Y,trace.Y_vs_W);   
  }
  cout << "INFO specex::TraceSet::ReadSDSS_FullSet_Fits done" << endl;
  return true;
  
}
*/
