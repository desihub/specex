#include <iostream>
#include <fstream>

#include "specex_linalg.h"
#include "specex_legendre.h"
#include "specex_message.h"
#include "specex_fits.h" // for debugging


//! base of (not-normalized) Legendre polynomials
static double LegendrePol(const int Degree, const double &X)
{
  switch (Degree)
    {
    case 0: return 1; break;
    case 1: return X; break;
    case 2: return 0.5*(3*X*X-1.); break;
    case 3: return 0.5*(5*X*X-3)*X; break;
    case 4: return 0.125*((35*X*X-30)*X*X+3); break;
    case 5: return 0.125*(((63*X*X-70)*X*X+15)*X); break;
    case 6: return (((231*X*X-315)*X*X+105)*X*X-5)/16.; break;
    default : 
      return ((2*Degree-1)*X*LegendrePol(Degree-1, X) - 
	      (Degree-1)*LegendrePol(Degree-2, X))/double(Degree);
      /*
      FatalError(" consider implementing Legendre pol of abitrary degree");
      return 1; 
      */
      break;
    }
}


specex::Legendre1DPol::Legendre1DPol(int i_deg, const double& i_xmin, const double& i_xmax) : 
  name("undefined"),
  deg(i_deg),
  xmin(i_xmin),
  xmax(i_xmax)
{
  coeff.resize(deg+1);
  coeff.clear();
}



unhrp::vector_double specex::Legendre1DPol::Monomials(const double &x) const {

  // range is -1,1 if  xmin<x<xmax
  double rx= 2*(x-xmin)/(xmax-xmin)-1;
  unhrp::vector_double m(deg+1);
  for(int i=0;i<=deg;i++) {
    m[i]=LegendrePol(i,rx);
  }
  return m;
}


double specex::Legendre1DPol::Value(const double &x) const {
  return specex::dot(coeff,Monomials(x));
}

bool specex::Legendre1DPol::Fit(const unhrp::vector_double& X, const unhrp::vector_double& Y, const unhrp::vector_double* Yerr, bool set_range) {
   // fit x
  
  if(X.size() != Y.size()) SPECEX_ERROR("Legendre1DPol::Fit, not same size X:" << X.size() << " Y:" << Y.size());
  if(Yerr!=0 && X.size() != Yerr->size()) SPECEX_ERROR("Legendre1DPol::Fit, not same size");
  
  int npar = deg+1;
  int ndata = X.size();
  
  if(set_range) {
    specex::minmax(X,xmin,xmax);
  }
  
  unhrp::matrix_double A(npar,npar); A.clear();
  unhrp::vector_double B(npar); B.clear();
  
  
  
  for(int i=0;i<ndata;i++) {
    double w=1;
    if(Yerr) {
      w=1./square((*Yerr)[i]);
    }
    
    unhrp::vector_double h=Monomials(X[i]);
    specex::syr(w,h,A); // A += w*Mat(h)*h.transposed();
    specex::axpy(double(w*Y[i]),h,B); //B += (w*Y[i])*h;
  }

  //unhrp::matrix_double As=A;
  int status = cholesky_solve(A,B);
  if(status != 0) {
    //write_new_fits_image("A.fits",As);
    
    if(0) {
      for(int i=0;i<ndata;i++) {
	double w=1;
	if(Yerr) {
	  w=1./square((*Yerr)[i]);
	}
	cout << i << " " << X[i] << " " << Y[i] << " " << w << endl;
      }
    }
    SPECEX_ERROR("Legendre1DPol::Fit cholesky_solve failed with status " << status << " deg= " << deg << " xmin=" << xmin << " xmax=" << xmax);
  }  
  coeff=B;
  
  // SPECEX_INFO("successful Legendre1DPol::Fit");

  return true;
}

specex::Legendre1DPol specex::Legendre1DPol::Invert(int add_degree) const {
    
  specex::Legendre1DPol inverse;
  inverse.deg = deg+add_degree;
  inverse.coeff.resize(inverse.deg+1);
  inverse.coeff.clear();
  int npar = inverse.deg + 1;
  int ndata = npar*4;  // 
  double dx = (xmax-xmin)/ndata;
  unhrp::vector_double X(ndata);
  unhrp::vector_double Y(ndata);
  for(int i=0;i<ndata;i++) {
    X[i] = xmin+i*dx;
    Y[i] = Value(X[i]);
  }
  bool ok = inverse.Fit(Y,X,0,true);
  if(!ok) abort();
  return inverse;
}

specex::Legendre1DPol specex::composed_pol(const specex::Legendre1DPol& pol1, const specex::Legendre1DPol& pol2) {
  specex::Legendre1DPol composed;
  composed.deg = max(pol1.deg,pol2.deg);
  composed.coeff.resize(composed.deg+1);
  composed.coeff.clear();
  int npar = composed.deg + 1;
  int ndata = npar*4;  // 
  double dx = (pol2.xmax-pol2.xmin)/ndata;
  unhrp::vector_double X2(ndata);
  unhrp::vector_double Y1(ndata);
  for(int i=0;i<ndata;i++) {
    X2[i] = pol2.xmin+i*dx;
    Y1[i] = pol1.Value(pol2.Value(X2[i]));
  }
  bool ok = composed.Fit(X2,Y1,0,true);
  if(!ok) abort();
  return composed;
}

//============================

specex::Legendre2DPol::Legendre2DPol(int i_xdeg, const double& i_xmin, const double& i_xmax,
			     int i_ydeg, const double& i_ymin, const double& i_ymax
			     ) : 
  name("undefined"),
  xdeg(i_xdeg),
  xmin(i_xmin),
  xmax(i_xmax),
  ydeg(i_ydeg),
  ymin(i_ymin),
  ymax(i_ymax)
{
  Fill();
}

void specex::Legendre2DPol::Fill() {
  coeff.resize((xdeg+1)*(ydeg+1));
  coeff.clear();
}
 
unhrp::vector_double specex::Legendre2DPol::Monomials(const double &x, const double &y) const {
  
  // range is -1,1 if  xmin<x<xmax
  double rx= 2*(x-xmin)/(xmax-xmin)-1;
  double ry= 2*(y-ymin)/(ymax-ymin)-1;
  
  unhrp::vector_double mx(xdeg+1); mx.clear();
  for(int i=0;i<=xdeg;i++)
    mx[i]=LegendrePol(i,rx);
  
  unhrp::vector_double m((xdeg+1)*(ydeg+1));
  int index=0;
  for(int j=0;j<=ydeg;j++) {
    double myj = LegendrePol(j,ry);
    for(int i=0;i<=xdeg;i++,index++) {
      m[index]=mx[i]*myj;
    }
  }
  return m;
}


double specex::Legendre2DPol::Value(const double &x,const double &y) const {
  return specex::dot(coeff,Monomials(x,y));
}

//============================

specex::SparseLegendre2DPol::SparseLegendre2DPol(int i_xdeg, const double& i_xmin, const double& i_xmax,
			     int i_ydeg, const double& i_ymin, const double& i_ymax
			     ) : 
  name("undefined"),
  xdeg(i_xdeg),
  xmin(i_xmin),
  xmax(i_xmax),
  ydeg(i_ydeg),
  ymin(i_ymin),
  ymax(i_ymax)
{
}


void specex::SparseLegendre2DPol::Add(int i,int j) {
  if(i<0 || i>xdeg || j<0 || j>ydeg) {SPECEX_ERROR("in SparseLegendre2DPol::Add not valid indices " << i << " " << j);}
  
  int index=i+j*(xdeg+1);
  // checking
  for(std::vector<int>::const_iterator k = non_zero_indices.begin(); k!=  non_zero_indices.end(); k++) {
    if(*k == index) {
      SPECEX_WARNING("in SparseLegendre2DPol::Add, indices " << i << " " << j << " already set");
      return;
    }
  }
  non_zero_indices.push_back(index);
  coeff.resize(non_zero_indices.size());
  coeff.clear();
}
 
void specex::SparseLegendre2DPol::Fill(bool sparse) {
  non_zero_indices.clear();
  if ( ! sparse ) {
    for(int k=0;k<(xdeg+1)*(ydeg+1);k++) non_zero_indices.push_back(k);
    coeff.resize(non_zero_indices.size());
    coeff.clear();
  }else{
    for(int i=0;i<=xdeg;i++) { // x coordinate
      for(int j=0;j<=ydeg;j++) { // wave coordinate	      
	if(i==0) {
	  Add(i,j); // full wavelength resolution
	}else if(i==1) {
	  if(j<2) Add(i,j); // only first x * wavelength cross-term
	}else{
	  if(j==0) Add(i,j); // only  x terms
	}	    
      }
    }
  }
}

void specex::SparseLegendre2DPol::Clear() {
  non_zero_indices.clear();
  coeff.resize(0);
}


unhrp::vector_double specex::SparseLegendre2DPol::Monomials(const double &x, const double &y) const {
  
  // range is -1,1 if  xmin<x<xmax
  double rx= 2*(x-xmin)/(xmax-xmin)-1;
  double ry= 2*(y-ymin)/(ymax-ymin)-1;
  
  
  unhrp::vector_double m(non_zero_indices.size()); m.clear();
  int index=0;
  for(std::vector<int>::const_iterator k = non_zero_indices.begin(); k!=  non_zero_indices.end(); k++, index++) {
    int i = (*k)%(xdeg+1);
    int j = (*k)/(xdeg+1);
    m[index]=LegendrePol(i,rx)*LegendrePol(j,ry);
  }
  return m;
}


double specex::SparseLegendre2DPol::Value(const double &x,const double &y) const {
  return specex::dot(coeff,Monomials(x,y));
}

