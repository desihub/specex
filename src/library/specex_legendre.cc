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
  deg(i_deg),
  xmin(i_xmin),
  xmax(i_xmax)
{
  coeff.resize(deg+1);
}



harp::vector_double specex::Legendre1DPol::Monomials(const double &x) const {

  // range is -1,1 if  xmin<x<xmax
  double rx= 2*(x-xmin)/(xmax-xmin)-1;
  
  harp::vector_double m(deg+1);
  for(int i=0;i<=deg;i++) {
    m(i)=LegendrePol(i,rx);
  }
  return m;
}


double specex::Legendre1DPol::Value(const double &x) const {
  return specex::dot(coeff,Monomials(x));
}

/*
void specex::Legendre1DPol::write(std::ostream &os) const {
  os << "BeginLegendre1DPol" << endl;
  os << deg << " " << xmin << " " << xmax << endl;
  coeff.writeASCII(os);
  os << "EndLegendre1DPol" << endl;
  
}
#define ERROR_MESSAGE {cout << "ERROR Legendre1DPol::read failed" << endl; abort();}

bool specex::Legendre1DPol::read(std::istream &is) {
  string label;
  if(! (is >> label)) ERROR_MESSAGE;
  if(! (label=="BeginLegendre1DPol")) ERROR_MESSAGE;

  if(! (is >> deg >> xmin >> xmax)) ERROR_MESSAGE;
  if(coeff.readASCII(is)!=0) ERROR_MESSAGE;
  
  if(int(coeff.size()) != (deg+1)) {
    cout << "ERROR Legendre1DPol::read, coef size doesn't match degree" << endl;
    ERROR_MESSAGE;
  }
  
  if(! (is >> label)) ERROR_MESSAGE;
  if(! (label=="EndLegendre1DPol")) ERROR_MESSAGE;
  
  return true;
}
#undef ERROR_MESSAGE

void specex::Legendre1DPol::write(const std::string &FileName) const {
  ofstream os(FileName.c_str());
  write(os);
  os.close();
}

bool specex::Legendre1DPol::read(const std::string &FileName) {
  ifstream is(FileName.c_str());
  bool ok = read(is);
  is.close();
  return ok;
}
*/

bool specex::Legendre1DPol::Fit(const harp::vector_double& X, const harp::vector_double& Y, const harp::vector_double* Yerr, bool set_range) {
   // fit x
  
  if(X.size() != Y.size()) SPECEX_ERROR("Legendre1DPol::Fit, not same size X:" << X.size() << " Y:" << Y.size());
  if(Yerr!=0 && X.size() != Yerr->size()) SPECEX_ERROR("Legendre1DPol::Fit, not same size");
  
  int npar = deg+1;
  int ndata = X.size();
  
  if(set_range) {
    specex::minmax(X,xmin,xmax);
  }
  
  harp::matrix_double A(npar,npar); A *= 0;
  harp::vector_double B(npar); B *= 0;
  
  
  
  for(int i=0;i<ndata;i++) {
    double w=1;
    if(Yerr) {
      w=1./square((*Yerr)[i]);
    }
    
    harp::vector_double h=Monomials(X[i]);
    specex::syr(w,h,A); // A += w*Mat(h)*h.transposed();
    specex::axpy(double(w*Y[i]),h,B); //B += (w*Y[i])*h;
  }

  //harp::matrix_double As=A;
  int status = cholesky_solve(A,B);
  if(status != 0) {
    //write_new_fits_image("A.fits",As);
    SPECEX_ERROR("Legendre1DPol::Fit cholesky_solve failed with status " << status);
  }  
  coeff=B;
  
  //SPECEX_INFO("successful Legendre1DPol::Fit");

  return true;
}

specex::Legendre1DPol specex::Legendre1DPol::Invert(int add_degree) const {
    
  specex::Legendre1DPol inverse;
  inverse.deg = deg+add_degree;
  inverse.coeff.resize(inverse.deg+1);
  int npar = inverse.deg + 1;
  int ndata = npar*4;  // 
  double dx = (xmax-xmin)/ndata;
  harp::vector_double X(ndata);
  harp::vector_double Y(ndata);
  for(int i=0;i<ndata;i++) {
    X(i) = xmin+i*dx;
    Y(i) = Value(X(i));
  }
  bool ok = inverse.Fit(Y,X,0,true);
  if(!ok) abort();
  return inverse;
}

specex::Legendre1DPol specex::composed_pol(const specex::Legendre1DPol& pol1, const specex::Legendre1DPol& pol2) {
  specex::Legendre1DPol composed;
  composed.deg = max(pol1.deg,pol2.deg);
  composed.coeff.resize(composed.deg+1);
  int npar = composed.deg + 1;
  int ndata = npar*4;  // 
  double dx = (pol2.xmax-pol2.xmin)/ndata;
  harp::vector_double X2(ndata);
  harp::vector_double Y1(ndata);
  for(int i=0;i<ndata;i++) {
    X2(i) = pol2.xmin+i*dx;
    Y1(i) = pol1.Value(pol2.Value(X2(i)));
  }
  bool ok = composed.Fit(X2,Y1,0,true);
  if(!ok) abort();
  return composed;
}

//============================

specex::Legendre2DPol::Legendre2DPol(int i_xdeg, const double& i_xmin, const double& i_xmax,
			     int i_ydeg, const double& i_ymin, const double& i_ymax
			     ) : 
  xdeg(i_xdeg),
  xmin(i_xmin),
  xmax(i_xmax),
  ydeg(i_ydeg),
  ymin(i_ymin),
  ymax(i_ymax)
{
  coeff.resize((xdeg+1)*(ydeg+1));
  coeff *= 0;
}
 
harp::vector_double specex::Legendre2DPol::Monomials(const double &x, const double &y) const {
  
  // range is -1,1 if  xmin<x<xmax
  double rx= 2*(x-xmin)/(xmax-xmin)-1;
  double ry= 2*(y-ymin)/(ymax-ymin)-1;
  
  harp::vector_double mx(xdeg+1);
  for(int i=0;i<=xdeg;i++)
    mx(i)=LegendrePol(i,rx);
  
  harp::vector_double my(ydeg+1);
  for(int j=0;j<=ydeg;j++)
    my(j)=LegendrePol(j,ry);
  
  harp::vector_double m((xdeg+1)*(ydeg+1));
  for(int j=0;j<=ydeg;j++) {
    for(int i=0;i<=xdeg;i++) {
      m(i+j*(xdeg+1))=mx(i)*my(j);
    }
  }
  return m;
}


double specex::Legendre2DPol::Value(const double &x,const double &y) const {
  return specex::dot(coeff,Monomials(x,y));
}

