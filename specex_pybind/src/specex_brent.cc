#include <iostream>
#include <cstdlib> // for abort

#include "specex_brent.h"

#define CGOLD 0.3819660
#define ZEPS 1.e-60
/*
Here ITMAX is the maximum allowed number of iterations; CGOLD is the golden ratio; ZEPS is
a small number that protects against trying to achieve fractional accuracy for a minimum that
happens to be exactly zero.
*/

#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

double sign(const double& a, const double& b) {
  if(b>0) return fabs(a);
  return -fabs(a);
}

double brent(AnalyticFunction * f, double ax, double bx, double cx, double tol, const void *context, double& fval, int& status, int ITMAX)
	   
/* Given a function f, and given a bracketing triplet of abscissas ax,
   bx, cx (such that bx is between ax and cx, and f(bx) is less than
   both f(ax) and f(cx)), this routine isolates the minimum to a
   fractional precision of about tol using Brent method. The abscissa
   of the minimum is returned as the function value (brent), and the
   minimum function value is returned as fval.
*/
{
  status=0;
  int iter;
  double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  d=0.;
  double e=0.0;  // This will be the distance moved on the step before last.
  a=(ax < cx ? ax : cx); //a and b must be in ascending order,

  b=(ax > cx ? ax : cx); //but input abscissas need  not be.
  x=w=v=bx; // Initializations...
  fw=fv=fx=(*f)(x,context);
  for (iter=1;iter<=ITMAX;iter++) { // Main program loop.
    xm=0.5*(a+b);
    tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) { // Test for done here.
      fval = fx;
      return x;
    }
    if (fabs(e) > tol1) { // Construct a trial parabolic fit.
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if (q > 0.0) p = -p;
      q=fabs(q);
      etemp=e;
      e=d;
      if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	d=CGOLD*(e=(x >= xm ? a-x : b-x));
      /* The above conditions determine the acceptability of the
	 parabolic fit. Here we take the golden section step into the
	 larger of the two segments.
      */ 
      else {
	d=p/q; //Take the parabolic step.
	u=x+d;
	if (u-a < tol2 || b-u < tol2)
	  d=sign(tol1,xm-x);
      }
    } else {
      d=CGOLD*(e=(x >= xm ? a-x : b-x));
    }
    u=(fabs(d) >= tol1 ? x+d : x+sign(tol1,d));
    fu=(*f)(u,context);
    /*This is the one function evaluation per iteration.*/
    if (fu <= fx) { // Now decide what to do with our function evaluatiON
      if (u >= x) a=x; else b=x; 
      
      SHFT(v,w,x,u); // Housekeeping follows:
      SHFT(fv,fw,fx,fu);
    } else {
      if (u < x) a=u; else b=u;
      if (fu <= fw || w == x) {
	v=w;
	w=u;
	fv=fw;
	fw=fu;
      } else if (fu <= fv || v == x || v == w) {
	v=u;
	fv=fu;
      }
    } //Done with housekeeping. Back for
  } //another iteration.
  std::cerr << "WARNING Too many iterations in brent" << std::endl;
  status=1;
  fval = fx;
  return x; // never get here
}


double extended_brent(AnalyticFunction * f, double ax, double bx, double cx, double tol, const void *context, double& fval, int& status, int ITMAX) {

  if( ! (ax<bx && bx<cx ) ) {
    std::cerr << "need ax<bx<cx" << std::endl;
    abort();
  }
  
  // doesnt need chi2(bx)<chi2(ax) &&  chi2(bx)<chi2(cx)
  double fa,fb,fc;
  fa = f(ax,context);
  fb = f(bx,context);
  fc = f(cx,context);
  
  if(fb>fa && fb>fc) {
    std::cout << "extended_brent : f has a maximum, better exit" << std::endl;
    status=1;
    fval=f(0,context);
    return 0;
  }

  double x_min = bx;
  double f_min = fb; // keep track of best value if failure to bound minimum
  
  if(fa<f_min) {f_min = fa; x_min=ax;}
  if(fb<f_min) {f_min = fb; x_min=bx;}
  
  
  int count = 0;
  while(count<ITMAX) {
    if(0) {
      std::cout << "eb " << count
	   << " f(" << ax << ")=" << fa
	   << " f(" << bx << ")=" << fb
	   << " f(" << cx << ")=" << fc
	   << std::endl;
    }
    if(fa>fb && fb<fc)
      return brent(f,ax,bx,cx,tol,context,fval,status,ITMAX);
    if(fc<fb) {
      fb = fc;
      bx = cx;
      cx += (cx-ax);
      fc = f(cx,context);
      if(fc<f_min) {f_min = fc; x_min=cx;}
    }else{
      fb = fa;
      bx = ax;
      ax -= (cx-ax);
      fa = f(ax,context);
      if(fa<f_min) {f_min = fa; x_min=ax;}
    }
    count ++;
  }
  std::cout << "extended_brent : failed to bound minimum, x_min = " << x_min << std::endl;
  fval = f_min;
  status = 0;
  return x_min;
}


