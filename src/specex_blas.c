#ifdef USE_MKL
#include <mkl.h>
#endif

double specex_dot(int n, const double *x, const double *y){
  
  return cblas_ddot(n, x, 1, y, 1);

}

