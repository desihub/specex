#ifdef USE_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif

// contains all calls to C-interface LAPACK functions (LAPACKE)
// https://www.netlib.org/lapack/lapacke.html

// returns the solution, x, to a real system of linear equations
//   A * x = b,
// solution is returned in b, i.e. b --> x, for return value 0
int specex_posv(int n, const double *A, const double *b){
  
  return LAPACKE_dposv(LAPACK_COL_MAJOR, 'L', n, 1, A, n, b, n);
  
}

// invert matrix A in place; A := inv(A)
// http://www.netlib.org/lapack/explore-html/d8/d63/dpotri_8f_source.html
int specex_potri(int n, const double *A){
  
  return LAPACKE_dpotri(LAPACK_COL_MAJOR, 'L', n, A, n);
  
}


