#ifdef USE_MKL
#include <mkl.h>
#else
#include <cblas.h>
#endif

double specex_dot(int n, const double *x, const double *y){
  
  return cblas_ddot(n, x, 1, y, 1);

}

void specex_syr(int n, int m, const double *w, const double *h, const double *A){

  CBLAS_LAYOUT layout;
  CBLAS_UPLO   uplo;

  layout = CblasColMajor;
  uplo   = CblasLower;
  
  cblas_dsyr(layout, uplo, n, *w, h, 1, A, m);
  
}

void specex_axpy(int n, const double *a, const double *h, const double *B){

  cblas_daxpy(n, *a, h, 1, B, 1);
  
}

void specex_gemv(int n, int m, const double *alpha, const double *A, const double *h,
		 const double *beta, const double *B){

  CBLAS_LAYOUT    layout;
  CBLAS_TRANSPOSE transa;

  layout = CblasColMajor;
  transa = CblasNoTrans;

  cblas_dgemv(layout, transa, n, m, *alpha, A, m, h, 1, *beta, B, 1);
  
}
