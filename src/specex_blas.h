#ifndef SPECEX_BLAS_H_
#define SPECEX_BLAS_H_

extern "C" {
  
  double specex_dot(int, const double *, const double *);
  
  void specex_axpy(int, const double *, const double *, const double *);
  void specex_syr(int, const double *, const double *, const double *);
  void specex_syrk(int, int, const double *, const double *, const double *,
		   const double *);  
  void specex_gemv(int, int, const double *, const double *, const double *,
		   const double *, const double*);
  void specex_gemm(int, int, int, const double *, const double *, const double *,
		   const double *, const double*);
  
}

#endif
