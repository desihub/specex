#ifndef SPECEX_LINALG__H
#define SPECEX_LINALG__H

#include <unbls.h>

namespace specex {
    
  // ! scalar product
  double dot(const unbls::vector_double&, const unbls::vector_double&);
  double dot(const unbls::vector_double&, int, int, const unbls::vector_double&);
  
  // ! B += a*h
  void axpy(const double&, const unbls::vector_double&, unbls::vector_double&);
  void axpy(const double&, const unbls::vector_double&, int, int, unbls::vector_double&);
  
  // !  A += w*h*h.transposed(), where A is a symmetric matrx (only lower half is filled!)
  void syr(const double&, const unbls::vector_double&, unbls::matrix_double&);
  void syr(const double&, const unbls::vector_double&, int, int, unbls::matrix_double&);

  // ! B := alpha*A*h + beta*B (B += alpha*A*h for beta=1)
  void gemv(const double &alpha,  const unbls::matrix_double &A,  const unbls::vector_double& h, const double &beta, unbls::vector_double& B);

  // ! C = alpha*A*B + beta*C if side='L' , C = alpha*B*A + beta*C if side='R', where A is a symmetric matrix
  // void symm(const char side, const double& alpha, const unbls::matrix_double &A, const unbls::matrix_double &B, const double& beta, unbls::matrix_double &C);
  
  // ! C = alpha*A*B + beta*C if side='L' , C = alpha*B*A + beta*C if side='R', where A is a symmetric matrix
  void gemm(const double& alpha, const unbls::matrix_double &A, const unbls::matrix_double &B, const double& beta, unbls::matrix_double &C);
  
  // ! C = alpha*A*At + beta*C
  void syrk(const double& alpha, const unbls::matrix_double &A, const double& beta, unbls::matrix_double &C);

  // must exist somewhere in boost
  void minmax(const unbls::vector_double& v, double& minv, double& maxv);
  
  // !x*x
  double square(const double& x);
  
  int cholesky_solve(unbls::matrix_double& A, unbls::vector_double& B);
  
  // ! assumes A has been through cholesky_solve before
  int cholesky_invert_after_decomposition(unbls::matrix_double& A);
  
}
#endif
