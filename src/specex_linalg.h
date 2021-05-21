#ifndef SPECEX_LINALG__H
#define SPECEX_LINALG__H

#include <unhrp.h>

namespace specex {
    
  // ! scalar product
  double dot(const unhrp::vector_double&, const unhrp::vector_double&);
  double dot(const unhrp::vector_double&, int, int, const unhrp::vector_double&);
  
  // ! B += a*h
  void axpy(const double&, const unhrp::vector_double&, unhrp::vector_double&);
  void axpy(const double&, const unhrp::vector_double&, int, int, unhrp::vector_double&);
  
  // !  A += w*h*h.transposed(), where A is a symmetric matrx (only lower half is filled!)
  void syr(const double&, const unhrp::vector_double&, unhrp::matrix_double&);
  void syr(const double&, const unhrp::vector_double&, int, int, unhrp::matrix_double&);

  // ! B := alpha*A*h + beta*B (B += alpha*A*h for beta=1)
  void gemv(const double &alpha,  const unhrp::matrix_double &A,  const unhrp::vector_double& h, const double &beta, unhrp::vector_double& B);

  // ! C = alpha*A*B + beta*C if side='L' , C = alpha*B*A + beta*C if side='R', where A is a symmetric matrix
  // void symm(const char side, const double& alpha, const unhrp::matrix_double &A, const unhrp::matrix_double &B, const double& beta, unhrp::matrix_double &C);
  
  // ! C = alpha*A*B + beta*C if side='L' , C = alpha*B*A + beta*C if side='R', where A is a symmetric matrix
  void gemm(const double& alpha, const unhrp::matrix_double &A, const unhrp::matrix_double &B, const double& beta, unhrp::matrix_double &C);
  
  // ! C = alpha*A*At + beta*C
  void syrk(const double& alpha, const unhrp::matrix_double &A, const double& beta, unhrp::matrix_double &C);

  // must exist somewhere in boost
  void minmax(const unhrp::vector_double& v, double& minv, double& maxv);
  
  // !x*x
  double square(const double& x);
  
  int cholesky_solve(unhrp::matrix_double& A, unhrp::vector_double& B);
  
  // ! assumes A has been through cholesky_solve before
  int cholesky_invert_after_decomposition(unhrp::matrix_double& A);
  
}
#endif
