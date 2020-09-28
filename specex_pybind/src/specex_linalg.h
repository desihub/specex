#ifndef SPECEX_LINALG__H
#define SPECEX_LINALG__H

#include <harp.hpp>

#include <boost/numeric/ublas/io.hpp>
#include <cmath>
#include <algorithm>

// in harp :
// typedef boost::numeric::ublas::matrix < double, boost::numeric::ublas::column_major > matrix_double;


namespace ublas  = boost::numeric::ublas;
namespace blas   = boost::numeric::bindings::blas;
namespace lapack = boost::numeric::bindings::lapack;

namespace specex {
  
  
  //typedef boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::upper, boost::numeric::ublas::column_major> symmetric_matrix_double;
  
  // ! scalar product
  double dot(const harp::vector_double& v1, const harp::vector_double& v2);
  
  // !  A += w*h*h.transposed(), where A is a symmetric matrx (only lower half is filled!)
  void syr(const double& w, const harp::vector_double& h, harp::matrix_double& A);

  // ! B += a*h
  void axpy(const double &a, const harp::vector_double& h,  harp::vector_double& B);
  
  // ! B := alpha*A*h + beta*B (B += alpha*A*h for beta=1)
  void gemv(const double &alpha,  const harp::matrix_double &A,  const harp::vector_double& h, const double &beta, harp::vector_double& B);

  // ! C = alpha*A*B + beta*C if side='L' , C = alpha*B*A + beta*C if side='R', where A is a symmetric matrix
  // void symm(const char side, const double& alpha, const harp::matrix_double &A, const harp::matrix_double &B, const double& beta, harp::matrix_double &C);
  
  // ! C = alpha*A*B + beta*C if side='L' , C = alpha*B*A + beta*C if side='R', where A is a symmetric matrix
  void gemm(const double& alpha, const harp::matrix_double &A, const harp::matrix_double &B, const double& beta, harp::matrix_double &C);
  
  // ! C = alpha*A*At + beta*C
  void syrk(const double& alpha, const harp::matrix_double &A, const double& beta, harp::matrix_double &C);

  // must exist somewhere in boost
  void minmax(const harp::vector_double& v, double& minv, double& maxv);
  
  // !x*x
  double square(const double& x);
  
  int cholesky_solve(harp::matrix_double& A, harp::vector_double& B);
  
  // ! assumes A has been through cholesky_solve before
  int cholesky_invert_after_decomposition(harp::matrix_double& A);
  
}
#endif
