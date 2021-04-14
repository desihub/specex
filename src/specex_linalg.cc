#include <specex_linalg.h>
#include <specex_message.h>
#include <specex_blas.h>

// contains all calls to C-wrappers (specex_*) calling C-interface BLAS and LAPACK functions 

// returns the dot product of x and y
double specex::dot(const harp::vector_double& x, const harp::vector_double& y) {  
  return specex_dot(x.size(), &x[0], &y[0]);  
}
  
// y += alpha*x
void specex::axpy(const double &alpha, const harp::vector_double& x,  harp::vector_double& y) {  
  specex_axpy(x.size(), &alpha, &x[0], &y[0]);
}

// A += alpha*x*x**T, where A is a symmetric matrx (only lower half is filled)
void specex::syr(const double& alpha, const harp::vector_double& x, harp::matrix_double& A) {  
  specex_syr(x.size(), &alpha, &x[0], &A(0,0)); 
}

// C = alpha*A*A**T + beta*C
void specex::syrk(const double& alpha, const harp::matrix_double &A, const double& beta, harp::matrix_double &C) {
  specex_syrk(A.size1(), A.size2(), &alpha, &A(0,0), &beta, &C(0,0));
}

// y = alpha*A*x + beta*y 
void specex::gemv(const double &alpha,  const harp::matrix_double &A,  const harp::vector_double& x, const double &beta, harp::vector_double& y) {  
  specex_gemv(A.size1(),A.size2(), &alpha, &A(0,0), &x(0), &beta, &y[0]);  
}

// C = alpha*A*B + beta*C
void specex::gemm(const double& alpha, const harp::matrix_double &A, const harp::matrix_double &B, const double& beta, harp::matrix_double &C) {
  specex_gemm(A.size1(), B.size2(), A.size2(), &alpha, &A(0,0), &B(0,0), &beta, &C(0,0));  
}

// returns 
int specex::cholesky_solve(harp::matrix_double& A, harp::vector_double& B) { 
  return lapack::posv(boost::numeric::bindings::lower(A),B);
}

// ! assumes A has been through cholesky_solve before
int specex::cholesky_invert_after_decomposition(harp::matrix_double& A) {
  return lapack::potri(boost::numeric::bindings::lower(A));
}

// min and max of vector
void specex::minmax(const harp::vector_double& v, double& minv, double& maxv) {
  harp::vector_double::const_iterator it=v.begin();
  minv=maxv=(*it);
  for(; it!=v.end() ; ++it) {
    if(*it<minv) minv=*it;
    if(*it>maxv) maxv=*it;
  }
}

double specex::square(const double& x) {
  return x*x;
}

