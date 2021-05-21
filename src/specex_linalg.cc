#include <specex_linalg.h>
#include <specex_message.h>
#include <specex_blas.h>
#include <specex_lapack.h>
#include <vector>

// contains all calls to C-wrappers (specex_*) calling C-interface BLAS and LAPACK functions 

// returns the dot product of x and y
double specex::dot(const unhrp::vector_double& x, const unhrp::vector_double& y) {  
  return specex_dot(x.size(), &x[0],  &y[0]);  
}
double specex::dot(const unhrp::vector_double& x, int i0, int i1, const unhrp::vector_double& y) {  
  return specex_dot(i1-i0,    &x[i0], &y[0]);  
}
  
// y += alpha*x
void specex::axpy(const double &alpha, const unhrp::vector_double& x,
		  unhrp::vector_double& y) {  
  specex_axpy(x.size(), &alpha, &x[0],  &y[0]);
}
void specex::axpy(const double &alpha, const unhrp::vector_double& x, int i0, int i1,
		  unhrp::vector_double& y) {  
  specex_axpy(i1-i0,    &alpha, &x[i0], &y[0]);
}

// A += alpha*x*x**T, where A is a symmetric matrx (only lower half is filled)
void specex::syr(const double& alpha, const unhrp::vector_double& x,
		 unhrp::matrix_double& A) {  
  specex_syr(x.size(), &alpha, &x[0],  &A(0,0)); 
}
void specex::syr(const double& alpha, const unhrp::vector_double& x, int i0, int i1,
		 unhrp::matrix_double& A) {  
  specex_syr(i1-i0,    &alpha, &x[i0], &A(0,0)); 
}

// C = alpha*A*A**T + beta*C
void specex::syrk(const double& alpha, const unhrp::matrix_double &A, const double& beta, unhrp::matrix_double &C) {
  specex_syrk(A.size1(), A.size2(), &alpha, &A(0,0), &beta, &C(0,0));
}

// y = alpha*A*x + beta*y 
void specex::gemv(const double &alpha,  const unhrp::matrix_double &A,  const unhrp::vector_double& x, const double &beta, unhrp::vector_double& y) {  
  specex_gemv(A.size1(),A.size2(), &alpha, &A(0,0), &x[0], &beta, &y[0]);  
}

// C = alpha*A*B + beta*C
void specex::gemm(const double& alpha, const unhrp::matrix_double &A, const unhrp::matrix_double &B, const double& beta, unhrp::matrix_double &C) {
  specex_gemm(A.size1(), B.size2(), A.size2(), &alpha, &A(0,0), &B(0,0), &beta, &C(0,0));  
}

// returns the solution, x, to a real system of linear equations
//   A * x = b,
// solution is returned in b, i.e. b --> x, for return value 0
int specex::cholesky_solve(unhrp::matrix_double& A, unhrp::vector_double& b) {
  return specex_posv(b.size(),&A(0,0),&b[0]);
}

// invert matrix A in place; A := inv(A)
int specex::cholesky_invert_after_decomposition(unhrp::matrix_double& A) {
  return specex_potri(A.size1(),&A(0,0));
}

// min and max of vector
void specex::minmax(const unhrp::vector_double& v, double& minv, double& maxv) {
  unhrp::vector_double::const_iterator it=v.begin();
  minv=maxv=(*it);
  for(; it!=v.end() ; ++it) {
    if(*it<minv) minv=*it;
    if(*it>maxv) maxv=*it;
  }
}

// square of scalar
double specex::square(const double& x) {
  return x*x;
}

