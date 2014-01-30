#include <specex_linalg.h>


//#define CHECK_BOUNDS

double specex::dot(const harp::vector_double& v1, const harp::vector_double& v2) {
#ifdef CHECK_BOUNDS
  if(v1.size() != v2.size())
    HARP_THROW("not same size");
#endif
  return blas::dot(v1,v2);
}
  



// !  A += w*h*h.transposed(), where A is a symmetric matrx (only lower half is filled!)
// see http://svn.boost.org/svn/boost/sandbox/numeric_bindings/libs/numeric/bindings/doc/html/boost_numeric_bindings/reference/blas/level_2_blas/syr.html
void specex::syr(const double& w, const harp::vector_double& h, harp::matrix_double& A) {
  blas::syr(w,h,boost::numeric::bindings::lower(A));
}

// ! B += a*h;
// see http://svn.boost.org/svn/boost/sandbox/numeric_bindings/libs/numeric/bindings/doc/html/boost_numeric_bindings/reference/blas/level_1_blas/axpy.html
void specex::axpy(const double &a, const harp::vector_double& h,  harp::vector_double& B) {
  blas::axpy(a,h,B);
}
// ! B := alpha*A*h + beta*B (B += alpha*A*h for beta=1)
void specex::gemv(const double &alpha,  const harp::matrix_double &A,  const harp::vector_double& h, const double &beta, harp::vector_double& B) {
  blas::gemv(alpha,A,h,beta,B);
}
// ! C = alpha*A*B + beta*C if side='L' , C = alpha*B*A + beta*C if side='R' , where A is a symmetric matrix
// see vn.boost.org/svn/boost/sandbox/numeric_bindings/libs/numeric/bindings/doc/html/boost_numeric_bindings/reference/blas/level_3_blas/symm.html
void specex::symm(const char side, const double& alpha, const harp::matrix_double &A, const harp::matrix_double &B, const double& beta, harp::matrix_double &C) {
#warning NEED TO FIX THIS !
  HARP_THROW("specex::symm not implemented");
  //blas::symm('L',alpha,boost::numeric::bindings::lower(A),B,beta,C);
}

// ! C = alpha*A*B + beta*C 
// see vn.boost.org/svn/boost/sandbox/numeric_bindings/libs/numeric/bindings/doc/html/boost_numeric_bindings/reference/blas/level_3_blas/gemm.html
void specex::gemm(const double& alpha, const harp::matrix_double &A, const harp::matrix_double &B, const double& beta, harp::matrix_double &C) {
  blas::gemm(alpha,A,B,beta,C);
}
// ! C = alpha*A*At + beta*C
void specex::syrk(const double& alpha, const harp::matrix_double &A, const double& beta, harp::matrix_double &C) {
  blas::syrk(alpha,A,beta,boost::numeric::bindings::lower(C));
}
// must exist somewhere in boost
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

int specex::cholesky_solve(harp::matrix_double& A, harp::vector_double& B) { 
  return lapack::posv(boost::numeric::bindings::lower(A),B);
}

// ! assumes A has been through cholesky_solve before
int specex::cholesky_invert_after_decomposition(harp::matrix_double& A) {
  return lapack::potri(boost::numeric::bindings::lower(A));
}
