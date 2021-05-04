#ifndef SPECEX_UNHRP__H
#define SPECEX_UNHRP__H

#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/storage.hpp>

#include <boost/numeric/bindings/ublas.hpp>
#include <boost/numeric/bindings/blas.hpp>
#include <boost/numeric/bindings/lapack.hpp>
#include <boost/numeric/bindings/views.hpp>

// specex::unhrp
// definitions, methods and classes necessary for SPECEX calculations that were originally
// part of HARP but have been modified to
//   - remove dependencies on boost
//   - be a part of specex source
// see https://github.com/tskisner/HARP for the original HARP package developed by
// Ted Kisner

namespace specex::unhrp {

  typedef enum {
    EIG_NONE,
    EIG_SQRT,
    EIG_INVSQRT,
    EIG_INV
  } eigen_op;

  typedef boost::numeric::ublas::matrix < double, boost::numeric::ublas::column_major > matrix_double;

  typedef boost::numeric::ublas::compressed_matrix < double, boost::numeric::ublas::row_major > matrix_double_sparse;

  typedef boost::numeric::ublas::vector < double > vector_double;

  typedef boost::numeric::ublas::matrix < uint8_t, boost::numeric::ublas::column_major > matrix_mask;

  typedef boost::numeric::ublas::vector < uint8_t > vector_mask;

  typedef boost::numeric::ublas::matrix < float, boost::numeric::ublas::column_major > matrix_float;

  typedef boost::numeric::ublas::compressed_matrix < float, boost::numeric::ublas::row_major > matrix_float_sparse;

  typedef boost::numeric::ublas::vector < float > vector_float;


}

#endif
