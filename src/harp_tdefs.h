#ifndef HARP_TDEFS
#define HARP_TDEFS

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

namespace harp {

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
