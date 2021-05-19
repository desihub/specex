#ifndef SPECEX_UNHRP__H
#define SPECEX_UNHRP__H

/*
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/bindings/blas.hpp>
#include <boost/numeric/bindings/lapack.hpp>
#include <boost/numeric/bindings/views.hpp>
*/

#include <harp.hpp>
#include <boost/numeric/bindings/ublas.hpp>

// unhrp:
// definitions, methods and classes necessary for SPECEX calculations that were originally
// part of HARP but have been modified to
//   - remove dependencies on boost
//   - be a part of specex source
// see https://github.com/tskisner/HARP for the original HARP package developed by
// Ted Kisner

namespace unhrp {
  
  typedef boost::numeric::ublas::vector < int > vector_int;
  typedef boost::numeric::ublas::vector < double > vector_double;
  typedef boost::numeric::ublas::matrix < double, boost::numeric::ublas::column_major > matrix_double;

}

template < class T >
inline std::ostream& operator << (std::ostream& os, const std::vector<T>& v) 
{
    os << "[";
    for (int i = 0; i<v.size(); i++)
    {
        os << " " << v[i];
    }
    os << " ]";
    return os;
}

#endif
