#ifndef SPECEX_UNHRP__H
#define SPECEX_UNHRP__H

#include <boost/numeric/bindings/ublas.hpp>

namespace unbls {
  
  typedef boost::numeric::ublas::vector < int > vector_int;
  typedef boost::numeric::ublas::vector < double > vector_double;
  typedef boost::numeric::ublas::matrix < double, boost::numeric::ublas::column_major > matrix_double;
  typedef boost::numeric::ublas::vector < uint8_t > vector_mask;
  typedef boost::numeric::ublas::matrix < uint8_t, boost::numeric::ublas::column_major > matrix_mask;
  
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
