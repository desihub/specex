#ifndef SPECEX_UNHRP__H
#define SPECEX_UNHRP__H

#include <boost/numeric/bindings/ublas.hpp>

namespace unbls {

  //typedef boost::numeric::ublas::vector < int >     vector_int;
  //typedef boost::numeric::ublas::vector < double >  vector_double;
  //typedef boost::numeric::ublas::vector < uint8_t > vector_mask;
  
  typedef std::vector < int >     vector_int;
  typedef std::vector < double >  vector_double;
  typedef std::vector < uint8_t > vector_mask;
  
  typedef boost::numeric::ublas::matrix < double, boost::numeric::ublas::column_major > matrix_double;
  typedef boost::numeric::ublas::matrix < uint8_t, boost::numeric::ublas::column_major > matrix_mask;

  inline void zero(vector_int &v){
    std::fill(v.begin(),v.end(),0);
  }

  inline void zero(vector_double &v){
    std::fill(v.begin(),v.end(),0.0);
  }

  inline void zero(matrix_double &m){

    for (unsigned i = 0; i < m.size1 (); ++ i)
      for (unsigned j = 0; j < m.size2 (); ++ j)
	m(i,j) = 0.;
  }

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
