#ifndef SPECEX_UNHRP__H
#define SPECEX_UNHRP__H

#include <vector>
#include <cstdint>
#include <ostream>
#include <memory>
#include <map>

namespace unbls {

  typedef std::vector < int >     vector_int;
  typedef std::vector < double >  vector_double;
  typedef std::vector < uint8_t > vector_mask;

  // matrix class
  class matrix {
    
  protected :
    size_t _nrows;
    size_t _ncols;
    
  public :

    matrix ();
    
    size_t size1() { return _nrows; }
    size_t size2() { return _ncols; }
      
    const size_t size1() const { return _nrows; }
    const size_t size2() const { return _ncols; }
      
  };
  
  // matrix_double class
  class matrix_double : public matrix {

  public :
    
    std::vector < double >  vals;
    
    matrix_double ( );
    matrix_double ( size_t nrows, size_t ncols);
    matrix_double ( size_t nrows, size_t ncols, const unbls::vector_double& i_data);
    
    void resize( size_t nrows, size_t ncols);
    double* address();
    
    double& operator()(const int i, const int j) { return vals[i+j*_nrows]; }    
    const double& operator()(const int i, const int j) const { return vals[i+j*_nrows]; }
      
  };
  
  // matrix_mask class
  class matrix_mask : public matrix {

  public :

    std::vector < uint8_t > vals;

    matrix_mask ( );
    matrix_mask ( size_t nrows, size_t ncols);
    matrix_mask ( size_t nrows, size_t ncols, const std::vector<uint8_t>& i_data);
    
    void resize( size_t nrows, size_t ncols);
    uint8_t* address();

    uint8_t& operator()(const int i, const int j) { return vals[i+j*_nrows]; }    
    const uint8_t& operator()(const int i, const int j) const { return vals[i+j*_nrows]; }
      
  };
  
  template <class T>
  void zero(std::vector<T>& v){
    std::fill(v.begin(),v.end(),0);
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
