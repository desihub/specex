#ifndef SPECEX_UNBST__H
#define SPECEX_UNBST__H
#include <harp.hpp>
#include <vector>
#include <type_traits>

// unbst: functions to replace basic boost/numeric/ublas vector operations

namespace ublas  = boost::numeric::ublas;

namespace specex::unbst {

  void subcopy(const harp::vector_double&, int, int, harp::vector_double&, int);
  void subcopy(const harp::vector_double&, harp::vector_double&, int);
  void subcopy(const harp::vector_double&, harp::vector_double&, int, double);
  
  void subadd(const  harp::vector_double&, int, int, harp::vector_double&, int, double);
  void subadd(const  harp::vector_double&, harp::vector_double&, int, double);
  
  template<class vtype1, class vtype2>
  static void subadd(vtype1 &vin, vtype2 &vout, int i0){
    for(int i=0; i<vin.size(); i++) vout[i+i0] += vin[i];
  }
  
  harp::vector_double scalevec(const harp::vector_double&, double);
  
  harp::vector_double subrange(const harp::vector_double&, int, int);

}

#endif
