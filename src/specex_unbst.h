#ifndef SPECEX_UNBST__H
#define SPECEX_UNBST__H
#include <unhrp.h>
#include <vector>
#include <type_traits>

// unbst: functions to replace basic boost vector operations

namespace specex::unbst {

  void subcopy(const unhrp::vector_double&, int, int, unhrp::vector_double&, int);
  void subcopy(const unhrp::vector_double&, unhrp::vector_double&, int);
  void subcopy(const unhrp::vector_double&, unhrp::vector_double&, int, double);
  
  void subadd(const  unhrp::vector_double&, int, int, unhrp::vector_double&, int, double);
  void subadd(const  unhrp::vector_double&, unhrp::vector_double&, int, double);
  
  template<class vtype1, class vtype2>
  static void subadd(vtype1 &vin, vtype2 &vout, int i0){
    for(int i=0; i<vin.size(); i++) vout[i+i0] += vin[i];
  }
  
  unhrp::vector_double scalevec(const unhrp::vector_double&, double);
  unhrp::vector_double subrange(const unhrp::vector_double&, int, int);

}

#endif
