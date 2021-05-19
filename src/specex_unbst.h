#ifndef SPECEX_UNBST__H
#define SPECEX_UNBST__H
#include <unbls.h>
#include <vector>
#include <type_traits>

// unbst: functions to replace basic boost vector operations

namespace specex::unbst {

  void subcopy(const unbls::vector_double&, int, int, unbls::vector_double&, int);
  void subcopy(const unbls::vector_double&, unbls::vector_double&, int);
  void subcopy(const unbls::vector_double&, unbls::vector_double&, int, double);
  
  void subadd(const  unbls::vector_double&, int, int, unbls::vector_double&, int, double);
  void subadd(const  unbls::vector_double&, unbls::vector_double&, int, double);
  
  template<class vtype1, class vtype2>
  static void subadd(vtype1 &vin, vtype2 &vout, int i0){
    for(int i=0; i<vin.size(); i++) vout[i+i0] += vin[i];
  }
  
  unbls::vector_double scalevec(const unbls::vector_double&, double);
  unbls::vector_double subrange(const unbls::vector_double&, int, int);

}

#endif
