#ifndef SPECEX_UNBST__H
#define SPECEX_UNBST__H
#include <vector>

// unbst: functions to replace basic boost/numeric/ublas vector operations

namespace specex::unbst {

  // returns vin[i0:i1]
  template<class vtype>
  static vtype subrange(vtype vin, int i0, int i1){
    int n = i1-i0;
    vtype vout(n);
    std::copy(vin.begin()+i0,vin.begin()+i1,vout.begin());
    return vout;
  }
    
  // returns alpha*vin[i0:i1]
  template<class vtype1, class vtype2>
    static vtype1 subrange(vtype1 vin, int i0, int i1, vtype2 alpha){
    int n = i1-i0;
    vtype1 vout(n);
    std::copy(vin.begin()+i0,vin.begin()+i1,vout.begin());
    for(int i=0; i<vout.size();i++) vout[i] *= alpha;
    return vout;
  }
    
  // assigns vout[i0:] = vin
  template<class vtype1, class vtype2>
  static void subcopy(vtype1 vin, vtype2 &vout, int i0){
    std::copy(vin.begin(),vin.end(),vout.begin()+i0);
  }

  // assigns vout[i0:] += vin
  template<class vtype1, class vtype2>
  static void subadd(vtype1 vin, vtype2 &vout, int i0){
    for(int i=0; i<vin.size(); i++) vout[i+i0] += vin[i];
  }

  // assigns vout[i0:] *= alpha*vin
  template<class vtype1, class vtype2, class vtype3>
  static void subadd(vtype1 vin, vtype2 &vout, int i0, vtype3 alpha){
    for(int i=0; i<vin.size(); i++) vout[i+i0] += alpha*vin[i];
  }
  
}

#endif