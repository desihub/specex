#ifndef SPECEX_UNBST__H
#define SPECEX_UNBST__H
#include <harp.hpp>
#include <vector>
#include <type_traits>

// unbst: functions to replace basic boost/numeric/ublas vector operations

namespace ublas  = boost::numeric::ublas;

namespace specex::unbst {

  void subcopyr(const harp::vector_double&, int, int, harp::vector_double&, int);
  void subcopyr(const harp::vector_double&, harp::vector_double&, int);
  
  void subaddr(const  harp::vector_double&, int, int, harp::vector_double&, int, double);
  void subaddr(const  harp::vector_double&, harp::vector_double&, int, double);

  harp::vector_double scalevec(const harp::vector_double&, double);
  
  harp::vector_double subrangerr(const harp::vector_double&, int, int);

  // returns vin[i0:i1]
  template<typename vtype>
  static vtype subrange(vtype& vin, int i0, int i1){    
    //return ublas::project(vin,ublas::range(i0,i1));
    int n = i1-i0;
    vtype vout(n);
    //for(int i=0; i<n; i++) vout[i]=vin[i+i0];
    std::copy(vin.begin()+i0,vin.begin()+i1,vout.begin());	
    return vout;
  }

  static harp::vector_double subranger(harp::vector_double vin, int i0, int i1){    
    //return ublas::project(vin,ublas::range(i0,i1));
    int n = i1-i0;
    harp::vector_double vout(n);
    for(int i=0; i<n; i++) vout[i]=vin[i+i0];
    //std::copy(vin.begin()+i0,vin.begin()+i1,vout.begin());	
    return vout;
  }

  // returns alpha*vin[i0:i1]
  template<class vtype1, class vtype2>
  static vtype1 subrange(vtype1 &vin, int i0, int i1, vtype2 alpha){
    //return ublas::project(alpha*vin,ublas::range(i0,i1));
    int n = i1-i0;
    vtype1 vout(n);
    //std::copy(vin.begin()+i0,vin.begin()+i1,vout.begin());
    for(int i=0; i<vout.size();i++) vout[i] = alpha*vin[i+i0];
    return vout;
  }
    
  // assigns vout[i0:] = vin
  template<class vtype1, class vtype2>
  static void subcopyt(vtype1 &vin, vtype2 &vout, int i0){
    //size_t m_size = vin.size();
    //ublas::noalias(ublas::project(vout,ublas::range(i0,i0+m_size))) = vin;
    //std::copy(vin.begin(),vin.end(),vout.begin()+i0);
  }

  // assigns vout[i0:] = vin
  static void subcopy(boost::numeric::ublas::vector<double>& vin, harp::vector_double &vout, int i0){
    //size_t m_size = vin.size();
    //ublas::noalias(ublas::project(vout,ublas::range(i0,i0+m_size))) = vin;
    std::copy(vin.begin(),vin.end(),vout.begin()+i0);
  }

    // assigns vout[i0:] = alpha*vin
  template<class vtype1, class vtype2, class vtype3>
  static void subcopy(vtype1 &vin, vtype2 &vout, int i0, vtype3 alpha){
    size_t m_size = vin.size();
    ublas::noalias(ublas::project(vout,ublas::range(i0,i0+m_size))) = alpha*vin;
    //for(int i=0; i<vin.size(); i++) vin[i] *= alpha;
    //std::copy(vin.begin(),vin.end(),vout.begin()+i0);
  }

  // assigns vout[i0:] += vin
  template<class vtype1, class vtype2>
  static void subadd(vtype1 &vin, vtype2 &vout, int i0){
    for(int i=0; i<vin.size(); i++) vout[i+i0] += vin[i];
  }

  // assigns vout[i0:] *= alpha*vin
  template<class vtype1, class vtype2, class vtype3>
  static void subadd(vtype1 &vin, vtype2 &vout, int i0, vtype3 alpha){
    for(int i=0; i<vin.size(); i++) vout[i+i0] += alpha*vin[i];
  }
  
}

#endif
