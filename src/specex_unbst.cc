#include "specex_unbst.h"

// vout = vin
void specex::unbst::subcopy(const harp::vector_double& vin, int in0, int in1,
			     harp::vector_double& vout, int out0){
  for(int i=0; i<in1-in0; i++) vout[i+out0] = vin[i+in0];
}
void specex::unbst::subcopy(const harp::vector_double& vin, 
			     harp::vector_double& vout, int out0){
  for(int i=0; i<vin.size(); i++) vout[i+out0] = vin[i];
}
void specex::unbst::subcopy(const harp::vector_double& vin, 
			     harp::vector_double& vout, int out0, double alpha){
  for(int i=0; i<vin.size(); i++) vout[i+out0] = alpha*vin[i];
}

// vout += vin
void specex::unbst::subadd(const harp::vector_double& vin, int in0, int in1, harp::vector_double& vout, int out0, double alpha){
  for(int i=0; i<in1-in0; i++) vout[i+out0] += alpha*vin[i+in0];
}
void specex::unbst::subadd(const harp::vector_double& vin, harp::vector_double& vout, int out0, double alpha){
  for(int i=0; i<vin.size(); i++) vout[i+out0] += alpha*vin[i];
}

harp::vector_double specex::unbst::scalevec(const harp::vector_double& vin, double alpha){
  harp::vector_double vout(vin.size());
  for(int i=0; i<vin.size(); i++) vout[i]=alpha*vin[i];
  return vout;
}

harp::vector_double specex::unbst::subrange(const harp::vector_double& vin, int in0, int in1){
  harp::vector_double vout(in1-in0);
  for(int i=0; i<in1-in0; i++) vout[i]=vin[i+in0];
  return vout;
}


