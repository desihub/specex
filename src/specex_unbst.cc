#include "specex_unbst.h"

// vout = vin
void specex::unbst::subcopy(const unhrp::vector_double& vin, int in0, int in1,
			     unhrp::vector_double& vout, int out0){
  for(int i=0; i<in1-in0; i++) vout[i+out0] = vin[i+in0];
}
void specex::unbst::subcopy(const unhrp::vector_double& vin, 
			     unhrp::vector_double& vout, int out0){
  for(int i=0; i<vin.size(); i++) vout[i+out0] = vin[i];
}
void specex::unbst::subcopy(const unhrp::vector_double& vin, 
			     unhrp::vector_double& vout, int out0, double alpha){
  for(int i=0; i<vin.size(); i++) vout[i+out0] = alpha*vin[i];
}

// vout += vin
void specex::unbst::subadd(const unhrp::vector_double& vin, int in0, int in1, unhrp::vector_double& vout, int out0, double alpha){
  for(int i=0; i<in1-in0; i++) vout[i+out0] += alpha*vin[i+in0];
}
void specex::unbst::subadd(const unhrp::vector_double& vin, unhrp::vector_double& vout, int out0, double alpha){
  for(int i=0; i<vin.size(); i++) vout[i+out0] += alpha*vin[i];
}

unhrp::vector_double specex::unbst::scalevec(const unhrp::vector_double& vin, double alpha){
  unhrp::vector_double vout(vin.size());
  for(int i=0; i<vin.size(); i++) vout[i]=alpha*vin[i];
  return vout;
}

unhrp::vector_double specex::unbst::subrange(const unhrp::vector_double& vin, int in0, int in1){
  unhrp::vector_double vout(in1-in0);
  for(int i=0; i<in1-in0; i++) vout[i]=vin[i+in0];
  return vout;
}


