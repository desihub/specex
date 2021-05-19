#include "specex_unbst.h"

// vout = vin
void specex::unbst::subcopy(const unbls::vector_double& vin, int in0, int in1,
			     unbls::vector_double& vout, int out0){
  for(int i=0; i<in1-in0; i++) vout[i+out0] = vin[i+in0];
}
void specex::unbst::subcopy(const unbls::vector_double& vin, 
			     unbls::vector_double& vout, int out0){
  for(int i=0; i<vin.size(); i++) vout[i+out0] = vin[i];
}
void specex::unbst::subcopy(const unbls::vector_double& vin, 
			     unbls::vector_double& vout, int out0, double alpha){
  for(int i=0; i<vin.size(); i++) vout[i+out0] = alpha*vin[i];
}

// vout += vin
void specex::unbst::subadd(const unbls::vector_double& vin, int in0, int in1, unbls::vector_double& vout, int out0, double alpha){
  for(int i=0; i<in1-in0; i++) vout[i+out0] += alpha*vin[i+in0];
}
void specex::unbst::subadd(const unbls::vector_double& vin, unbls::vector_double& vout, int out0, double alpha){
  for(int i=0; i<vin.size(); i++) vout[i+out0] += alpha*vin[i];
}

unbls::vector_double specex::unbst::scalevec(const unbls::vector_double& vin, double alpha){
  unbls::vector_double vout(vin.size());
  for(int i=0; i<vin.size(); i++) vout[i]=alpha*vin[i];
  return vout;
}

unbls::vector_double specex::unbst::subrange(const unbls::vector_double& vin, int in0, int in1){
  unbls::vector_double vout(in1-in0);
  for(int i=0; i<in1-in0; i++) vout[i]=vin[i+in0];
  return vout;
}


