#include <string>
#include <iostream>
#include <specex_psfpy.h>
using namespace std;

void specex::PSFPy::SetCoeff2d(specex::image_data coeff2d, bool is_x){
  if(is_x){
    this->coeff2d_x = coeff2d;
  } else{
    this->coeff2d_y = coeff2d;
  }
}

specex::image_data specex::PSFPy::GetCoeff2d(bool is_x){
  if(is_x){
    return coeff2d_x;
  }else{
    return coeff2d_y;
  }
}
