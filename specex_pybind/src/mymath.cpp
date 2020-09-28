// Licensed under a 3-clause BSD style license - see LICENSE.rst

#include <mymath.h>
#include <ext.h>

namespace spx = specex;

spx::MyMath::MyMath(int id){
  
  this->id = id;
  
}

int spx::MyMath::add(int a, int b){

  std::cout << "id in add = " << this->id << " " << std::endl;

  return ext_add(a, b);
}
