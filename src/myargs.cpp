// Licensed under a 3-clause BSD style license - see LICENSE.rst

#include <myargs.h>

#include <sstream>

namespace spx = specex;

spx::MyArgs::MyArgs(int id){
  
  this->ids = id;

  return;
}

void spx::MyArgs::get_args(int argc, char **argv){
  std::vector<std::string> argvc(argv, argv + argc);
  this->argc = argc;
  this->argv = argvc;

  return;
}

void spx::MyArgs::print_args(){
  std::cout << "id in print_args = " << this->ids << std::endl;
  for (auto iarg = 0; iarg < this->argc; iarg++)
    std::cout << "arg " << iarg << " is: "
	      << this->argv[iarg] << std::endl;
    
  return;
}

