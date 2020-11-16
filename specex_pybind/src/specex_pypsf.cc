#include <string>
#include <iostream>
#include <specex_pypsf.h>
using namespace std;

specex::image_data specex::PyPSF::get_trace(std::string axis){

  ///specex::PSF_p psf = this->psf;
  
  if(axis == "x"){
    return this->psf->pydata.coeff2d_x;
  }
  else{
    return this->psf->pydata.coeff2d_y;
  }
  
}

std::vector<std::string> specex::PyPSF::get_tablestring(){

  int status = 0;
  int nrows  = this->psf->pydata.table.data.size();
  int ncols  = this->psf->pydata.table.columns.size();

  cout << "nrows in get_tablestring = " << nrows << endl;
  cout << "ncols in get_tablestring = " << ncols << endl;

  std::vector<std::string> dummy;
  return dummy; 
  
}

std::vector<std::vector<int>> specex::PyPSF::get_tableint(){

  std::vector<std::vector<int>> dummy;
  return dummy;
  
}

std::vector<std::vector<double>> specex::PyPSF::get_tabledouble(){

  std::vector<std::vector<double>> dummy;
  return dummy; 
  
}

