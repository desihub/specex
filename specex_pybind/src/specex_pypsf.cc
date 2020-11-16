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

void specex::PyPSF::get_table(std::vector<std::string>         &table_string,
			      std::vector<std::vector<int>>    &table_int,
			      std::vector<std::vector<double>> &table_double){

  int status = 0;
  int nrows  = this->psf->pydata.table.data.size();
  int ncols  = this->psf->pydata.table.columns.size();

  cout << "nrows in get_tablestring = " << nrows << endl;
  cout << "ncols in get_tablestring = " << ncols << endl;

  
  
}
