#include <iostream>
#include <cstdio>
#include <string>

#include <unbls.h>

#include <specex_image_data.h>
#include <specex_linalg.h>

using namespace std;

specex::image_data::image_data() : 
  unhrp::image()
{
 rows_=0;
 cols_=0;
}

specex::image_data::image_data(size_t ncols, size_t nrows) : 
  unhrp::image()
{
  resize(ncols,nrows);
}

specex::image_data::image_data(size_t ncols, size_t nrows, const unbls::vector_double& i_data) : 
  unhrp::image()
{
  resize(ncols,nrows);
  data = i_data; // copy
}

void specex::image_data::resize(size_t ncols, size_t nrows) {
  rows_ = nrows;
  cols_ = ncols;
  data.resize(rows_*cols_);
  data.clear();
}




