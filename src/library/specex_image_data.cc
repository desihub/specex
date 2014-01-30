#include <iostream>
#include <cstdio>
#include <string>

#include <harp.hpp>

#include <specex_image_data.h>
#include <specex_linalg.h>

using namespace std;

specex::image_data::image_data() : 
  harp::image()
{
 rows_=0;
 cols_=0;
}

specex::image_data::image_data(size_t ncols, size_t nrows) : 
  harp::image()
{
  resize(ncols,nrows);
}

specex::image_data::image_data(size_t ncols, size_t nrows, const harp::vector_double& i_data) : 
  harp::image()
{
  resize(ncols,nrows);
  if(i_data.size() != ncols*nrows) 
    HARP_THROW("Invalid input vector size");
  data = i_data; // copy
}

void specex::image_data::resize(size_t ncols, size_t nrows) {
  rows_ = nrows;
  cols_ = ncols;
  data.resize(rows_*cols_);
  data.clear();
}




