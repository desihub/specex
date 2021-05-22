#include <iostream>
#include <cstdio>
#include <string>

#include <specex_unbls.h>

#include <specex_image_data.h>
#include <specex_linalg.h>

using namespace std;

// matrix methods
unbls::matrix::matrix() 
{
 _nrows=0;
 _ncols=0;
}

// matrix_double methods
unbls::matrix_double::matrix_double() : unbls::matrix() {}
unbls::matrix_double::matrix_double(size_t nrows, size_t ncols) :
  unbls::matrix()
{
  resize(nrows, ncols);
}
unbls::matrix_double::matrix_double(size_t nrows, size_t ncols, const std::vector<double>& i_vals) :
  unbls::matrix()
{
  resize(nrows, ncols);
  vals = i_vals;
}
void unbls::matrix_double::resize(size_t nrows, size_t ncols) 
{
  _nrows = nrows;
  _ncols = ncols;
  vals.resize(_nrows*_ncols);
  unbls::zero(vals);
}
double* unbls::matrix_double::address()
{
  return &vals[0];
}

// matrix_mask methods
unbls::matrix_mask::matrix_mask(size_t nrows, size_t ncols) :
  unbls::matrix()
{
  resize(nrows, ncols);
}
unbls::matrix_mask::matrix_mask(size_t nrows, size_t ncols, const std::vector<uint8_t>& i_vals) :
  unbls::matrix()
{
  resize(nrows, ncols);
  vals = i_vals;
}
void unbls::matrix_mask::resize(size_t nrows, size_t ncols) 
{
  _nrows = nrows;
  _ncols = ncols;
  vals.resize(_nrows*_ncols);
  unbls::zero(vals);
}
uint8_t* unbls::matrix_mask::address() 
{
  return &vals[0];
}


