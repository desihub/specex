// Base image class for specex adapted from the HARP package developed
// by Ted Kisner, see https://github.com/tskisner/HARP

#include <unbls.h>
#include <specex_image_data_base.h>

using namespace std;

image_data_base::image_data_base ( ) {
  type_ = "";
}

image_data_base::~image_data_base ( ) {
}

size_t image_data_base::n_rows ( ) const {
  return 0;
}


size_t image_data_base::n_cols ( ) const {
  return 0;
}


void image_data_base::values ( unbls::vector_double & data ) const {
  return;
}


void image_data_base::inv_variance ( unbls::vector_double & invvar ) const {
  return;
}

void image_data_base::mask ( unbls::vector_mask & msk ) const {
  return;
}

void image_data_base::values ( unbls::matrix_double & data ) const {

  size_t imgrows = n_rows();
  size_t imgcols = n_cols();

  size_t nelem = imgrows * imgcols;

  data.resize ( imgrows, imgcols );

  unbls::vector_double tempdata ( nelem );

  values ( tempdata );

  for ( size_t i = 0; i < imgcols; ++i ) {
    for ( size_t j = 0; j < imgrows; ++j ) {
      data( j, i ) = tempdata[ i * imgrows + j ];
    }
  }

  return;
}


void image_data_base::inv_variance ( unbls::matrix_double & invvar ) const {

  size_t imgrows = n_rows();
  size_t imgcols = n_cols();

  size_t nelem = imgrows * imgcols;

  invvar.resize ( imgrows, imgcols );

  unbls::vector_double tempvar ( nelem );

  inv_variance ( tempvar );

  for ( size_t i = 0; i < imgcols; ++i ) {
    for ( size_t j = 0; j < imgrows; ++j ) {
      invvar( j, i ) = tempvar[ i * imgrows + j ];
    }
  }

  return;
}

void image_data_base::mask ( unbls::matrix_mask & msk ) const {

  size_t imgrows = n_rows();
  size_t imgcols = n_cols();

  size_t nelem = imgrows * imgcols;

  msk.resize ( imgrows, imgcols );

  unbls::vector_mask tempmask ( nelem );

  mask ( tempmask );

  for ( size_t i = 0; i < imgcols; ++i ) {
    for ( size_t j = 0; j < imgrows; ++j ) {
      msk( j, i ) = tempmask[ i * imgrows + j ];
    }
  }

  return;
}

