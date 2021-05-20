#ifndef SPECEX_IMAGE_DATA__H
#define SPECEX_IMAGE_DATA__H

#define CHECK_BOUNDS

#include <specex_message.h>
#include <specex_image_data_base.h>

namespace specex {
  class image_data : public image_data_base {

  protected :
    size_t rows_;
    size_t cols_;
    
  public :

    unbls::vector_double data;
    
    image_data ();
    image_data ( size_t ncols, size_t nrows);
    image_data ( size_t ncols, size_t nrows, const unbls::vector_double& i_data);
    void resize( size_t ncols, size_t nrows); 
    size_t n_rows ( ) const { return rows_; }
    size_t n_cols ( ) const { return cols_; }
    void values ( unbls::vector_double & i_data ) const {i_data=data;}
    size_t Ny ( ) const { return rows_; }
    size_t Nx ( ) const { return cols_; }
    
    
    inline double& operator()(const int i, const int j) {
#ifdef CHECK_BOUNDS
      if (i<0 || i>=cols_ || j<0 || j>=rows_)
	SPECEX_ERROR("Out of range");
#endif
      return data[i+j*cols_]; // "STANDARD" PACKING (FITSIO)      
      
    }
    
    inline const double& operator()(const int i, const int j) const {
#ifdef CHECK_BOUNDS
      if (i<0 || i>=cols_ || j<0 || j>=rows_)
	SPECEX_ERROR("Out of range");
#endif
      
      return data[i+j*cols_];// "STANDARD" PACKING (FITSIO)
      
    }

  };
};

#endif
