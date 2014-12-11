#ifndef SPECEX_IMAGE_DATA__H
#define SPECEX_IMAGE_DATA__H

#define CHECK_BOUNDS

#include <harp.hpp>
#include <specex_message.h>

namespace specex {
  class image_data : public harp::image {
    friend class boost::serialization::access;

  protected :
    size_t rows_;
    size_t cols_;
    

  public :

    harp::vector_double data;
    
    image_data ();
    image_data ( size_t ncols, size_t nrows);
    image_data ( size_t ncols, size_t nrows, const harp::vector_double& i_data);
    void resize( size_t ncols, size_t nrows); 
    size_t n_rows ( ) const { return rows_; }
    size_t n_cols ( ) const { return cols_; }
    void values ( harp::vector_double & i_data ) const {i_data=data;}
    size_t Ny ( ) const { return rows_; }
    size_t Nx ( ) const { return cols_; }
    
    
    inline double& operator()(const int i, const int j) {
#ifdef CHECK_BOUNDS
      if (i<0 || i>=cols_ || j<0 || j>=rows_)
	SPECEX_ERROR("Out of range");
#endif
      return data[i+j*cols_]; // "STANDARD" PACKING (FITSIO)
      // return data[j+i*rows_]; // "REVERSE" PACKING (HARP)
      
      
    }
    
    inline const double& operator()(const int i, const int j) const {
#ifdef CHECK_BOUNDS
      if (i<0 || i>=cols_ || j<0 || j>=rows_)
	SPECEX_ERROR("Out of range");
#endif
      
      return data[i+j*cols_];// "STANDARD" PACKING (FITSIO)
      // return data[j+i*rows_]; // "REVERSE" PACKING (HARP)
      
    }

  };
};
BOOST_SERIALIZATION_SHARED_PTR(image_data)

#endif
