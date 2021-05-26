#ifndef SPECEX_IMAGE_DATA_BASE__H
#define SPECEX_IMAGE_DATA_BASE__H

// Base image class for specex adapted from the HARP package developed
// by Ted Kisner, see https://github.com/tskisner/HARP

#include <unbls.h>

class image_data_base {  

 public :
  
  image_data_base ( );
  
  virtual ~image_data_base ( );
  
  virtual size_t n_rows ( ) const;
  
  virtual size_t n_cols ( ) const;
  
  virtual void values ( unbls::vector_double & data ) const;
  
  virtual void inv_variance ( unbls::vector_double & invvar ) const;
  
  virtual void mask ( unbls::vector_mask & msk ) const;
  
  void values ( unbls::matrix_double & data ) const;
  
  void inv_variance ( unbls::matrix_double & invvar ) const;
  
  void mask ( unbls::matrix_mask & msk ) const;
  
 private :
  
  std::string type_;
    
};

#endif

