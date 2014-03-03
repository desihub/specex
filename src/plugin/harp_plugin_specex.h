// @COPYRIGHT@

// To get the usual type definitions from HARP...
#include <harp.hpp>

// This header must be included by external plugins
#include <harp/plugin.hpp>


#include <iostream>


#include <specex_psf.h>


namespace harp {

  // example spec class which does nothing.  This is just to show how to compile
  // an external plugin.
  
  class specex_psf : public psf {

    friend class boost::serialization::access;

    
    public :

      specex_psf ( ) : psf () { }

      specex_psf ( boost::property_tree::ptree const & props );

      ~specex_psf ( ) { }

      // overloaded virtual methods from base class

      boost::property_tree::ptree metadata ( ) const { return boost::property_tree::ptree(); }

      size_t n_spec ( ) const { return nspec_; }

      size_t n_lambda ( ) const { return lambda_.size(); }
      
      size_t img_rows ( ) const { return rows_; }
      
      size_t img_cols ( ) const { return cols_; }
      
      vector_double lambda ( ) const { return lambda_; }

      void response ( size_t spec, size_t lambda, size_t & x_offset, size_t & y_offset, matrix_double & patch ) const ;

      size_t response_nnz_estimate ( ) const;

    private :

      template < class Archive >
      void serialize ( Archive & ar, const unsigned int version ) {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(psf);
        ar & BOOST_SERIALIZATION_NVP(nspec_);
        ar & BOOST_SERIALIZATION_NVP(lambda_);
        ar & BOOST_SERIALIZATION_NVP(rows_);
        ar & BOOST_SERIALIZATION_NVP(cols_);
        return;
      }

      vector_double lambda_;
      size_t nspec_;
      size_t rows_;
      size_t cols_;
      specex::PSF_p actual_specex_psf;
      std::map<int,int> bundle_; // bundles of fibers;
  };
  BOOST_SERIALIZATION_SHARED_PTR(specex_psf)

  psf * specex_psf_create ( boost::property_tree::ptree const & props );

}

// This global function in the shared object must be unique and outside the harp namespace.
// this is the function actually called by the plugin registry.

extern "C" {
  void initialize ( void * registry );
}

