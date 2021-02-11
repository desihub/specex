#ifndef SPECEX_PYOPTIONS__H
#define SPECEX_PYOPTIONS__H

#include <boost/program_options.hpp>
#include <vector>
#include <string>

#include <harp.hpp>

#include <specex_psf_fitter.h>

namespace popts = boost::program_options;

namespace specex {
  
  class PyOptions : public std::enable_shared_from_this <PyOptions> {

  public :

    typedef std::shared_ptr <PyOptions> pshr;

    std::string arc_image_filename;
    std::string input_psf_filename; 
    std::string output_xml_filename; 
    std::string output_fits_filename; 
    std::string output_spots_filename; 
    int flux_hdu; 
    int ivar_hdu; 
    int mask_hdu; 
    int header_hdu; 
    int xtrace_hdu; 
    int ytrace_hdu; 
    
    std::string psf_model; 
    int    first_fiber_bundle; 
    int    last_fiber_bundle; 
    int    first_fiber; 
    int    last_fiber; 
    int    half_size_x; 
    int    half_size_y; 
    int gauss_hermite_deg; 
    int gauss_hermite_deg2; 
    int legendre_deg_wave; 
    int legendre_deg_x; 
    int trace_deg_wave; 
    int trace_deg_x; 
    int trace_prior_deg; 
    double psf_error; 
    double psf_core_wscale; 
    int max_number_of_lines; 
    
    std::string broken_fibers_string;

    std::string lamp_lines_filename;

    std::vector<std::string> argurment_priors;

    popts::variables_map vm;

    popts::options_description desc;

    int parse(int argc, char *argv[] ); 

    PyOptions()
      : desc(popts::options_description("Options")) 
      {
        
      arc_image_filename="";
      input_psf_filename="";
      output_xml_filename="";
      output_fits_filename="";
      output_spots_filename="";
      flux_hdu=0;
      ivar_hdu=0;
      mask_hdu=0;
      header_hdu=0;
      xtrace_hdu=0;
      ytrace_hdu=0;
    
      psf_model="";
      first_fiber_bundle=0;
      last_fiber_bundle=0;
      first_fiber=0;
      last_fiber=0;
      half_size_x=0;
      half_size_y=0;  
      gauss_hermite_deg=0;
      gauss_hermite_deg2=0; 
      legendre_deg_wave=0;
      legendre_deg_x=0;
      trace_deg_wave=0;
      trace_deg_x=0;  
      trace_prior_deg=0;
      psf_error=0;
      psf_core_wscale=0;
      max_number_of_lines=0; 
    
      broken_fibers_string="";

      lamp_lines_filename="";
      
      return;
      
    }
    
  };
  
}

#endif
