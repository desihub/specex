#ifndef SPECEX_PYOPTIONS__H
#define SPECEX_PYOPTIONS__H

#define GETOPT

#include <vector>
#include <string>
#include <specex_unbls.h>
#include <specex_psf_fitter.h>


#ifndef GETOPT
#include <boost/program_options.hpp>
namespace popts = boost::program_options;
#endif

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
        
#ifndef GETOPT
    popts::variables_map vm;    
    popts::options_description desc;
#endif
    
    bool fit_traces;
    bool fit_sigmas;
    bool fit_thepsf;
    bool single_bundle;
    bool write_tmp_results;
    bool fit_psf_tails;
    bool fit_continuum;
    bool use_variance_model;
    bool fit_individual_spots_position;
    
    int parse(int argc, char *argv[] ); 

    PyOptions()
#ifndef GETOPT
      : desc(popts::options_description("Options"))
#endif
      {	
      arc_image_filename="";
      input_psf_filename="";
      output_xml_filename="";
      output_fits_filename="";
      output_spots_filename="";
      flux_hdu=1;
      ivar_hdu=2;
      mask_hdu=3;
      header_hdu=1;
      xtrace_hdu=-1;
      ytrace_hdu=-1;
    
      psf_model="GAUSSHERMITE";
      first_fiber_bundle=0;
      last_fiber_bundle=0;
      first_fiber=0;
      last_fiber=100000;
      half_size_x=8;
      half_size_y=5;  
      gauss_hermite_deg=6;
      gauss_hermite_deg2=2; 
      legendre_deg_wave=3;
      legendre_deg_x=1;
      trace_deg_wave=6;
      trace_deg_x=6;  
      trace_prior_deg=0;
      psf_error=0;
      psf_core_wscale=0;
      max_number_of_lines=200; 
    
      broken_fibers_string="";

      lamp_lines_filename="";

      fit_traces = true;
      fit_sigmas = true;
      fit_thepsf = true;
      single_bundle = false;
      write_tmp_results = false;
      fit_psf_tails = false;
      fit_continuum = false;
      use_variance_model = false;
      fit_individual_spots_position = false;
      
      return;
      
    }
    
  };
  
}

#endif
