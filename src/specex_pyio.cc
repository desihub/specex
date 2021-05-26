#include <string>
#include <iostream>
#include <specex_message.h>
#include <specex_spot.h>

#include <specex_pyio.h>
#include <specex_pyimage.h>
#include <specex_psf_io.h>
#include <specex_psf_proc.h>

using namespace std;

int specex::PyIO::read_preproc(specex::PyOptions opts, specex::PyImage& pyimg){

  read_DESI_preprocessed_image(opts.arc_image_filename,pyimg.image,pyimg.weight,
			       pyimg.mask,pyimg.rdnoise,pyimg.header);

  return EXIT_SUCCESS;
}

int specex::PyIO::read_img_datam(specex::PyOptions opts,
				 specex::image_data& image,
				 specex::image_data& weight,
				 specex::image_data& mask,
				 specex::image_data& rdnoise,
				 std::map<std::string,std::string>& header
				 )
				 {

  read_DESI_preprocessed_image(opts.arc_image_filename,image,weight,
			       mask,rdnoise,header);

  return EXIT_SUCCESS;
}

int specex::PyIO::load_psf(specex::PyOptions opts, specex::PyPSF& pypsf){

  load_psf_work(pypsf.psf);  
  
  return EXIT_SUCCESS;

}

int specex::PyIO::write_spots(specex::PyOptions opts, specex::PyPSF& pypsf){

  vector <Spot_p> fitted_spots = pypsf.fitted_spots;

  // xml serialization is removed nothing done here
  //if(opts.output_spots_filename != "")
  //  write_spots_xml(fitted_spots,opts.output_spots_filename);    
  
  return EXIT_SUCCESS;

}

int specex::PyIO::read_psf(specex::PyOptions opts, specex::PyPSF& pypsf){

  if( ! use_input_specex_psf ) {
    std::cout << "read_psf 1" << std::endl;
    SPECEX_INFO("Initializing a " << opts.psf_model << " PSF");
    pypsf.psf = PSF_p(new specex::GaussHermitePSF(opts.gauss_hermite_deg));
    
    pypsf.psf->hSizeX = opts.half_size_x;
    pypsf.psf->hSizeY = opts.half_size_y;
    SPECEX_INFO("trace_deg_x=" << opts.trace_deg_x << " trace_deg_wave=" <<
		opts.trace_deg_wave);
    read_traceset_fits(pypsf.psf,opts.input_psf_filename,opts.trace_deg_x,
		       opts.trace_deg_wave);     
    
  }else{ // use_input_specex_psf
    std::cout << "read_psf 2" << std::endl;
    read_psf_fits(pypsf.psf,opts.input_psf_filename);
  }

  return EXIT_SUCCESS;
  
}

int specex::PyIO::set_inputpsf(specex::PyOptions opts){
  
  use_input_specex_psf = true;

  psf_change_req |= (! opts.vm["half-size-x"].defaulted());
  psf_change_req |= (! opts.vm["half-size-y"].defaulted());
  psf_change_req |= (! opts.vm["gauss-hermite-deg"].defaulted());
  psf_change_req |= (! opts.vm["gauss-hermite-deg2"].defaulted());
  psf_change_req |= (! opts.vm["legendre-deg-wave"].defaulted());
  psf_change_req |= (! opts.vm["legendre-deg-x"].defaulted());
  psf_change_req |= (! opts.vm["trace-deg-wave"].defaulted());
  psf_change_req |= (! opts.vm["trace-deg-x"].defaulted());
  
  if(psf_change_req && use_input_specex_psf) {
    SPECEX_WARNING("option(s) were given to specify the psf properties, so we cannot use the input PSF parameters as a starting point (except for the trace coordinates)");
    use_input_specex_psf = false;
  }
  if(use_input_specex_psf) {
    SPECEX_INFO("will start from input psf parameters")
  }

  return EXIT_SUCCESS;
  
}

