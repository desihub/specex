#ifndef SPECEX_PYFITTING_H
#define SPECEX_PYFITTING_H

#include <specex_pyoptions.h>
#include <specex_pyio.h>
#include <specex_pyprior.h>

namespace specex {

  class PyFitting : public std::enable_shared_from_this <PyFitting> {

  public :

    typedef std::shared_ptr <PyFitting> pshr;
    
    int fit_psf(
		specex::PyOptions,
		specex::PyIO,
		specex::PyPrior,
		specex::PyImage,
		specex::PyPSF&
		);
    
    // int specex_desi_psf_fit_main(int argc, char *argv[]);
    // int specex_desi_psf_fit_main();

    PyFitting(){}

  };

}

#endif

