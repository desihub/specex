
#ifndef SPECEX_DESI_H
#define SPECEX_DESI_H

#include <specex_pyoptions.h>
#include <specex_pyio.h>
#include <specex_pyprior.h>

int specex_desi_psf_fit_main(
			     specex::PyOptions,
			     specex::PyIO,
			     specex::PyPrior,
			     specex::PyImage,
			     specex::PyPSF
			      );
// int specex_desi_psf_fit_main(int argc, char *argv[]);
//int specex_desi_psf_fit_main();

#endif

