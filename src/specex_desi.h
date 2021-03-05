
#ifndef SPECEX_DESI_H
#define SPECEX_DESI_H

#include <specex_pyoptions.h>
#include <specex_pyio.h>
#include <specex_pyprior.h>

int fit_psf(
			     specex::PyOptions,
			     specex::PyIO,
			     specex::PyPrior,
			     specex::PyImage,
			     specex::PyPSF
			      );

#endif

