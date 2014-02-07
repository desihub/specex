#ifndef SPECEX_SERIALIZATION__H
#define SPECEX_SERIALIZATION__H

// this file needs to be included only before main() 

#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>

#include <specex_psf.h>
#include <specex_gauss_hermite_psf.h>
#include <specex_hat_hermite_psf.h>



BOOST_CLASS_EXPORT(specex::PSF_Params)
BOOST_CLASS_EXPORT(specex::PSF)
BOOST_CLASS_EXPORT(specex::GaussHermitePSF)
BOOST_CLASS_EXPORT(specex::HatHermitePSF)

#endif
