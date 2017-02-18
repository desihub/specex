#ifndef SPECEX_SERIALIZATION__H
#define SPECEX_SERIALIZATION__H

// This file includes all archives used for serialization
// before the includes for the classes that need to be serialized.

#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <specex_trace.h>
#include <specex_psf.h>
#include <specex_gauss_hermite_psf.h>
#include <specex_gauss_hermite_two_psf.h>
#include <specex_hat_hermite_psf.h>
#include <specex_hat_moffat_psf.h>
#include <specex_disk_moffat_psf.h>


#endif
