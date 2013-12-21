#ifndef SPECEX_SERIALIZATION__H
#define SPECEX_SERIALIZATION__H

// this file needs to be included only before main() 

#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>

BOOST_CLASS_EXPORT(specex::PSF_Params)
BOOST_CLASS_EXPORT(specex::PSF)
BOOST_CLASS_EXPORT(specex::GaussHermitePSF)

#endif
