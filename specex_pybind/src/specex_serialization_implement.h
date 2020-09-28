#ifndef SPECEX_SERIALIZATION_IMPLEMENT__H
#define SPECEX_SERIALIZATION_IMPLEMENT__H

// This file registers the derived classes that need to be
// serialized.  This should be put at the end of the source
// file containing main(), or in the main entry point of
// the loadable module.

BOOST_CLASS_EXPORT(specex::PSF_Params)
BOOST_CLASS_EXPORT(specex::PSF)
BOOST_CLASS_EXPORT(specex::GaussHermitePSF)
BOOST_CLASS_EXPORT(specex::GaussHermite2PSF)
BOOST_CLASS_EXPORT(specex::HatHermitePSF)
BOOST_CLASS_EXPORT(specex::HatMoffatPSF)

#endif
