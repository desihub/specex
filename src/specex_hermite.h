#ifndef SPECEX_HERMITE__H
#define SPECEX_HERMITE__H

#include <specex_unbls.h>

namespace specex {
  double HermitePol(const int Degree, const double &x);
  double HermitePolDerivative(const int Degree, const double &x);
}

#endif
