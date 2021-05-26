#ifndef SPECEX_HERMITE__H
#define SPECEX_HERMITE__H

#include <unbls.h>

namespace specex {
  double HermitePol(const int Degree, const double &x);
  double HermitePolDerivative(const int Degree, const double &x);
  void HermitePols(unbls::vector_double& H, const int Degree, const double &x);
  void HermitePolsAndDerivatives(unbls::vector_double& H, unbls::vector_double& dHdx, const int Degree, const double &x);

}

#endif
