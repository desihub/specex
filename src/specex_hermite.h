#ifndef SPECEX_HERMITE__H
#define SPECEX_HERMITE__H

#include <unhrp.h>

namespace specex {
  double HermitePol(const int Degree, const double &x);
  double HermitePolDerivative(const int Degree, const double &x);
  void HermitePols(unhrp::vector_double& H, const int Degree, const double &x);
  void HermitePolsAndDerivatives(unhrp::vector_double& H, unhrp::vector_double& dHdx, const int Degree, const double &x);

}

#endif
