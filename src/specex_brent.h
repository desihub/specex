#ifndef BRENT__H
#define BRENT__H

#include <math.h>

typedef double (AnalyticFunction)(const double &, const void*);

double brent(AnalyticFunction * f, double ax, double bx, double cx, double tol, const void *context, double& fval, int& status, int maxiter=100);

#endif
