#ifndef VUTILS__H
#define VUTILS__H

/*! \file 
   \brief utilities around vector and matrix algebra
   */

//! return the median of the array. The array is mixed up. 
/*!  Use next function if you want to keep your array untouched  */
double DArrayMedian(double *array, const int size);

//! same as above but does not mix up the array. (calls previous one on a copy of its input ... ). 
double DConstArrayMedian(const double *array, const int size);

#endif
