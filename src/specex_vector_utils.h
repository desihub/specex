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

//!
float FArrayMedian(float *array, const int size);

/*! fit a gaussian on a region about mean of half width k-sigma, return mean
 * if first_evalutation, mean and sigma are guessed with  DConst_mean_median_sigma (using median)
 */
double gaussianfit(const double *values , int nval, double &mean, double &sigma, const double &k=3.5, bool first_evalutation=true);


//template<class T>  T ScalProd(const T A[], const T B[], int N);
#endif /* VUTILS__H */
