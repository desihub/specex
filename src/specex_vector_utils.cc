#include <math.h>
#include <algorithm>
#include <string.h> // memcpy 

#include "specex_vector_utils.h"

using namespace std;

double DArrayMedian(double *array, const int size)
{
  sort(array, array+size);
  return size&1? array[size/2] : (array[size/2-1] + array[size/2])*0.5;
}

double DConstArrayMedian(const double *array, const int size)
{
  double *to_delete = new double[size*sizeof(double)];
  memcpy(to_delete, array, size*sizeof(double));
  double median = DArrayMedian(to_delete, size);
  delete [] to_delete;
  return median;
}
 

float FArrayMedian(float *array, const int size)
{
  sort(array, array+size);
  return size&1? array[size/2] : (array[size/2-1] + array[size/2])*0.5;
}


void Dmean_median_sigma(double *values, const int nval, double &mean,  double &median, double &sigma)
{
  mean =0;
  sigma =0;
  median = DArrayMedian(values, nval);
  for (int i=nval-1; i >= 0 ; --i)
    {
      mean += values[i];
      sigma += values[i]*values[i];
    }
  mean /= double(nval);
  sigma = sigma/double(nval) - mean*mean;
  if (sigma>0)  sigma = sqrt(sigma); else sigma = 0;
}

float Fmedian_sigma(float *values, const int nval, float &sigma)
{
  double dmean =0;
  double dsigma =0;
  double median = FArrayMedian(values, nval);
  for (int i=nval-1; i>=0 ; --i)
    {
      dmean += values[i];
      dsigma += values[i]*values[i];
    }
  dmean /= double(nval);
  dsigma = dsigma/double(nval) - dmean * dmean;
  if (dsigma>0)  sigma = sqrt(dsigma); else sigma = 0;
  return median;
}





