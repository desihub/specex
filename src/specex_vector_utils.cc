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





