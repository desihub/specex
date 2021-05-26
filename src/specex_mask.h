#ifndef SPECEX_MASK__H
#define SPECEX_MASK__H

#include <vector>

namespace specex {

class PSF;
class image_data;

class Interval {
 public :
  double min;
  double max;
  Interval(const double& imin, const double& imax) :
    min(imin),
    max(imax)
    {
    }
};

class Mask {

 public :
  
  std::vector<Interval> WaveIntervals; // do not fit psf in those wavelength intervals, primarily because of missing lines
  
  Mask();
  void ApplyMaskToImage(image_data& img, const PSF& psf, const double& value=0) const;
  void Clear() { WaveIntervals.clear();}
};

}


#endif
