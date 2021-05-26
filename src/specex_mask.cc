
#include <specex_unbls.h>

#include <specex_mask.h>
#include <specex_psf.h>
#include <specex_image_data.h>


specex::Mask::Mask() {
}

void specex::Mask::ApplyMaskToImage(specex::image_data& img, const specex::PSF& psf, const double& value) const{

  int hsize = 5;
  for(size_t m=0;m<WaveIntervals.size();m++) {
    const Interval& inter=WaveIntervals[m];
    
    for(map<int,specex::Trace>::const_iterator it = psf.FiberTraces.begin(); it!=psf.FiberTraces.end(); it++) {
      
      int i1 = int(it->second.X_vs_W.Value(inter.min));
      int j1 = int(it->second.Y_vs_W.Value(inter.min));
      int i2 = int(it->second.X_vs_W.Value(inter.max));
      int j2 = int(it->second.Y_vs_W.Value(inter.max));
      
      int imin = max(0,min(i1,i2)-hsize);
      int imax = min(int(img.Nx())-1,max(i1,i2)+hsize);
      int jmin = max(0,min(j1,j2));
      int jmax = min(int(img.Ny())-1,max(j1,j2));
      
      
      //if(it==PSF.FiberTraces.begin()) {
      
      //cout << "DEBUG specex::PSF_Fitter::FitSeveralSpots mask #" << m << " fiber #" << it->first << " : region [" << imin << ":" << imax << "," << jmin << ":" << jmax << "] for first fiber" << endl;
      
      
	
      for(int j=jmin;j<=jmax;j++)
	for(int i=imin;i<=imax;i++)
	  img(i,j)=value; 
    }
  }
}
