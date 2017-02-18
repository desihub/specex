#include "specex_trace.h"
#include "specex_spectrograph.h"
#include "specex_vector_utils.h"
#include "specex_message.h"


using namespace std;

void specex::Spectrograph::AutoConfigure(const specex::TraceSet& traceset) {
  SPECEX_INFO("Guess number of bundles from traces");
  int nfibers = traceset.size();
  
  double* central_waves = new double[nfibers];
  for(int f=0;f<nfibers;f++)
    central_waves[f]=(traceset.find(f)->second.X_vs_W.xmin+traceset.find(f)->second.X_vs_W.xmax)/2.;
  double central_wave=DConstArrayMedian(central_waves,nfibers);
  SPECEX_INFO("Central wavelength             = " << central_wave);
  delete [] central_waves;
  
  double* spacing = new double[nfibers-1];
  for(int f=0;f<nfibers-1;f++) {
    spacing[f] = traceset.find(f+1)->second.X_vs_W.Value(central_wave)-traceset.find(f)->second.X_vs_W.Value(central_wave);
    //SPECEX_INFO("Spacing=" << spacing[f]);
  }
  double median_spacing = DConstArrayMedian(spacing,nfibers-1);
  SPECEX_INFO("Median distance between fibers = " << median_spacing);
  
  int number_of_bundles=0;
  int bundle_size=0;
  int first_fiber=0;
  for(int f=0;f<nfibers-1;f++) {
    if(spacing[f]>median_spacing*1.5) {
      // we have a bundle
      int current_bundle_size = f-first_fiber+1;
      number_of_bundles += 1;
      SPECEX_INFO("Bundle of size " << current_bundle_size);
      
      if(bundle_size==0) {
	bundle_size = current_bundle_size;
      }else{
	if(current_bundle_size != bundle_size) {
	  SPECEX_ERROR("cannot deal with varying bundle size");
	}
      }
      first_fiber=f+1;
    }
  }
  number_of_bundles += 1;
  if(number_of_bundles==1) { // there is only one 
    bundle_size = nfibers;
  }
  // save result
  number_of_fiber_bundles_per_ccd = number_of_bundles;
  number_of_fibers_per_bundle = bundle_size;

  SPECEX_INFO("number of fibers in CCD        = " << nfibers);
  SPECEX_INFO("number of fiber bundles in CCD = " << number_of_bundles);
  SPECEX_INFO("number of fibers per bundle    = " << number_of_fibers_per_bundle);
  

  delete[] spacing;
  
}
