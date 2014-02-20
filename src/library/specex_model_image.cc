#include <specex_model_image.h>

using namespace std;


specex::Stamp specex::compute_stamp(const specex::image_data& model_image, const specex::PSF_p psf, std::vector<specex::Spot_p>& spots) {
  
  Stamp global_stamp(model_image);
  global_stamp.begin_i = 10000;
  global_stamp.end_i   = 0;
  global_stamp.begin_j = 10000;
  global_stamp.end_j   = 0;
  
  for(size_t s=0;s<spots.size();s++) {
    
    specex::Spot_p spot = spots[s];
      
    Stamp stamp(model_image);
    psf->StampLimits(spot->xc,spot->yc,stamp.begin_i,stamp.end_i,stamp.begin_j,stamp.end_j);
    stamp.begin_i = max(0,stamp.begin_i);
    stamp.end_i   = min(stamp.Parent_n_cols(),stamp.end_i);
    stamp.begin_j = max(0,stamp.begin_j);
    stamp.end_j   = min(stamp.Parent_n_rows(),stamp.end_j);
    
    global_stamp.begin_i = min(global_stamp.begin_i,stamp.begin_i);
    global_stamp.end_i   = max(global_stamp.end_i,stamp.end_i);
    global_stamp.begin_j = min(global_stamp.begin_j,stamp.begin_j);
    global_stamp.end_j   = max(global_stamp.end_j,stamp.end_j);
  }
  return global_stamp;
}

specex::Stamp specex::compute_model_image(specex::image_data& model_image, const specex::PSF_p psf, std::vector<specex::Spot_p>& spots, const int first_fiber, const int last_fiber, bool only_on_spots, bool only_psf_core, bool only_positive) {
  
  Stamp global_stamp = compute_stamp(model_image,psf,spots);
  
  model_image.clear();
  
  vector<specex::Stamp> spot_stamps;
  image_data spot_stamp_footprint; // because can overlap
  if(only_on_spots) {
    spot_stamp_footprint = image_data(model_image.n_cols(),model_image.n_rows());
  }

  for(size_t s=0;s<spots.size();s++) {
    
    specex::Spot_p spot = spots[s];
      
    Stamp stamp(model_image);
    psf->StampLimits(spot->xc,spot->yc,stamp.begin_i,stamp.end_i,stamp.begin_j,stamp.end_j);
    stamp.begin_i = max(0,stamp.begin_i);
    stamp.end_i   = min(stamp.Parent_n_cols(),stamp.end_i);
    stamp.begin_j = max(0,stamp.begin_j);
    stamp.end_j   = min(stamp.Parent_n_rows(),stamp.end_j);
    spot_stamps.push_back(stamp);
   
    if(only_on_spots) {
      for(int j=stamp.begin_j;j<stamp.end_j;j++)
	for(int i=stamp.begin_i;i<stamp.end_i;i++)
	  spot_stamp_footprint(i,j)=1;
    }
  }
  
  
#ifdef CONTINUUM
  for (int j=global_stamp.begin_j; j <global_stamp.end_j; ++j) {  
    
    
    for(std::map<int,PSF_Params>::iterator bundle_it = psf->ParamsOfBundles.begin(); bundle_it != psf->ParamsOfBundles.end(); bundle_it++) {
      
      const PSF_Params * psf_params = &(bundle_it->second);
      
      int begin_i = max(global_stamp.begin_i, int(floor(psf->GetTrace(first_fiber).X_vs_Y.Value(double(j))+0.5))-psf->hSizeX-1);
      int end_i   = min(global_stamp.end_i  , int(floor(psf->GetTrace(last_fiber).X_vs_Y.Value(double(j))+0.5))+psf->hSizeX+2);
      
      for(int fiber=first_fiber; fiber<=last_fiber; fiber++) {
	double x = psf->GetTrace(fiber).X_vs_Y.Value(double(j));
	double w = psf->GetTrace(fiber).W_vs_Y.Value(double(j));
	double continuum_flux = psf->ContinuumPol.Value(w);
	double expfact_for_continuum=continuum_flux/(2*M_PI*square(psf->continuum_sigma_x));
	if(expfact_for_continuum!=0) {
	  for (int i=begin_i ; i <end_i; ++i) {    

	    if(weight(i,j)<=0) continue;
	    if(only_on_spots && spot_stamp_footprint(i,j)==0) continue;
	    
	    double val = expfact_for_continuum*exp(-0.5*square((i-x)/psf->continuum_sigma_x));
	    model(i,j) += val;
	    if(val>0) model_for_var(i,j) += val;
	  } // end of loop on i
	} // expfact
      } // end of loop on fiber
    } // end of loop on bundles
  } // end of loop on j   
#endif  
  
  
  for(size_t s=0;s<spots.size();s++) {
    specex::Spot_p spot = spots[s];
    harp::vector_double spot_params = psf->AllLocalParamsXW(spot->xc,spot->wavelength,spot->fiber_bundle);
    
    for (int j=global_stamp.begin_j; j <global_stamp.end_j; ++j) {  
      
      int begin_i = max(global_stamp.begin_i, int(floor(psf->GetTrace(first_fiber).X_vs_Y.Value(double(j))+0.5))-psf->hSizeX-1);
      int end_i   = min(global_stamp.end_i  , int(floor(psf->GetTrace(last_fiber).X_vs_Y.Value(double(j))+0.5))+psf->hSizeX+2);
      
      for(int i=begin_i; i<end_i;i++) {
	
	if(weight(i,j)<=0) continue;
		
	bool in_core =  spot_stamps[s].Contains(i,j);
	if(!in_core && only_psf_core) continue;
	
	if(only_on_spots && spot_stamp_footprint(i,j)==0) continue;
	
	double val = spot->flux*psf->PSFValueWithParamsXY(spot->xc,spot->yc, i, j, spot_params, 0, 0, in_core, true); // compute CPU expensive PSF core only if needed

	if(val>0 || (!only_positive))
	  model_image(i,j) += val; // this now includes core and tails
	
	
      }
    }
  }
}
