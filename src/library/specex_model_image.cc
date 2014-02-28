#include <specex_model_image.h>
#include <specex_message.h>
#include <specex_psf.h>

using namespace std;


specex::Stamp specex::compute_stamp(const specex::image_data& model_image, const specex::PSF_p psf, const std::vector<specex::Spot_p>& spots, int x_margin, int y_margin, int only_this_bundle) {
  
  Stamp global_stamp(model_image);
  global_stamp.begin_i = 10000;
  global_stamp.end_i   = 0;
  global_stamp.begin_j = 10000;
  global_stamp.end_j   = 0;
  
  for(size_t s=0;s<spots.size();s++) {
    
    specex::Spot_p spot = spots[s];

    if((only_this_bundle>=0) && (spot->fiber_bundle != only_this_bundle)) continue;

    Stamp stamp(model_image);
    psf->StampLimits(spot->xc,spot->yc,stamp.begin_i,stamp.end_i,stamp.begin_j,stamp.end_j);
    stamp.begin_i = max(0,stamp.begin_i);
    stamp.end_i   = min(stamp.Parent_n_cols(),stamp.end_i);
    stamp.begin_j = max(0,stamp.begin_j);
    stamp.end_j   = min(stamp.Parent_n_rows(),stamp.end_j);
    
    global_stamp.begin_i = min(global_stamp.begin_i,stamp.begin_i-x_margin);
    global_stamp.end_i   = max(global_stamp.end_i,stamp.end_i+x_margin);
    global_stamp.begin_j = min(global_stamp.begin_j,stamp.begin_j-y_margin);
    global_stamp.end_j   = max(global_stamp.end_j,stamp.end_j+y_margin);
  }
  global_stamp.begin_i = max(0,global_stamp.begin_i);
  global_stamp.end_i = min(global_stamp.end_i,global_stamp.Parent_n_cols());
  global_stamp.begin_j = max(0,global_stamp.begin_j);
  global_stamp.end_j = min(global_stamp.end_j,global_stamp.Parent_n_rows());
  

  //SPECEX_INFO("stamp [" << global_stamp.begin_i << ":" << global_stamp.end_i << "," << global_stamp.begin_j << ":" << global_stamp.end_j << "]");

  return global_stamp;
}

void specex::parallelized_compute_model_image(specex::image_data& model_image, const specex::image_data& weight, const specex::PSF_p psf, const std::vector<specex::Spot_p>& spots, bool only_on_spots, bool only_psf_core, bool only_positive, int x_margin, int y_margin, int only_this_bundle) {

  SPECEX_INFO("parallelized_compute_model_image");
  
  int number_of_image_chuncks = 1;
  char* OMP_NUM_THREADS = getenv("OMP_NUM_THREADS");
  if(OMP_NUM_THREADS) {
    number_of_image_chuncks = atoi(OMP_NUM_THREADS);
    SPECEX_INFO("Using " << number_of_image_chuncks << " image chunks equal to value of OMP_NUM_THREADS");
  }

#ifdef EXTERNAL_TAIL  
  // precompute tail profile
  psf->TailProfile(0,0,psf->AllLocalParamsFW(spots[0]->fiber,spots[0]->wavelength,spots[0]->fiber_bundle));
#endif

  //SPECEX_INFO("specex::parallelized_compute_model_image");
  Stamp stamp = compute_stamp(model_image,psf,spots,x_margin,y_margin);
  int step_j  = (stamp.end_j-stamp.begin_j)/number_of_image_chuncks;
  
  model_image.data.clear();
  
  int chunk=0;
  
#pragma omp parallel for 
  for(chunk=0; chunk<number_of_image_chuncks; chunk++) {

    int begin_j = stamp.begin_j + chunk*step_j;
    int end_j   = stamp.begin_j + (chunk+1)*step_j;
    if(chunk==number_of_image_chuncks-1) end_j = stamp.end_j;
    
    if(end_j>begin_j)
      compute_model_image(model_image,weight,psf,spots,only_on_spots,only_psf_core,only_positive,begin_j,end_j,x_margin,y_margin,only_this_bundle);
  } 
}



void specex::compute_model_image(specex::image_data& model_image, const specex::image_data& weight, const specex::PSF_p psf, const std::vector<specex::Spot_p>& spots, bool only_on_spots, bool only_psf_core, bool only_positive, int predefined_begin_j, int predefined_end_j, int x_margin, int y_margin, int only_this_bundle) {
  
  
  Stamp global_stamp = compute_stamp(model_image,psf,spots,x_margin,y_margin,only_this_bundle);
  if(global_stamp.end_i==0) {
    SPECEX_WARNING("empty global stamp (can occur in parallel processing)");
    return;
  }
  //SPECEX_INFO("specex::compute_model_image [" << global_stamp.begin_i << ":" << global_stamp.end_i << "," << predefined_begin_j << ":" << predefined_end_j << "]");
  
  
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

    if(only_this_bundle>=0 && (spot->fiber_bundle != only_this_bundle)) continue;

    if(only_on_spots) {
      for(int j=stamp.begin_j;j<stamp.end_j;j++)
	for(int i=stamp.begin_i;i<stamp.end_i;i++)
	  spot_stamp_footprint(i,j)=1;
    }
  }
  
  int begin_j = global_stamp.begin_j;
  int end_j   = global_stamp.end_j;
  
  if(predefined_begin_j>=0)
    begin_j = max(begin_j,predefined_begin_j);
  if(predefined_end_j>=0)
    end_j = min(end_j,predefined_end_j);
  			     
  
  
#ifdef EXTERNAL_TAIL  
  int psf_tail_index = psf->ParamIndex("TAILAMP");
#endif
  double eps=1.e-20;

  // process per bundle
  int count=0;
  for(std::map<int,PSF_Params>::iterator bundle_it = psf->ParamsOfBundles.begin(); bundle_it != psf->ParamsOfBundles.end(); bundle_it++ , count++) {
    
    const PSF_Params& params_of_bundle = bundle_it->second;
    if(only_this_bundle>=0 && (params_of_bundle.bundle_id != only_this_bundle)) continue;
    
    int margin = min(MAX_X_MARGIN,psf->hSizeX); // 7 is half distance between center of ext. fibers of adjacent bundles
    if(x_margin>0) margin=x_margin;
    
    for (int j=begin_j; j <end_j; ++j) { 
      
      int begin_i = max(global_stamp.begin_i, int(floor(psf->GetTrace(params_of_bundle.fiber_min).X_vs_Y.Value(double(j))+0.5))-margin);
      int end_i   = min(global_stamp.end_i  , int(floor(psf->GetTrace(params_of_bundle.fiber_max).X_vs_Y.Value(double(j))+0.5))+margin+1);

      // zero the frame
      for (int i=begin_i ; i <end_i; ++i) {
	model_image(i,j) = 0;
      }
      
      
#ifdef CONTINUUM
      
      if(params_of_bundle.ContinuumPol.coeff(0)!=0) {
	for(int fiber=params_of_bundle.fiber_min; fiber<=params_of_bundle.fiber_max; fiber++) {
	  double x = psf->GetTrace(fiber).X_vs_Y.Value(double(j));
	  double w = psf->GetTrace(fiber).W_vs_Y.Value(double(j));
	  double continuum_flux = params_of_bundle.ContinuumPol.Value(w);
	  double expfact_for_continuum=continuum_flux/(sqrt(2*M_PI)*params_of_bundle.continuum_sigma_x);
	  if(expfact_for_continuum!=0) {
	    for (int i=begin_i ; i <end_i; ++i) {    
	      if(weight(i,j)<=0) continue;
	      if(only_on_spots && spot_stamp_footprint(i,j)==0) continue;
	      
	      double val = expfact_for_continuum*exp(-0.5*square((i-x)/params_of_bundle.continuum_sigma_x));
	      if(val>0 || (!only_positive))
		model_image(i,j) += val;
	      else if(model_image(i,j)==0) 
		model_image(i,j) = eps;
	      
	    } // end of loop on i
	  } // expfact
	} // end of loop on fiber
      } // end of test on continuum
#endif  
    
    } // end of loop on j
      
    int nspots_in_stamp = 0;
    for(size_t s=0;s<spots.size();s++) {
      if(spots[s]->fiber_bundle == bundle_it->first) nspots_in_stamp ++;
    }

    // loop on spots
    int sok=0;
    for(size_t s=0;s<spots.size();s++) {

      

      specex::Spot_p spot = spots[s];
      if(spot->fiber_bundle != bundle_it->first) continue; // keep only spots of this bundle
      
      //if(sok%100==0) SPECEX_INFO("spot " << sok << "/" << nspots_in_stamp);
      sok++;
      
      const Stamp& spot_stamp=spot_stamps[s];

      harp::vector_double spot_params = psf->AllLocalParamsXW(spot->xc,spot->wavelength,spot->fiber_bundle);
      bool has_tail  = spot_params(psf_tail_index)!=0;
      bool only_core = ( only_psf_core || (!has_tail) );
      
      for (int j=begin_j; j <end_j; ++j) { 
	
	int margin = min(MAX_X_MARGIN,psf->hSizeX); 
 	if(x_margin>0) margin=x_margin;
	
	int begin_i = max(global_stamp.begin_i, int(floor(psf->GetTrace(params_of_bundle.fiber_min).X_vs_Y.Value(double(j))+0.5))-margin);
	int end_i   = min(global_stamp.end_i  , int(floor(psf->GetTrace(params_of_bundle.fiber_max).X_vs_Y.Value(double(j))+0.5))+margin+1);
	
	if(only_core) {
	  if(j<spot_stamp.begin_j || j>=spot_stamp.end_j ) continue;
	  begin_i = max(begin_i,spot_stamp.begin_i);
	  end_i   = min(end_i,spot_stamp.end_i);	  
	}

	for(int i=begin_i; i<end_i;i++) {
	  
	  if(weight(i,j)<=0) continue;
	  
	  bool in_core =  spot_stamp.Contains(i,j);
	  if(!in_core && only_psf_core) continue;
	  
	  if(only_on_spots && spot_stamp_footprint(i,j)==0) continue;
	  
	  double val = spot->flux*psf->PSFValueWithParamsXY(spot->xc,spot->yc, i, j, spot_params, 0, 0, in_core, has_tail); // compute CPU expensive PSF core only if needed
	  
	  if(val>0 || (!only_positive))
	    model_image(i,j) += val; // this now includes core and tails
	  else if(model_image(i,j)==0) 
	    model_image(i,j) = eps;
	  
	} // end of loop on i  
      } // end of loop on j
    } // end of loop on spots
  } // end of loop on bundles
} // end of func
