#include <specex_model_image.h>
#include <specex_message.h>

using namespace std;


specex::Stamp specex::compute_stamp(const specex::image_data& model_image, const specex::PSF_p psf, const std::vector<specex::Spot_p>& spots) {
  
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

void specex::parallelized_compute_model_image(specex::image_data& model_image, const specex::image_data& weight, const specex::PSF_p psf, const std::vector<specex::Spot_p>& spots, const int first_fiber, const int last_fiber, bool only_on_spots, bool only_psf_core, bool only_positive, double minimal_tail_amp) {

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


  Stamp stamp = compute_stamp(model_image,psf,spots);
  int step_j  = (stamp.end_j-stamp.begin_j)/number_of_image_chuncks;
  

  vector<specex::image_data> image_chunks;
  
  for(int c=0;c<number_of_image_chuncks;c++) {
    //image_data toto; 
    //toto.resize(model_image.Nx(),model_image.Ny());
    image_chunks.push_back(image_data(model_image.Nx(),model_image.Ny()));
  }
  
  int chunk=0;
  
#pragma omp parallel for 
  for(chunk=0; chunk<number_of_image_chuncks; chunk++) {

    int begin_j = stamp.begin_j + chunk*step_j;
    int end_j   = stamp.begin_j + (chunk+1)*step_j;
    if(chunk==number_of_image_chuncks-1) end_j = stamp.end_j;
    
    compute_model_image(image_chunks[chunk],weight,psf,spots,first_fiber,last_fiber,only_on_spots,only_psf_core,only_positive,minimal_tail_amp,begin_j,end_j);
  }
  
  model_image.data.clear();
  for(chunk=0; chunk<number_of_image_chuncks; chunk++) {
    model_image.data += image_chunks[chunk].data;
  }
  
  
}



void specex::compute_model_image(specex::image_data& model_image, const specex::image_data& weight, const specex::PSF_p psf, const std::vector<specex::Spot_p>& spots, const int first_fiber, const int last_fiber, bool only_on_spots, bool only_psf_core, bool only_positive, double minimal_tail_amp, int predefined_begin_j, int predefined_end_j) {
  
  Stamp global_stamp = compute_stamp(model_image,psf,spots);
  
  model_image.data.clear();
  
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
  
  int begin_j = global_stamp.begin_j;
  int end_j   = global_stamp.end_j;
  
  if(predefined_begin_j>=0)
    begin_j = max(begin_j,predefined_begin_j);
  if(predefined_end_j>=0)
    end_j = min(end_j,predefined_end_j);
  			     
  
  double eps=1.e-20;

#ifdef CONTINUUM
  for (int j=begin_j; j <end_j; ++j) {  
    
    
    for(std::map<int,PSF_Params>::iterator bundle_it = psf->ParamsOfBundles.begin(); bundle_it != psf->ParamsOfBundles.end(); bundle_it++) {
      
      int begin_i = max(global_stamp.begin_i, int(floor(psf->GetTrace(first_fiber).X_vs_Y.Value(double(j))+0.5))-psf->hSizeX-1);
      int end_i   = min(global_stamp.end_i  , int(floor(psf->GetTrace(last_fiber).X_vs_Y.Value(double(j))+0.5))+psf->hSizeX+2);
      
      for(int fiber=first_fiber; fiber<=last_fiber; fiber++) {
	double x = psf->GetTrace(fiber).X_vs_Y.Value(double(j));
	double w = psf->GetTrace(fiber).W_vs_Y.Value(double(j));
	double continuum_flux = bundle_it->second.ContinuumPol.Value(w);
	double expfact_for_continuum=continuum_flux/(2*M_PI*square(bundle_it->second.continuum_sigma_x));
	if(expfact_for_continuum!=0) {
	  for (int i=begin_i ; i <end_i; ++i) {    

	    if(weight(i,j)<=0) continue;
	    if(only_on_spots && spot_stamp_footprint(i,j)==0) continue;
	    
	    double val = expfact_for_continuum*exp(-0.5*square((i-x)/bundle_it->second.continuum_sigma_x));
	    if(val>0 || (!only_positive))
	      model_image(i,j) += val;
	    else if(model_image(i,j)==0) 
	      model_image(i,j) = eps;
	      
	  } // end of loop on i
	} // expfact
      } // end of loop on fiber
    } // end of loop on bundles
  } // end of loop on j   
#endif  

#ifdef EXTERNAL_TAIL  
  int psf_tail_index = psf->ParamIndex("TAILAMP");
#endif

  for(size_t s=0;s<spots.size();s++) {
    specex::Spot_p spot = spots[s];
    harp::vector_double spot_params = psf->AllLocalParamsXW(spot->xc,spot->wavelength,spot->fiber_bundle);
    
#ifdef EXTERNAL_TAIL 
    if(minimal_tail_amp) {
      if(spot_params(psf_tail_index)<minimal_tail_amp) spot_params(psf_tail_index)=minimal_tail_amp;
    }
#endif   

    for (int j=begin_j; j <end_j; ++j) {  
      
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
	else if(model_image(i,j)==0) 
	  model_image(i,j) = eps;
	     
	
      }
    }
  }
}
