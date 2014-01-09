#include <assert.h>

#include "harp.hpp"

//#include "specex_base_analytic_psf.h"
#include "specex_psf_fitter.h"
#include "specex_brent.h"
#include "specex_legendre.h"
#include "specex_spot.h"
#include "specex_spot_array.h"
#include "specex_linalg.h"
#include "specex_image_data.h"
#include "specex_fits.h"
#include "specex_message.h"


using namespace std;

void specex::PSF_Fitter::SelectFiberBundle(int bundle) {
  std::map<int,PSF_Params>::iterator it = psf->ParamsOfBundles.find(bundle);
  if(it==psf->ParamsOfBundles.end()) SPECEX_ERROR("no such bundle #" << bundle_id);
  bundle_id = bundle;
  psf_global_params = & (it->second.Polynomials);
}

void specex::PSF_Fitter::SetStampLimitsFromPSF(specex::Stamp& stamp, const specex::PSF_p psf, const double &X, const double &Y) {
  psf->StampLimits(X,Y,stamp.begin_i,stamp.end_i,stamp.begin_j,stamp.end_j);
  // now check image bounds
  stamp.begin_i = max(0,stamp.begin_i);
  stamp.end_i   = min(stamp.Parent_n_cols(),stamp.end_i);
  stamp.begin_j = max(0,stamp.begin_j);
  stamp.end_j   = min(stamp.Parent_n_rows(),stamp.end_j);
}
 
void specex::PSF_Fitter::SetStampLimitsFromPSF(specex::Stamp& stamp, const specex::PSF_p psf, const double &xc_min, const double &xc_max, const double &yc_min, const double &yc_max) {
  int k,p;
  psf->StampLimits(xc_min,yc_min,stamp.begin_i,k,stamp.begin_j,p);
  psf->StampLimits(xc_max,yc_max,k,stamp.end_i,p,stamp.end_j);
  
  // now check image bounds
  stamp.begin_i = max(0,stamp.begin_i);
  stamp.end_i   = min(stamp.Parent_n_cols(),stamp.end_i);
  stamp.begin_j = max(0,stamp.begin_j);
  stamp.end_j   = min(stamp.Parent_n_rows(),stamp.end_j);
}

struct BrentBox { // handler to data and params needed to compute chi2 in brent routine
  specex::PSF_Fitter& fitter;
  harp::vector_double& delta_P;
  vector<specex::Spot_p>& spots;
  BrentBox(specex::PSF_Fitter& i_fitter,  harp::vector_double& i_delta_P, vector<specex::Spot_p>& i_spots) :
    fitter(i_fitter),
    delta_P(i_delta_P),
    spots(i_spots)
  {}
};



double compute_chi2_for_a_given_step(const double &current_step, BrentBox* bbox) {
  bbox->fitter.Params += current_step*bbox->delta_P;
  double chi2 = bbox->fitter.ComputeChi2AB(false);
  bbox->fitter.Params -= current_step*bbox->delta_P; // go back after testing
  if(bbox->fitter.verbose) {
    SPECEX_INFO("brent step=" << current_step << " chi2=" << chi2);
  }
  return chi2;
}



int specex::PSF_Fitter::NPar(int nspots) const {
  int npar = 0; // fluxes
  if(fit_psf) npar += psf->VaryingCoordNPar(bundle_id);
  if(fit_trace) npar += psf->TracesNPar();
  if(fit_flux) npar += nspots;
  if(fit_position) npar += 2*nspots;
  return npar;
}


double specex::PSF_Fitter::ComputeChi2AB(bool compute_ab) 
{
  
  
  // update spot_tmp_data (spots are called several times because we loop on pixels)
  for(size_t s=0;s<spot_tmp_data.size();s++) {
    
    specex::SpotTmpData &tmp = spot_tmp_data[s];
    
    if(fit_flux) tmp.flux = Params(tmp.flux_parameter_index);
    if(fit_trace) {
      tmp.x = specex::dot(ublas::project(Params,ublas::range(tmp.trace_x_parameter_index,tmp.trace_x_parameter_index+tmp.trace_x_monomials.size())),tmp.trace_x_monomials);
      tmp.y = specex::dot(ublas::project(Params,ublas::range(tmp.trace_y_parameter_index,tmp.trace_y_parameter_index+tmp.trace_y_monomials.size())),tmp.trace_y_monomials);
    }
    if(fit_position) {
      tmp.x = Params(tmp.x_parameter_index);
      tmp.y = Params(tmp.y_parameter_index);
    }
    if(fit_psf) {
      tmp.psf_params = psf->FixedCoordParamsXW(tmp.x,tmp.wavelength,bundle_id,Params);
    }
    if(fit_psf && ( fit_trace || fit_position ) && compute_ab) { // need to update at each step monomials
      int index=0;
      for(int p=0;p<npar_fixed_coord;p++) {
	harp::vector_double legendre_monomials_for_this_psf_parameter = (*psf_global_params)[p].Monomials(tmp.x,tmp.wavelength);
	size_t m_size = legendre_monomials_for_this_psf_parameter.size();
	ublas::project(tmp.psf_monomials,ublas::range(index,index+m_size))=legendre_monomials_for_this_psf_parameter;
	index += m_size;
      }
    }
  }
  
  
  //#define USE_SPARSE_VECTOR
#ifndef USE_SPARSE_VECTOR 
  harp::vector_double H;
  if(compute_ab) {
    H.resize(nparTot);
  }
#endif
  
  double chi2 = 0;
  
  harp::matrix_double &Anc  = const_cast<specex::PSF_Fitter*>(this)->A;
  harp::vector_double &Bnc  = const_cast<specex::PSF_Fitter*>(this)->B;
  harp::vector_double gradPar,gradPos;
  harp::vector_double *gradPar_pointer = 0;
  harp::vector_double *gradPos_pointer = 0;
  
  if(compute_ab) {
    Anc *= 0;
    Bnc *= 0;
    if(fit_psf) {gradPar.resize(npar_fixed_coord); gradPar_pointer = &gradPar;}// will remain zero if (!fit_psf)
    if(fit_position || fit_trace) {gradPos.resize(2); gradPos_pointer = &gradPos;}
  }
  
  bool use_footprint = (spot_tmp_data.size()>1 && footprint_weight.Nx()>0);
  
  for (int j=stamp.begin_j; j <stamp.end_j; ++j) {
  
    for (int i=stamp.begin_i ; i < stamp.end_i; ++i) {
      
      double w = 1;
      if(use_footprint) { 
	w=footprint_weight(i,j); 
      } else{
	w=weight(i,j);
      }
      if (w<=0) continue;
      
      double res = double(image(i,j));
      
#ifdef USE_SPARSE_VECTOR
      ublas::coordinate_vector<double> H;// slow , why??
      if(compute_ab) H.resize(nparTot,10); 
#else
      if(compute_ab)
	H *= 0;
#endif
      
      int nspots_in_pix = 0;
      for(size_t s=0;s<spot_tmp_data.size();s++) { // loop on tmp spots data
	
	specex::SpotTmpData &tmp = spot_tmp_data[s];
	
	if( ! tmp.stamp.Contains(i,j)) continue;
	
	nspots_in_pix++;
	
	double psfVal =  psf->PSFValueWithParamsXY(tmp.x,tmp.y, i, j, tmp.psf_params, gradPos_pointer, gradPar_pointer);
	
	if(psfVal==PSF_NAN_VALUE) SPECEX_ERROR("PSF value returns NAN");

	res -= tmp.flux*psfVal;
	
	if (compute_ab) {
	  
	  if(fit_psf) {
	    int index = 0;
	    for(int p=0;p<npar_fixed_coord;p++) {
	      size_t m_size = (*psf_global_params)[p].coeff.size();
	      //blas::axpy(tmp.flux*gradPar[p],ublas::project(tmp.psf_monomials,ublas::range(index,index+m_size)),ublas::project(H,ublas::range(index,index+m_size))); // doesnt compile
	      ublas::project(H,ublas::range(index,index+m_size)) += (tmp.flux*gradPar[p])*ublas::project(tmp.psf_monomials,ublas::range(index,index+m_size));
	      index += m_size;
	    }
	  }
	  if(fit_trace) {
	    ublas::project(H,ublas::range(tmp.trace_x_parameter_index,tmp.trace_x_parameter_index+tmp.trace_x_monomials.size())) 
	      += (gradPos(0) * tmp.flux)*tmp.trace_x_monomials;
	    ublas::project(H,ublas::range(tmp.trace_y_parameter_index,tmp.trace_y_parameter_index+tmp.trace_y_monomials.size())) 
	      += (gradPos(1) * tmp.flux)*tmp.trace_y_monomials;
	  }
	  
	  if(fit_flux) H[tmp.flux_parameter_index] += psfVal;
	  
	  if(fit_position) {
	    H[tmp.x_parameter_index] += gradPos(0) * tmp.flux;
	    H[tmp.y_parameter_index] += gradPos(1) * tmp.flux;
	  }







	} // end of test compute_ab
      } // end of loop on spots
      
      chi2 += w*res*res;

      /////////////////////////////////////////////////////
      // this is for debugging when developping the code
      //#define TESTING_H_VECTOR
#ifdef TESTING_H_VECTOR

      if(compute_ab && fit_psf && specex::dot(H,H)>0) {
	cout << "DEBUGGING TESTING_H_VECTOR" << endl;
	
	footprint_weight.data *= 0;
	footprint_weight(i,j) = w; // only this pixel
	
	int nok = 0;
	for (int j=stamp.begin_j; j <stamp.end_j; ++j) {
	  for (int i=stamp.begin_i ; i < stamp.end_i; ++i) {
	    double w = 1;
	    if(use_footprint) { 
	      w=footprint_weight(i,j); 
	    } else{
	      w=weight(i,j);
	    }
	    if (w<=0) continue;
	    nok++;
	  }
	}
	cout << "nok = " << nok << endl;
	
	double chi2_0 = ComputeChi2AB(false);
	
      
	double eps = 1.e-6;
	for(size_t p=0;p<Params.size();p++) {
	  Params(p) += 0.5*eps;
	  double chi2_plus = ComputeChi2AB(false);
	  Params(p) -= eps;
	  double chi2_minus = ComputeChi2AB(false);
	  
	  double chi2_der_numeric  = (chi2_plus-chi2_minus)/eps;
	  if(chi2_der_numeric) {
	    double chi2_der_analytic = -2*w*res*H(p); 
	  
	    cout << "DEBUGGING chi2 derivative p=" << p << " numeric=" <<  chi2_der_numeric << " analytic=" << chi2_der_analytic << " diff=" << chi2_der_numeric-chi2_der_analytic;
	    if(chi2_der_analytic!=0) cout << " ratio " << chi2_der_numeric/chi2_der_analytic -1; 
	    cout << endl;
	  }
	}
	exit(12);
      }
#endif
      /////////////////////////////////////////////////////


      if (compute_ab) { 	
	specex::syr(w,H,A); //  A += w*h*h.transposed();
	specex::axpy(w*res,H,B); // B += w*res*h;
      }
      
    } // end of loop on pix coord. i
  } // end of loop on pix coord. j
  
  
  // psf priors
#ifdef STORAGE
  if(fit_psf && !(psf->Priors.empty())) {

    for(int s=0;s<nspots;s++) {
      const harp::vector_double& SpotParams=psf_Params_of_spots[s];
      
      int index=0;
      for(int p=0;p<npar_fixed_coord;p++) {
	
	const double& par = SpotParams(p);
	
	int ncoef=(*psf_global_params)[p].coeff.size();
	std::map<int,Prior*>::const_iterator it = psf->Priors.find(p);
	
	if(it==psf->Priors.end()) {index += ncoef; continue;} // no prior for this psf parameter
	const Prior* prior = it->second;
	
	if (compute_ab) {
	  
	  const harp::vector_double& sM = psf_Monomials_of_spots[s];
	  
	  for (int c=0; c<ncoef; c++, index++) {
	    
	    const double& monomial_val = sM(index);
	    Bnc(index)       += monomial_val * prior->hdChi2dx(par);
	    Anc(index,index) += square(monomial_val) * prior->hd2Chi2dx2(par);
	    
	  }
	}
	//cout << "spot=" << s << " p=" << p << " chi2=" << prior->Chi2(par) << endl;
	chi2 += prior->Chi2(par);
      }
      
    }

  }
#endif
  
  return chi2;
}

   
static double sign(const double& a, const double& b) {
  if(b>0) return fabs(a);
  return -fabs(a);
}

void specex::PSF_Fitter::ComputeWeigthImage(vector<specex::Spot_p>& spots, int* npix) {
  /* 
     footprint weight includes :
     - weight set to zero for pixels unaffected by any spot (depends on flux and psf footprint)
     - weight from input image, just to set weight = 0 to bad pixels, we ignore variance that include signal because biasing
     - if include_signal_in_weighg, weight account for Poisson signal to noise.
  */
  
  if(spots.size()>1) {
    // compute psf footprint
    footprint_weight.resize(weight.Nx(),weight.Ny());
    footprint_weight.data *= 0;
    
    
    // if include_signal_in_weight, first store model in footprint_weight

    *npix = 0;
    if(include_signal_in_weight) {
      
      SPECEX_INFO(" specex::PSF_Fitter::FitSeveralSpots, computing image model for weights");
      
      for(size_t s=0;s<spots.size();s++) {
	specex::Spot_p& spot= spots[s];
	if(spot->flux<=0) continue;
	harp::vector_double psfParams = psf->FixedCoordParamsXW(spot->xc,spot->wavelength,bundle_id);
	
	Stamp spot_stamp(image);
	SetStampLimitsFromPSF(spot_stamp,psf,spot->xc,spot->yc);
	for(int j=spot_stamp.begin_j;j<spot_stamp.end_j;j++) {
	  for(int i=spot_stamp.begin_i;i<spot_stamp.end_i;i++) {
	    if(weight(i,j)==0) continue; // no need to compute anything
	    footprint_weight(i,j) += spot->flux*psf->PSFValueWithParamsXY(spot->xc,spot->yc,i,j,psfParams,0,0);
	  }
	}
      }
      cout << "for debugging, write image" << endl;
      specex::write_new_fits_image("debugging-model-for-weights.fits",footprint_weight);
      
      for(int j=0;j<footprint_weight.Ny();j++) {
	for(int i=0;i<footprint_weight.Nx();i++) {
	  double model_flux=footprint_weight(i,j);
	  if(model_flux!=0) { // has been computed, pixel participates to footprint
	    double var = square(readout_noise);
	    if(model_flux>0) { // else negative fluctuation
	      var += model_flux/gain;
	      var += square(psf_error*model_flux);
	    }
	    footprint_weight(i,j) = 1/var;
	    (*npix)++;
	  }
	}
      }
      
      


      
      //exit(12);
      
    }else{
      // only footprint
      for(size_t s=0;s<spots.size();s++) {
	specex::Spot_p& spot= spots[s];
	Stamp spot_stamp(image);
	SetStampLimitsFromPSF(spot_stamp,psf,spot->xc,spot->yc);
	for(int j=spot_stamp.begin_j;j<spot_stamp.end_j;j++)
	  for(int i=spot_stamp.begin_i;i<spot_stamp.end_i;i++) {
	    footprint_weight(i,j)=weight(i,j);
	    (*npix)++;
	  }
      }
    }
    // mask out some regions with mis-understood lines
    mask.ApplyMaskToImage(footprint_weight,*psf,0);
    
    
    if(include_signal_in_weight) 
    {
      cout << "for debugging, write image" << endl;
      write_new_fits_image("debugging-weights.fits",footprint_weight);
      //exit(12); // debug
    }
    
  } // end of test on spot size

}


bool specex::PSF_Fitter::FitSeveralSpots(vector<specex::Spot_p>& spots, double *chi2, int *n_pixels, int *n_iterations) {
  
  if(fit_trace && fit_position) {
    SPECEX_ERROR("specex::PSF_Fitter::FitSeveralSpots cannot fit traces and spot positions at the same time");
  }

  // chi2 , iterations ...
  // ----------------------------------------------------
  double chi2_memory_slot = 0;
  double *psfChi2 = &chi2_memory_slot;
  if(chi2) psfChi2 = chi2;
  int niter_memory_slot = 0;
  int *niter = &niter_memory_slot;
  if(n_iterations) niter = n_iterations;
  int npixels_memory_slot = 0;
  int *npix = &npixels_memory_slot;
  if(n_pixels) npix = n_pixels;
  *niter=0;
  int maxiter = 100; 
  double oldChi2=1e30;
  if(*psfChi2<=0) *psfChi2 = 1e20;
  // ----------------------------------------------------
  
  
  // definition of fitted region of image
  // ----------------------------------------------------
  double xc_min=1e20;
  double xc_max=-1e20;
  double yc_min=1e20;
  double yc_max=-1e20;
  for(size_t s=0;s<spots.size();++s) {
    const specex::Spot_p& spot = spots[s];
    if(spot->xc<xc_min) xc_min=spot->xc;
    if(spot->xc>xc_max) xc_max=spot->xc;
    if(spot->yc<yc_min) yc_min=spot->yc;
    if(spot->yc>yc_max) yc_max=spot->yc;
  }
  SetStampLimitsFromPSF(stamp,psf,xc_min,xc_max,yc_min,yc_max);
  
  if(spots.size()>1) {
    // remove boundaries of stamp to avoid contamination of spots not included in the fit
    stamp.begin_i += (psf->hSizeX-min(psf->hSizeX,5)); // leave 5 pix
    stamp.end_i   -= (psf->hSizeX-min(psf->hSizeX,5)); // leave 5 pix
    //stamp.begin_j += max(0,psf->hSizeY-2);
    //stamp.end_j   -= max(0,psf->hSizeY-2);
  }
  // ----------------------------------------------------
  
  int npar_psf = 0;
  if(fit_psf) npar_psf = psf->VaryingCoordNPar(bundle_id);
  int npar_trace = 0;
  if(fit_trace) npar_trace = psf->TracesNPar();
  
  npar_fixed_coord = psf->FixedCoordNPar();
  npar_varying_coord = psf->VaryingCoordNPar(bundle_id);
  nparTot  = NPar(spots.size());
  
  SPECEX_INFO("specex::PSF_Fitter::FitSeveralSpots npar_fixed_coord   = " << npar_fixed_coord);
  SPECEX_INFO("specex::PSF_Fitter::FitSeveralSpots npar_varying_coord = " << npar_varying_coord);

  
  if(specex_verbose()) {
    stringstream ss;
    ss << "specex::PSF_Fitter::FitSeveralSpots starting fitting ";
    if(fit_flux) ss << "+flux";
    if(fit_psf) ss << "+psf";
    if(fit_trace) ss << "+trace";
    if(fit_position) ss << "+position";
    if(fit_psf) ss << " npar_psf=" << npar_psf;
    if(fit_trace) ss << " npar_trace=" << npar_trace;
    specex_info(ss.str());
  }
  
  
  // allocation of parameters, setting Params 
  // and recording some indices
  // ----------------------------------------------------
  Params.resize(nparTot);
  map<int,int> tmp_trace_x_parameter;
  map<int,int> tmp_trace_y_parameter;
  int index_of_spot_parameters;
  
  {
    int index=0;
    if(fit_psf) {
      for(size_t p=0;p<psf_global_params->size();p++) {
	const harp::vector_double& coeff=(*psf_global_params)[p].coeff;
	size_t c_size = coeff.size();
	ublas::project(Params,ublas::range(index,index+c_size)) = coeff;
	index += c_size;
      }
    }
    
    if(fit_trace) {  
      for(std::map<int,specex::Trace>::const_iterator it=psf->FiberTraces.begin(); it!=psf->FiberTraces.end(); ++it) {

	tmp_trace_x_parameter[it->first] = index;

	{
	  const harp::vector_double& coeff = it->second.X_vs_W.coeff;
	  size_t c_size = coeff.size();
	  ublas::project(Params,ublas::range(index,index+c_size)) = coeff;
	  index += c_size;
	}
	
	tmp_trace_y_parameter[it->first] = index;
	
	{
	  const harp::vector_double& coeff = it->second.Y_vs_W.coeff;
	   size_t c_size = coeff.size();
	  ublas::project(Params,ublas::range(index,index+c_size)) = coeff;
	  index += c_size;
	}
      }
    }

    index_of_spot_parameters = index;
  }
  // ----------------------------------------------------
  

  
  ComputeWeigthImage(spots,npix);
  
  
 
  // load spot_tmp_data 
  spot_tmp_data.clear();
  
  {    

    int index = index_of_spot_parameters;

    for(size_t s=0;s<spots.size();s++) {
      const specex::Spot_p& spot=spots[s];
      
      SpotTmpData tmp;
      tmp.flux = spot->flux;
      //tmp.x    = spot->xc;
      //tmp.y    = spot->yc;
      tmp.x    = psf->Xccd(spot->fiber,spot->wavelength);
      tmp.y    = psf->Yccd(spot->fiber,spot->wavelength);
      tmp.wavelength    = spot->wavelength;

      // set parameters concerning spots
      if(fit_flux) {
	tmp.flux_parameter_index = index++; 
	Params(tmp.flux_parameter_index) = tmp.flux;
      }
      if(fit_position) {
	tmp.x_parameter_index = index++; 
	Params(tmp.x_parameter_index) = tmp.x;
	tmp.y_parameter_index = index++; 
	Params(tmp.y_parameter_index) = tmp.x;
      }
      
      // stamp
      tmp.stamp = Stamp(image);
      SetStampLimitsFromPSF(tmp.stamp,psf,tmp.x,tmp.y);
      tmp.stamp = tmp.stamp.Intersection(stamp);
      
      // psf parameters
      if(fit_psf) {
	tmp.psf_params = psf->FixedCoordParamsXW(tmp.x,tmp.wavelength,bundle_id,Params);
      }else{
	tmp.psf_params = psf->FixedCoordParamsXW(tmp.x,tmp.wavelength,bundle_id);
      }
      
      // psf parameters legendre monomials
      if(fit_psf) {
	tmp.psf_monomials.resize(npar_varying_coord);
	int index=0;
	for(int p=0;p<npar_fixed_coord;p++) {
	  harp::vector_double legendre_monomials_for_this_psf_parameter = (*psf_global_params)[p].Monomials(tmp.x,tmp.wavelength);
	  size_t m_size = legendre_monomials_for_this_psf_parameter.size();
	  ublas::project(tmp.psf_monomials,ublas::range(index,index+m_size))=legendre_monomials_for_this_psf_parameter;
	  index += m_size;
	}
      }

      // trace monomials
      if(fit_trace) {
	const specex::Trace& trace = psf->FiberTraces[spot->fiber];
	tmp.trace_x_monomials = trace.X_vs_W.Monomials(spot->wavelength);
	tmp.trace_y_monomials = trace.Y_vs_W.Monomials(spot->wavelength);
	tmp.trace_x_parameter_index = tmp_trace_x_parameter[spot->fiber];
	tmp.trace_y_parameter_index = tmp_trace_y_parameter[spot->fiber];
      }
      
      spot_tmp_data.push_back(tmp);
    }
  }
  
  if(verbose)
    SPECEX_INFO("specex::PSF_Fitter::FitSeveralSpots inc. signal in w=" << include_signal_in_weight << ", npix footprint = " << *npix);
  B.resize(nparTot);
  A.resize(nparTot, nparTot);
  
  
  
    
  while(true) { // minimization loop 
      
    oldChi2 = *psfChi2;
    if(verbose) SPECEX_INFO("specex::PSF_Fitter::FitSeveralSpots iter=" << *niter << " old chi2=" << oldChi2);
    
    if(verbose) SPECEX_INFO("specex::PSF_Fitter::FitSeveralSpots filling matrix ...");
    *psfChi2 = ComputeChi2AB(true);
    
    if(verbose) SPECEX_INFO("specex::PSF_Fitter::FitSeveralSpots solving ...");
    
    //harp::matrix_double As=A;
    //harp::vector_double Bs=B;
    int status = cholesky_solve(A,B);
    if (status != 0) {
      *psfChi2 = 1e30; 
      //As.writeFits("A.fits");
      //Bs.writeASCII("B.dat");
      if(fatal) {
	SPECEX_ERROR("cholesky_solve failed with status " << status);
      } else {
	SPECEX_WARNING("cholesky_solve failed with status " << status);
	return false;
      }
    }
    
    fitWeight = A; // to extract covariance once at minimum.
    
    bool linear = false;
    //#warning TEST_DE_JULIEN
     
    if( (fit_flux) && (!fit_position) && (!fit_trace) && (!fit_psf) ) {
      if( verbose) SPECEX_INFO("specex::PSF_Fitter::FitSeveralSpots linear because only fit flux");
      linear = true;
    }
    if( (fit_psf) && (!fit_position) && (!fit_trace) && (!fit_flux) && psf->IsLinear()) {
      if( verbose) SPECEX_INFO("specex::PSF_Fitter::FitSeveralSpots linear because only fit psf (linear wrt params)");
      linear = true;
    } 
    
    bool use_brent = true;
    if(1) {
      if(linear) {
	if( verbose) SPECEX_INFO("specex::PSF_Fitter::FitSeveralSpots no brent because linear");
	use_brent = false;
      }else{
	// check whether step indeed decreases chi2
	double chi2_0 = ComputeChi2AB(false);
	Params += B;
	double chi2_1 = ComputeChi2AB(false);
	Params -= B;
	if(chi2_1>chi2_0) {
	  use_brent = true;
	  if( verbose) SPECEX_INFO("specex::PSF_Fitter::FitSeveralSpots use brent because step 1 increases chi2");
	}else{
	  //#warning TEST_DE_JULIEN
	  use_brent = false;
	  if( verbose) SPECEX_INFO("specex::PSF_Fitter::FitSeveralSpots don't brent because chi2 decreases");	  
	}
      }
    }
    if(use_brent) {
      double brent_precision = chi2_precision/10.;
      if(brent_precision>0.1) brent_precision = 0.1;
      
      if(verbose) SPECEX_INFO("specex::PSF_Fitter::FitSeveralSpots starting brent with precision = " << brent_precision << " ...");
     
      // need to use brent here
      BrentBox bbox(*this,B,spots);
      
      double min_step=0;
      double max_step=2;
      double prefered_step=1;
      int status=0;
      
      // check the chi2 decrement is not good enought with step=1
      double best_chi2=compute_chi2_for_a_given_step(1,&bbox);      
      double best_step = 1;
      
      if(fabs(best_chi2-*psfChi2)>chi2_precision) { // really try brent now
	best_step = brent((AnalyticFunction*)(compute_chi2_for_a_given_step),
			  min_step,prefered_step,max_step,
			  brent_precision,&bbox,best_chi2,status,100);
      }else{
	if(verbose) SPECEX_INFO("brent not needed for required chi2 decrement dchi2=" << (-best_chi2+*psfChi2));
	if(best_chi2>*psfChi2) best_step=0; // check this anyway
      }
      
      if (status != 0) {
	SPECEX_WARNING("specex::PSF_Fitter::FitSeveralSpots brent fit didn't converge" 
		       << " best_step=" << best_step << " best_chi2=" << best_chi2);
	
	if(best_chi2>*psfChi2) {
	  if(fatal) {
	    SPECEX_ERROR("specex::PSF_Fitter::FitSeveralSpots brent doesn't improve things, this is a failure");
	  } else{
	    SPECEX_WARNING("specex::PSF_Fitter::FitSeveralSpots brent doesn't improve things, this is a failure");
	    return false;
	  }
	} else {
	  SPECEX_WARNING("specex::PSF_Fitter::FitSeveralSpots brent improves things so we continue");
	}
      }
      if(verbose) {
	double dchi2 = (-best_chi2+*psfChi2);
	if(dchi2>0 || fabs(dchi2)<brent_precision) {
	  if(verbose) SPECEX_INFO("specex::PSF_Fitter::FitSeveralSpots successful brent fit with step = " << best_step << " dchi2=" << dchi2
	       << " chi2pdf = " << best_chi2/(*npix-Params.size())
	      );
	}else{
	  SPECEX_WARNING("problem with brent dchi2 = " << dchi2);
	  best_step = 0;
	  best_chi2 = *psfChi2;
	}
      }
      
      Params += best_step*B;
      *psfChi2 = best_chi2;
	
    } else { // didn't use brent
      Params += B;
      *psfChi2 = ComputeChi2AB(false);
      if(verbose) SPECEX_INFO("specex::PSF_Fitter::FitSeveralSpots dchi2=" << oldChi2-*psfChi2);
      
    }
    
    // some sanity checks flux == nan is really bad
    for (unsigned k=0; k < nparTot; ++k) {
      if (isnan(Params(k))) {
	if(fatal) {
	  SPECEX_ERROR("specex::PSF_Fitter::FitSeveralSpots one parameter read nan");
	} else {
	  SPECEX_WARNING("specex::PSF_Fitter::FitSeveralSpots one parameter read nan");
	}
	*psfChi2 = 1e30;
	return false;
      }
    }
    
    (*niter) ++;
    
    // ending tests

    if( (*niter) > maxiter ) {
      SPECEX_WARNING("specex::PSF_Fitter::FitSeveralSpots quit loop because max iterations reached");
      break;
    }
    if( fabs(oldChi2-*psfChi2)<chi2_precision) {
      if( verbose) SPECEX_INFO("specex::PSF_Fitter::FitSeveralSpots quit loop because dchi2 small");
      break;
    }
    if( linear ) {
      if( verbose) SPECEX_INFO("specex::PSF_Fitter::FitSeveralSpots quit loop because linear");
      break;
    }
  } // end of minimization loop
    
  // We have to extract the weight matrix of the psf parameters (the first npar of params vector).
  //  This involves "marginalization" over the position and flux of the star. Compute 
  // full covariance matrix (from chi2 second derivatives), extract the psf params sub block,
  // and invert it back to get a weight matrix
  
  
  if (specex::cholesky_invert_after_decomposition(fitWeight) != 0) {
    SPECEX_ERROR("cholesky_invert_after_decomposition failed");
  }
  // fitWeight.Symmetrize("L");
  
  // fitWeight is now a covariance matrix
  // just to avoid confusion :
  harp::matrix_double& fitCovmat = fitWeight;

  int index=0;

  // copy fitted parameters in the right places:
  if (fit_psf) {
    // save params to psf
    
    for(size_t p=0;p<psf_global_params->size();p++) {
      harp::vector_double& coeff=(*psf_global_params)[p].coeff;
      for(size_t c=0;c<coeff.size();c++,index++)
	coeff(c)=Params(index);
    }
    
    
  }
  if (fit_trace) {
    for(std::map<int,specex::Trace>::iterator it=psf->FiberTraces.begin(); it!=psf->FiberTraces.end(); ++it) {
      harp::vector_double& PX = it->second.X_vs_W.coeff;
      for(size_t k=0;k<PX.size();k++,index++)
	PX(k)=Params(index);
      harp::vector_double& PY = it->second.Y_vs_W.coeff;
      for(size_t k=0;k<PY.size();k++,index++)
	PY(k)=Params(index);
    }
  }
  
  // save results for spots
  if(verbose) SPECEX_INFO("specex::PSF_Fitter::FitSeveralSpots saving spots fluxes");
  for(size_t s=0;s<spots.size();s++) {
    specex::Spot_p& spot= spots[s];
    specex::SpotTmpData& tmp = spot_tmp_data[s];
    if(fit_flux) {
      spot->flux = Params(tmp.flux_parameter_index); 
      spot->eflux = sqrt(fitCovmat(tmp.flux_parameter_index,tmp.flux_parameter_index));
    }
    if(fit_position) {
      spot->xc = Params(tmp.x_parameter_index);
      spot->yc = Params(tmp.y_parameter_index);
    }
    if(fit_trace) {
      spot->xc = psf->FiberTraces[spot->fiber].X_vs_W.Value(spot->wavelength);
      spot->yc = psf->FiberTraces[spot->fiber].Y_vs_W.Value(spot->wavelength);
    }
  }
  
  // spot->chi2=*psfChi2;
  
  bool ok=(*niter < maxiter);
  if(verbose) {
    
    if(ok) {
      if(specex_verbose()) {
	cout << "INFO specex::PSF_Fitter::FitSeveralSpots successful fit of ";
	if(fit_flux) cout << "+flux";
	if(fit_position) cout << "+pos";
	if(fit_psf) cout << "+psf";
	cout << " chi2= " << *psfChi2;
	cout << " niter=" << *niter;
	cout << endl;
      }
    }else{
      SPECEX_WARNING("specex::PSF_Fitter::FitSeveralSpots failed because reached max number of iterations");
    }

  }
  return ok;
};
 
bool specex::PSF_Fitter::FitOneSpot(specex::Spot_p& spot, double *chi2, int *n_iterations) {

  specex::Spot saved_spot = *spot;

  vector<specex::Spot_p> spots; 
  spots.push_back(spot);
  
  
  // int npar_psf = psf->FixedCoordNPar();
  // psf_global_params->clear();
  // for(int p=0;p<npar_psf;p++) {
  //   psf_global_params->push_back(specex::Legendre2DPol(0,0,image.Nx(),0,0,image.Ny()));
  //   (*psf_global_params)[p].coeff(0) = spot->PSFParams(p);
  // }
  
  
  bool ok = FitSeveralSpots(spots,&(spot->chi2),n_iterations);
  
  if(verbose) {
    
    if(ok) {
      cout << "INFO specex::PSF_Fitter::FitOneSpot successful fit of flux";
      if(fit_position) cout << "+pos";
      if(fit_psf) cout << "+psf";
      cout << " chi2= " << spot->chi2;
      cout << " dchi2= " <<  saved_spot.chi2-spot->chi2;
      if(saved_spot.flux>0)
	 cout << " dflux/flux=" << spot->flux/saved_spot.flux-1;
      cout << " dx=" << spot->xc-saved_spot.xc;
      cout << " dy=" << spot->yc-saved_spot.yc;
      if(n_iterations)
      cout << " niter=" << n_iterations;
      cout << endl;
    }else{
      cout << "WARNING specex::PSF_Fitter::FitOneSpot failed because reached max number of iterations" << endl;
    } 
  }
  return ok;
};

void specex::PSF_Fitter::SetPSFParams(const harp::vector_double &ParamsToSet) {
  // also, define psf polynomes (should be in fitter)
  psf_global_params->clear();
  for (size_t k =0; k < ParamsToSet.size(); ++k) {
    psf_global_params->push_back(specex::Legendre2DPol(0,0,image.Nx(),0,0,image.Ny()));
    psf_global_params->back().coeff(0)=ParamsToSet(k);
  }
}



bool specex::PSF_Fitter::FitTraces(vector<specex::Spot_p>& spots, int *n_fibers_fitted) {
  
  SPECEX_INFO("specex::PSF_Fitter::FitTraces starting");

  int nok_memory_slot;
  int *nok = &nok_memory_slot;
  if(n_fibers_fitted) nok = n_fibers_fitted;
  
  *nok = 0;
  for(map<int,specex::Trace>::iterator it=psf->FiberTraces.begin();
      it !=psf->FiberTraces.end(); ++it) {
    
    specex::Trace& trace=it->second;
    
    vector<specex::Spot_p> spots_of_fiber;
    for(size_t s=0;s<spots.size();s++) {
      specex::Spot_p& spot = spots[s];
      if(spot->fiber==trace.fiber)
	spots_of_fiber.push_back(spot);
    }
    bool ok = trace.Fit(spots_of_fiber);
    if(ok) { 
      (*nok)++;
      
      for(size_t s=0;s<spots_of_fiber.size();s++) {
	specex::Spot_p& spot = spots_of_fiber[s];
	spot->xc=trace.X_vs_W.Value(spot->wavelength);
	spot->yc=trace.Y_vs_W.Value(spot->wavelength);

      }


    }
  }
  SPECEX_INFO("specex::PSF_Fitter::FitTraces ended nok = " << (*nok) << "/" << psf->FiberTraces.size());
 
  return true;
}

bool specex::PSF_Fitter::FitIndividualSpotFluxes(std::vector<specex::Spot_p>& spots) {
  
  SPECEX_INFO("fitting independently the flux of each spot");
  include_signal_in_weight = false;
  fit_flux                 = true;
  fit_position             = false;
  fit_psf                  = false;
  fit_trace                = false;
  verbose                  = false; // debug
  fatal                    = false; // debug
    
  int nok=0;
  for(size_t s=0;s<spots.size();s++) {
    specex::Spot_p& spot = spots[s];    
    spot->eflux = 0;
    spot->flux = 0;	
    spot->status=-1;
    
    bool ok = FitOneSpot(spot);
    
    if(ok) {
      nok++;
      spot->status=1;
    }else{
      spot->status=0;
      SPECEX_WARNING("fit of flux of spot at x = " << spot->xc << " " << spot->yc << " failed");
    }
    if(int(s)%100==0 && s!=0) SPECEX_INFO("done " << s << "/" << spots.size() << " ...");
  }
  SPECEX_INFO("successful fit of each spot flux for " << nok << "/" << spots.size());
  return true;
}
bool specex::PSF_Fitter::FitIndividualSpotPositions(std::vector<specex::Spot_p>& spots) {
  
  SPECEX_INFO("fitting independently the flux+position of each spot");
  include_signal_in_weight = false;
  fit_flux                 = true;
  fit_position             = true;
  fit_psf                  = false;
  fit_trace                = false;
  verbose                  = false;
  fatal                    = false;
  
  int nok=0;
  for(size_t s=0;s<spots.size();s++) {
    specex::Spot_p& spot = spots[s];    
    spot->initial_xc = spot->xc;
    spot->initial_yc = spot->yc;
    spot->eflux = 0;
    spot->flux = 0;	
    spot->status=-1;
    
    bool ok = FitOneSpot(spot);
    
    if(ok) {
      nok++;
      spot->status=1;
    }else{
      spot->status=0;
    }
    if(int(s)%100==0 && s!=0) SPECEX_INFO("done " << s << "/" << spots.size() << " ...");
  }
  SPECEX_INFO("successful fit of each spot flux+pos for " << nok << "/" << spots.size());
  return true;
}

bool specex::PSF_Fitter::FitEverything(std::vector<specex::Spot_p>& input_spots, bool init_psf) {

  if(input_spots.size()==0) {
    SPECEX_ERROR("in fit_several_spots spotarrays is empty");
  }

  
  
  
  
  verbose = true;
  
  double chi2=1e30;
  int npix = 0;
  int niter = 0;
    
  SPECEX_INFO("starting to fit PSF with " <<  input_spots.size() << " spots");
 
  if(init_psf) {
    
    SPECEX_INFO("init PSF");
    
    FitTraces(input_spots);

    
    // here we need to count number of wavelength and spots per wavelength
    {
      std::vector<SpotArray> spot_arrays = find_spot_arrays(input_spots);
      int nspots_per_wave=0;
      
      double min_wave = 1e12;
      double max_wave = -1e12;
      double min_x = 1e12;
      double max_x = -1e12;
      
      for(size_t s=0;s<input_spots.size(); ++s) {
	const specex::Spot_p spot = input_spots[s];
	if(spot->xc < min_x ) min_x = spot->xc;
	if(spot->xc > max_x ) max_x = spot->xc;
	if(spot->wavelength < min_wave ) min_wave = spot->wavelength;
	if(spot->wavelength > max_wave ) max_wave = spot->wavelength;
	
      }

      
      for(size_t a=0;a<spot_arrays.size();a++) {
	int asize = spot_arrays[a].size();
	if(asize>nspots_per_wave) 
	  nspots_per_wave=asize;
      }
      if(polynomial_degree_along_x>nspots_per_wave-1) {
	polynomial_degree_along_x=nspots_per_wave-1;
	SPECEX_WARNING("Reducing polynomial degree along x to " << polynomial_degree_along_x << " because of number of fibers");
      }
      if(polynomial_degree_along_wave>int(spot_arrays.size())-1) {
	polynomial_degree_along_wave=int(spot_arrays.size())-1;
	SPECEX_WARNING("Reducing polynomial degree along wavelength to " << polynomial_degree_along_wave << " because of number of spots");
      }
      
      
      SPECEX_INFO("Setting PSF polynomial degrees " << polynomial_degree_along_x << " " << polynomial_degree_along_wave);
      
      psf_global_params->clear();

      int npar = psf->NPar();
      for(int p=0;p<npar;p++) {
	//if(p<first_index_for_tail)
	psf_global_params->push_back(specex::Legendre2DPol(polynomial_degree_along_x,min_x,max_x,polynomial_degree_along_wave,min_wave,max_wave));
	//else
	//psf_global_params->push_back(specex::Legendre2DPol(polynomial_degree_along_x_for_psf_tail,0,image.Nx(),polynomial_degree_along_y_for_psf_tail,0,image.Ny()));
      }
      
      harp::vector_double default_params = psf->DefaultParams();
      for(int p=0;p<npar;p++) {
	(*psf_global_params)[p].coeff *= 0;
	(*psf_global_params)[p].coeff[0] = default_params[0];
      }
      
    }

    

    /*
    bool ok = InterpolateSpotPSFs(input_spots,&chi2,&niter);
    if(ok) {
      SPECEX_INFO("InterpolateSpotPSFs successful");
      for(size_t s=0;s<input_spots.size();s++)
	input_spots[s]->status++;
    }
    else {
      SPECEX_ERROR("InterpolateSpotPSFs failed");
    }
    */

  }
  
  fatal = true;
  include_signal_in_weight = false;
  chi2_precision = 10;
  bool ok = true;
  

  if(psf->Name()=="NONE") {
    
    psf->ClearPriors();
    // from a not yet perfect fit of laser data in r1 CCD at lambda=8702 , expnum=00120726
    cout << "INFO Setting priors to PSF tails" << endl;
    cout << "DEBUG tail indices :" << endl;
    cout << "DEBUG TailNorm : " << psf->ParamIndex("TailNorm") << endl;
    cout << "DEBUG TailPower : " << psf->ParamIndex("TailPower") << endl;
    cout << "DEBUG TailXposScale : " << psf->ParamIndex("TailXposScale") << endl;
    cout << "DEBUG TailXnegScale : " << psf->ParamIndex("TailXnegScale") << endl;
    cout << "DEBUG TailYnegScale : " << psf->ParamIndex("TailYnegScale") << endl;
    cout << "DEBUG Npar : " << psf->NPar() << endl;
    
    double epsilon = 1.e-12;
    if(psf->HasParam("YTailNorm")) psf->SetPrior("YTailNorm",new GaussianPrior(0,epsilon)); // 
    if(psf->HasParam("TailNorm")) psf->SetPrior("TailNorm",new GaussianPrior(0,epsilon)); // 
    if(psf->HasParam("TailPower")) psf->SetPrior("TailPower",new GaussianPrior(2.,epsilon)); // 
    if(psf->HasParam("TailXposScale")) psf->SetPrior("TailXposScale",new GaussianPrior(1,epsilon));
    if(psf->HasParam("TailXnegScale")) psf->SetPrior("TailXnegScale",new GaussianPrior(1,epsilon));
    if(psf->HasParam("TailYnegScale")) psf->SetPrior("TailYnegScale",new GaussianPrior(1,epsilon));
    
    //exit(12);
  }

  
  
  SPECEX_INFO("Starting FitSeveralSpots FLUX");
  SPECEX_INFO("=======================================");
  fit_flux       = true;
  fit_position   = false;
  fit_psf        = false;
  fit_trace      = false;
  ok = FitSeveralSpots(input_spots,&chi2,&npix,&niter);
  if(!ok) SPECEX_ERROR("FitSeveralSpots failed for FLUX");
  
  
// selecting spots
  double minimum_signal_to_noise = 3;
  std::vector<specex::Spot_p> selected_spots;
  for(size_t s=0;s<input_spots.size();s++) {
    specex::Spot_p& spot = input_spots[s];
    if(spot->flux>0  && spot->eflux>0 && spot->flux/spot->eflux>minimum_signal_to_noise) {
      spot->status++;
      selected_spots.push_back(spot);
    }
    
  }
  SPECEX_INFO("selected " << selected_spots.size() << " spots out of " << input_spots.size() << " with S/N>" << minimum_signal_to_noise);

  
  // SPECEX_INFO("STOP HERE FOR DEBUG"); return ok;

  SPECEX_INFO("Starting FitSeveralSpots PSF");
  SPECEX_INFO("=======================================");
  fit_flux       = false;
  fit_position   = false;
  fit_psf        = true;
  fit_trace      = false;
  ok = FitSeveralSpots(selected_spots,&chi2,&npix,&niter);
  if(!ok) SPECEX_ERROR("FitSeveralSpots failed for PSF");
  
  // write_spots_list(spots,"spots-tmp-psf->list",PSF);
  
  // SPECEX_INFO("STOP HERE FOR DEBUGGING"); return ok;


  SPECEX_INFO("Starting FitSeveralSpots PSF+FLUX");
  SPECEX_INFO("=======================================");
  fit_flux       = true;
  fit_position   = false;
  fit_psf        = true;
  fit_trace      = false;
  ok = FitSeveralSpots(selected_spots,&chi2,&npix,&niter);
  if(!ok) SPECEX_ERROR("FitSeveralSpots failed for PSF+FLUX");

  SPECEX_INFO("Starting FitSeveralSpots TRACE");
  SPECEX_INFO("=======================================");
  fit_flux       = false;
  fit_position   = false;
  fit_psf        = false;
  fit_trace      = true;
  ok = FitSeveralSpots(selected_spots,&chi2,&npix,&niter);
  if(!ok) SPECEX_ERROR("FitSeveralSpots failed for TRACE");
  
  SPECEX_INFO("Starting FitSeveralSpots TRACE+FLUX");
  SPECEX_INFO("=======================================");
  fit_flux       = true;
  fit_position   = false;
  fit_psf        = false;
  fit_trace      = true;
  ok = FitSeveralSpots(selected_spots,&chi2,&npix,&niter);
  if(!ok) SPECEX_ERROR("FitSeveralSpots failed for TRACE+FLUX");
  
  SPECEX_INFO("Starting FitSeveralSpots PSF+FLUX");
  SPECEX_INFO("=======================================");
  fit_flux       = true;
  fit_position   = false;
  fit_psf        = true;
  fit_trace      = false;
  chi2_precision = 0.1;
  include_signal_in_weight = true;

  ok = FitSeveralSpots(selected_spots,&chi2,&npix,&niter);
  if(!ok) SPECEX_ERROR("FitSeveralSpots failed for PSF+FLUX");
  
  // write_spots_list(spots,"spots-tmp-PSF+FLUX.list",PSF);
  // psf->write("psf-tmp-PSF+FLUX.dat");
  
  
  SPECEX_INFO("Starting FitSeveralSpots FLUX (all spots)");
  SPECEX_INFO("=======================================");
  fit_flux       = true;
  fit_position   = false;
  fit_psf        = false;
  fit_trace      = false;
  ok = FitSeveralSpots(input_spots,&chi2,&npix,&niter);
  if(!ok) SPECEX_ERROR("FitSeveralSpots failed for FLUX");
  
  return ok;

#ifdef LA_SUITE

  cout << "==== Starting FitSeveralSpots TRACE+FLUX+PSF ==== " << endl;
  fit_flux       = true;
  fit_position   = false;
  fit_psf        = true;
  fit_trace      = true;
  chi2_precision = 10;
  ok = FitSeveralSpots(selected_spots,&chi2,&npix,&niter);
  if(!ok) {
    cout << "FitSeveralSpots failed for TRACE" << endl;
    return false;
  }
  
  write_spots_list(selected_spots,"spots-TRACE+FLUX+psf->list",PSF);
  write_spots_data(selected_spots,"spots-TRACE+FLUX+psf->dat");
  psf->write("psf-TRACE+FLUX+psf->dat");
  
  if(fit_psf_tails) {
    if(psf->Name()=="GAUSSHERMITE") {
      cout << "==== Starting FitSeveralSpots PSF TAILS ==== " << endl;
      psf->ClearPriors();
      
      
      // from a not yet perfect fit of laser data in r1 CCD at lambda=8702 , expnum=00120726
      double epsilon=1.e-12;
      //if(psf->HasParam("TailNorm")) psf->SetPrior("TailNorm",new GaussianPrior(1.8,epsilon)); // slope
      if(psf->HasParam("TailPower")) psf->SetPrior("TailPower",new GaussianPrior(2.,epsilon)); // slope
      if(psf->HasParam("TailXposScale")) psf->SetPrior("TailXposScale",new GaussianPrior(1.2,epsilon));
      if(psf->HasParam("TailXnegScale")) psf->SetPrior("TailXnegScale",new GaussianPrior(1.2,epsilon));
      if(psf->HasParam("TailYnegScale")) psf->SetPrior("TailYnegScale",new GaussianPrior(0.9,epsilon));
      
      include_signal_in_weight = false; // to deweight the core for the determination of tails
      fit_flux       = false;
      fit_position   = false;
      fit_psf        = true;
      fit_trace      = false;
      chi2_precision = 10;
      ok = FitSeveralSpots(selected_spots,&chi2,&npix,&niter);
      if(!ok) {
	cout << "FitSeveralSpots failed for PSF+FLUX" << endl;
	return false;
      }
      write_spots_data(selected_spots,"spots-PSF+TAILS.dat");
      write_spots_list(selected_spots,"spots-PSF+TAILS.list",PSF);
      psf->write("psf-PSF+TAILS.dat");
      
      /*
	include_signal_in_weight = true;
	cout << "==== Starting FitSeveralSpots PSF TAILS with signal in weights ==== " << endl;
	ok = FitSeveralSpots(spots,&chi2,&npix,&niter);
	if(!ok) {
	cout << "FitSeveralSpots failed for PSF+FLUX" << endl;
	return false;
	}
    */
      
    }
    
    
  
    
    //include_signal_in_weight = true;
    fit_flux       = true;
    fit_position   = false;
    fit_psf        = true;
    fit_trace      = false;
    chi2_precision = 1;
    
    for(int loop=0; loop<1; loop++) { // because change of weights
      cout << "==== Starting FitSeveralSpots PSF+TAILS+FLUX+TRACE"<<endl;
      ok = FitSeveralSpots(selected_spots,&chi2,&npix,&niter);
      if(!ok) {
	cout << "FitSeveralSpots failed for PSF+TAILS+FLUX+TRACE" << endl;
	return false;
      }
    }
    write_spots_data(selected_spots,"spots-PSF+TAILS+FLUX+TRACE.dat");
    write_spots_list(selected_spots,"spots-PSF+TAILS+FLUX+TRACE.list",PSF);
    psf->write("psf-PSF+TAILS+FLUX+TRACE.dat");
  }// end of test fit_psf_tails

#endif
  
  return ok;
}
