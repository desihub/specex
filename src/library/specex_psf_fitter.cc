#include <assert.h>
#include <time.h>

#include <harp.hpp>

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
#include "specex_psf_io.h"
#include "specex_model_image.h"
#include "specex_psf.h"

//#include <boost/archive/xml_oarchive.hpp>

#define SIDE_BAND_WEIGHT_SCALE 10.

using namespace std;

void specex::PSF_Fitter::SelectFiberBundle(int bundle) {
  std::map<int,PSF_Params>::iterator it = psf->ParamsOfBundles.find(bundle);
  if(it==psf->ParamsOfBundles.end()) SPECEX_ERROR("no such bundle #" << bundle);
  psf_params = & (it->second);
}

void specex::PSF_Fitter::SetStampLimitsFromPSF(specex::Stamp& stamp, const specex::PSF_p psf, const double &X, const double &Y) {
  stamp = Stamp(*(stamp.parent_image)); // reset
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
  double chi2 = 0;
  if(bbox->fitter.parallelized)
    chi2 = bbox->fitter.ParallelizedComputeChi2AB(false);
  else
    chi2 = bbox->fitter.ComputeChi2AB(false);

  bbox->fitter.Params -= current_step*bbox->delta_P; // go back after testing
  SPECEX_DEBUG("brent step=" << current_step << " chi2=" << chi2);
  
  return chi2;
}



int specex::PSF_Fitter::NPar(int nspots) const {
  int npar = 0; // fluxes
  if(fit_psf || fit_psf_tail) npar += psf->BundleNFitPar(psf_params->bundle_id);
  if(fit_trace) {
    //npar += psf->TracesNPar();
    for(std::map<int,specex::Trace>::const_iterator it=psf->FiberTraces.begin(); it!=psf->FiberTraces.end(); ++it) {
      if(it->first < psf_params->fiber_min || it->first > psf_params->fiber_max) continue;
      if(it->second.Off()) continue;
      
      npar += it->second.X_vs_W.coeff.size();
      npar += it->second.Y_vs_W.coeff.size();
      
    }
  }
  if(fit_flux) npar += nspots;
  if(fit_position) npar += 2*nspots;
#ifdef CONTINUUM
  if(fit_continuum) npar += psf_params->ContinuumPol.coeff.size();
#endif
  return npar;
}



double specex::PSF_Fitter::ParallelizedComputeChi2AB(bool compute_ab) {
  
  
  //int step_j  = (stamp.end_j-stamp.begin_j)/number_of_image_bands;
 
  
  //SPECEX_INFO("Begin parallelized ComputeChi2AB j range " << stamp.begin_j << " " << stamp.end_j);
  
  UpdateTmpData(compute_ab);
#ifdef EXTERNAL_TAIL  
  // precompute tail profile
  specex::SpotTmpData &tmp = spot_tmp_data[0];
  psf->TailProfile(0,0,psf->AllLocalParamsFW(tmp.fiber,tmp.wavelength,tmp.fiber_bundle));
#endif
  
  
  harp::vector_double chi2_of_band(number_of_image_bands);
  chi2_of_band.clear();
  
  // choosing bins
  harp::vector_double begin_j(number_of_image_bands);
  harp::vector_double end_j(number_of_image_bands);
  
  int band = 0;
  for(band=0; band<number_of_image_bands; band++) {
    begin_j(band)=0;
    end_j(band)=0;
  }
  

  band = 0;
  {
    int ref_fiber=(psf_params->fiber_min + psf_params->fiber_max)/2;
    while(psf->FiberTraces.find(ref_fiber)==psf->FiberTraces.end()) {ref_fiber++;} // we check we don't choose a missing one
    
    vector<int> spots_j;
    for(size_t s=0;s<spot_tmp_data.size();s++) {
      const specex::SpotTmpData &tmp = spot_tmp_data[s];
      if(tmp.fiber != ref_fiber) continue;
      if((tmp.y<0) || (tmp.y>=image.n_rows()))
	SPECEX_ERROR("Spot is outside of image ??? y=" << tmp.y);
      spots_j.push_back(int(tmp.y));
    }
    if(spots_j.empty()) {
      SPECEX_ERROR("No spots for ref fiber");
    }
    std::sort(spots_j.begin(),spots_j.end());
    
    //cout << "sorted spots yccd = ";
    //for(vector<int>::const_iterator j=spots_j.begin();j!=spots_j.end();j++) {
    //cout << " " << *j;
    //}

    //SPECEX_INFO("stamp begin " << stamp.begin_j << " end " << stamp.end_j);

    band = 0;
    int nspots_per_band = spots_j.size()/number_of_image_bands;
    vector<int>::const_iterator j=spots_j.begin();
    begin_j(0) = stamp.begin_j;
    int nspots=0;
    for(vector<int>::const_iterator j=spots_j.begin();j!=spots_j.end();j++) {
      if(band>0 && nspots==0) 
	begin_j(band) = end_j(band-1);
      
      end_j(band) = min(stamp.end_j,(*j)+psf->hSizeY+2);      
      nspots++;
      if(nspots>=nspots_per_band) {	
	band++;
	nspots=0;
	if(band>=number_of_image_bands) break;	
      }
    }
    
    // set last to bound of stamp
    for(band=number_of_image_bands-1;band>=0;band--)
      if(end_j(band)!=0) {end_j(band)=stamp.end_j; break;}
    
    
    //for(band=0; band<number_of_image_bands; band++) 
    //SPECEX_INFO("band " << band << " [" << begin_j(band) << ":" << end_j(band) << "]");
    
    
  }

  
  
#pragma omp parallel for 
  for(band=0; band<number_of_image_bands; band++) {
    if(end_j(band)>begin_j(band)) {
       chi2_of_band(band) = ComputeChi2AB(compute_ab,begin_j(band),end_j(band),& A_of_band[band], & B_of_band[band],false);
    }else if(compute_ab) {
      A_of_band[band].clear();
      B_of_band[band].clear();
    }
  }
  
  
  for(int band=1; band<number_of_image_bands; band++) {
    chi2_of_band(0) += chi2_of_band(band);
  }
  if(compute_ab) {
    for(int band=1; band<number_of_image_bands; band++) {
      A_of_band[0] += A_of_band[band];
      B_of_band[0] += B_of_band[band];
    }
  }
  //SPECEX_INFO("End of parallelized ComputeChi2AB chi2 = " << chi2_of_band(0));
  return chi2_of_band(0);
}

void specex::PSF_Fitter::InitTmpData(const vector<specex::Spot_p>& spots) {

  SPECEX_DEBUG("InitTmpData with " << spots.size() << " spots");
  
  // load spot_tmp_data 
  spot_tmp_data.clear();

  for(size_t s=0;s<spots.size();s++) {
    const specex::Spot_p spot=spots[s];
    
    SpotTmpData tmp;
    tmp.ignore = false;
    tmp.flux = spot->flux;
    tmp.x    = psf->Xccd(spot->fiber,spot->wavelength);
    tmp.y    = psf->Yccd(spot->fiber,spot->wavelength);

    // that's possible and not a bug 
    /*
    if ((tmp.y<0 || tmp.y>image.n_rows()) || (tmp.x<0 || tmp.x>image.n_cols())) {
      SPECEX_ERROR("spot outside image ?" 
		   << " t.x=" << tmp.x << " t.y=" << tmp.y 
		   << " s.x=" << spot->xc << " s.y=" << spot->yc
		   << " s.xi=" << spot->initial_xc << " s.yi=" << spot->initial_yc
		   << " s.fiber=" << spot->fiber << " s.wavelength=" << spot->wavelength << " s.flux=" << spot->flux);
    }
    */
    

    tmp.wavelength     = spot->wavelength;
    tmp.fiber          = spot->fiber;
    tmp.fiber_bundle   = spot->fiber_bundle;
        
    if( (fit_psf || fit_psf_tail) && tmp.flux<0) tmp.flux = 0.; // more robust
    tmp.frozen_flux = tmp.flux;
    
    



    // stamp
    tmp.stamp = Stamp(image);
    SetStampLimitsFromPSF(tmp.stamp,psf,tmp.x,tmp.y);
    tmp.stamp = tmp.stamp.Intersection(stamp);
    
    tmp.can_measure_flux = true;
    if(spots.size()>1) {
      // check whether we can rely on the flux measurement of this spot
      int i_center = int(floor(tmp.x)+0.5);
      int j_center = int(floor(tmp.y)+0.5);
      int nbad=0;
      for(int j=max(tmp.stamp.begin_j,j_center-2);j<min(tmp.stamp.end_j,j_center+3);j++) {
	for(int i=max(tmp.stamp.begin_i,i_center-2);i<min(tmp.stamp.end_i,i_center+3);i++) {
	  if(weight(i,j)==0) nbad++;
	}
      }
      tmp.can_measure_flux = (nbad<=5); // can survive one dead column = 5pix, not more
      if(!tmp.can_measure_flux) tmp.ignore = true; // ignore by default
      if(!tmp.can_measure_flux) SPECEX_WARNING("cannot measure flux of spot " << s << " at x=" << tmp.x << " y=" << tmp.y << ", will attach it to a neighbour");
    }
    
    // psf parameters
    tmp.psf_all_params = psf->AllLocalParamsXW(tmp.x,tmp.wavelength,psf_params->bundle_id);
    
    spot_tmp_data.push_back(tmp);
  }

  
#ifdef EXTERNAL_TAIL
  // compute this before parallel computing
  psf->TailProfile(0,0,psf->AllLocalParamsFW(spots[0]->fiber,spots[0]->wavelength,spots[0]->fiber_bundle));
#endif

}

void specex::PSF_Fitter::UpdateTmpData(bool compute_ab) {


#ifdef CONTINUUM
  if(fit_continuum)
    ublas::noalias(psf_params->ContinuumPol.coeff) = ublas::project(Params,ublas::range(continuum_index,continuum_index+psf_params->ContinuumPol.coeff.size()));
#endif
  
  // update spot_tmp_data (spots are called several times because we loop on pixels)
  for(size_t s=0;s<spot_tmp_data.size();s++) {
    
    specex::SpotTmpData &tmp = spot_tmp_data[s];
    
    if(fit_flux) {
      if(force_positive_flux) {
	tmp.flux = exp(min(max(Params(tmp.flux_parameter_index),-30.),+30.));
      }else{
	tmp.flux = Params(tmp.flux_parameter_index);
	if(tmp.flux<-100) {
	  SPECEX_WARNING("Ng. flux fiber=" << tmp.fiber << " wave=" << tmp.wavelength << " x=" << tmp.x << " y=" << tmp.y << " flux=" << tmp.flux);
	}
      }
    }
    if(fit_trace) {
      tmp.x = specex::dot(ublas::project(Params,ublas::range(tmp.trace_x_parameter_index,tmp.trace_x_parameter_index+tmp.trace_x_monomials.size())),tmp.trace_x_monomials);
      tmp.y = specex::dot(ublas::project(Params,ublas::range(tmp.trace_y_parameter_index,tmp.trace_y_parameter_index+tmp.trace_y_monomials.size())),tmp.trace_y_monomials);
    }
    if(fit_position) {
      tmp.x = Params(tmp.x_parameter_index);
      tmp.y = Params(tmp.y_parameter_index);
    }
    if(fit_psf || fit_psf_tail) {
      //harp::vector_double toto = tmp.psf_all_params;
      tmp.psf_all_params = psf->AllLocalParamsXW_with_FitBundleParams(tmp.x,tmp.wavelength,psf_params->bundle_id,Params);
      /*
      // or , test
      if(npar_fixed_coord>2 && tmp.psf_all_params(2)!=0) {
	
	for(int i=0;i<2;i++) {
	  cout << "TEST all param " << i << " " << tmp.psf_all_params(i) << " " << toto(i) << " " << tmp.psf_all_params(i)-toto(i) << endl;
	}
	int index = 0;
	for(int p=0;p<npar_fixed_coord;p++) {
	  int c_size = psf_params->FitParPolXW[p]->coeff.size();
	  double pval = specex::dot(ublas::project(tmp.psf_monomials,ublas::range(index,index+c_size)),ublas::project(Params,ublas::range(index,index+c_size)));
	  cout << "TEST param " << p << " " << pval << " " << tmp.psf_all_params(p+2) << " diff=" << pval-tmp.psf_all_params(p+2) << endl;
	  index += c_size;
	}
	exit(12);
      }
      */
    }
    if(fit_psf && ( fit_trace || fit_position ) && compute_ab) { // need to update at each step monomials
      int index=0;
      for(int p=0;p<npar_fixed_coord;p++) {
	harp::vector_double legendre_monomials_for_this_psf_parameter = psf_params->FitParPolXW[p]->Monomials(tmp.x,tmp.wavelength);
	size_t m_size = legendre_monomials_for_this_psf_parameter.size();
	ublas::noalias(ublas::project(tmp.psf_monomials,ublas::range(index,index+m_size)))=legendre_monomials_for_this_psf_parameter;
	index += m_size;
      }
    }
  }



}



double specex::PSF_Fitter::ComputeChi2AB(bool compute_ab, int input_begin_j, int input_end_j, harp::matrix_double* input_Ap, harp::vector_double* input_Bp, bool update_tmp_data) const  {

  
  
  
  
  int begin_j = input_begin_j;
  int end_j   = input_end_j;
  
  
  harp::matrix_double* Ap = input_Ap;
  harp::vector_double* Bp = input_Bp;
  
  if(begin_j==0) begin_j=stamp.begin_j;
  if(end_j==0) end_j=stamp.end_j;
  
  if(compute_ab) {
    if(Ap==0) Ap = & const_cast<specex::PSF_Fitter*>(this)->A_of_band[0];
    if(Bp==0) Bp = & const_cast<specex::PSF_Fitter*>(this)->B_of_band[0];
  }

  //===============================================
#define FASTER_THAN_SYR
#ifdef FASTER_THAN_SYR  
  
  vector<int> other_indices;
  harp::matrix_double Ablock;
  vector<harp::vector_double> Arect;
  bool do_faster_than_syr = compute_ab && fit_flux && spot_tmp_data.size()>1;
  if(do_faster_than_syr) {
    Ablock.resize(index_of_spots_parameters,index_of_spots_parameters);
    Ablock.clear();
    for(size_t s=0; s<spot_tmp_data.size(); s++) {
      //Arect.push_back(harp::vector_double(index_of_spots_parameters));
      //Arect.back().clear();    
      
      const SpotTmpData& tmp = spot_tmp_data[s];
      if(tmp.can_measure_flux && !tmp.ignore) {
	Arect.push_back(harp::vector_double(index_of_spots_parameters));
	Arect.back().clear();    
      }
    }
     
  }
#endif
  //=============================================== 
  
  if(update_tmp_data) const_cast<specex::PSF_Fitter*>(this)->UpdateTmpData(compute_ab);
  

  harp::vector_double H;  
  if(compute_ab) {
    H.resize(nparTot);
  }
  
  double chi2 = 0;
  int npix_in_chi2 = 0;
  double sum_flux = 0;
  
  harp::vector_double gradAllPar,gradPos;
  harp::vector_double *gradAllPar_pointer = 0;
  harp::vector_double *gradPos_pointer = 0;
  std::vector<int> indices_of_fitpar_in_allpar;
  
  if(compute_ab) {
    Ap->clear();
    Bp->clear();
    
    if(fit_psf || fit_psf_tail) {
      gradAllPar.resize(psf->LocalNAllPar()); 
      gradAllPar_pointer = &gradAllPar;

      indices_of_fitpar_in_allpar.resize(npar_fixed_coord); 

      const std::vector<Pol_p>& AP=psf_params->AllParPolXW;
      const std::vector<Pol_p>& FP=psf_params->FitParPolXW;

      // assert(npar_fixed_coord == FP.size()); // ok
      
      size_t fk=0;
      for (size_t ak =0; ak < AP.size(); ++ak) {
	const Pol_p FPk = FP[fk];
	const Pol_p APk = AP[ak];
	if(APk==FPk) {
	  indices_of_fitpar_in_allpar[fk]=int(ak);
	  fk++; // change free param index for next iteration
	  if(fk>=FP.size()) break;
	}
      }
    
    }// will remain zero if (!fit_psf)
    if(fit_position || fit_trace) {gradPos.resize(2); gradPos_pointer = &gradPos;}
  
    
    
  }
  
  bool use_footprint = (spot_tmp_data.size()>1 && footprint_weight.Nx()>0);
  

#ifdef CONTINUUM
  if(psf_params->fiber_min<psf_params->fiber_min) SPECEX_ERROR("fibers not defined");
  
  harp::vector_double continuum_params;
  harp::vector_double x_of_trace_for_continuum;
  //harp::vector_double w_of_trace_for_continuum;
  double expfact_for_continuum = 0;
  map<int,harp::vector_double> continuum_monomials;
  size_t np_continuum = psf_params->ContinuumPol.coeff.size();
  bool has_continuum  = fit_continuum;
  if(!has_continuum) for(size_t k=0; k<np_continuum; k++) if(psf_params->ContinuumPol.coeff[k]!=0) {has_continuum = true; break;}

  if(has_continuum) {
    if(fit_continuum)
      continuum_params = ublas::project(Params,ublas::range(continuum_index,continuum_index+np_continuum));
    else
      continuum_params = psf_params->ContinuumPol.coeff;
    
    expfact_for_continuum=1./(sqrt(2*M_PI)*psf_params->continuum_sigma_x);
    
    x_of_trace_for_continuum.resize(psf_params->fiber_max-psf_params->fiber_min+1);
  }
#endif

  for (int j=begin_j; j <end_j; ++j) {
    
#ifdef CONTINUUM
    if(has_continuum) {
      for(int fiber=psf_params->fiber_min;fiber<=psf_params->fiber_max;fiber++) {
	if(psf->GetTrace(fiber).Off()) continue;
	x_of_trace_for_continuum(fiber-psf_params->fiber_min) = psf->GetTrace(fiber).X_vs_Y.Value(j);
	continuum_monomials[fiber]=psf_params->ContinuumPol.Monomials(psf->GetTrace(fiber).W_vs_Y.Value(j));
      }
    }
#endif  

    int i1_side_band=0;
    int i2_side_band=0;
    if(increase_weight_of_side_bands && recompute_weight_in_fit) {
      i1_side_band=int(floor(psf->GetTrace(psf_params->fiber_min).X_vs_Y.Value(double(j))+0.5));
      i2_side_band=int(floor(psf->GetTrace(psf_params->fiber_max).X_vs_Y.Value(double(j))+0.5));
    }
    
    for (int i=stamp.begin_i ; i < stamp.end_i; ++i) {

      double w = 1;
      if(use_footprint) { 
	w=footprint_weight(i,j); 
      } else{
	w=weight(i,j);
      }
      if (w<=0) continue;

      bool in_side_band = false;
      if(increase_weight_of_side_bands && recompute_weight_in_fit) {
	if(i<=i1_side_band || i>=i2_side_band) in_side_band = true;
	
      }
      


      // w=1; // DEBUG

      double res = double(image(i,j));
      double signal = 0;

      if(compute_ab) {
	H.clear();
#ifdef FASTER_THAN_SYR
	other_indices.clear();
#endif
      }

#ifdef CONTINUUM
      if(has_continuum) {
	double continuum_value=0;
	for(int fiber=psf_params->fiber_min;fiber<=psf_params->fiber_max;fiber++) {
	  if(psf->GetTrace(fiber).Off()) continue;
	  double continuum_prof = expfact_for_continuum * exp(-0.5*square((i-x_of_trace_for_continuum(fiber-psf_params->fiber_min))/psf_params->continuum_sigma_x));
	  continuum_value += specex::dot(continuum_params,continuum_monomials[fiber])*continuum_prof;
	  if(compute_ab && fit_continuum) {
	    ublas::noalias(ublas::project(H,ublas::range(continuum_index,continuum_index+np_continuum))) += continuum_prof*continuum_monomials[fiber];
	  }
	}
	signal += continuum_value;
      }
#endif
      bool compute_tail = false;
#ifdef EXTERNAL_TAIL
      compute_tail = (fit_psf_tail || (psf_params->AllParPolXW[psf->ParamIndex("TAILAMP")]->coeff(0)!=0));
#endif      
      
      int nspots_in_pix = 0;

      

      for(size_t s=0;s<spot_tmp_data.size();s++) { // loop on tmp spots data
	
	const specex::SpotTmpData &tmp = spot_tmp_data[s];
	if (tmp.ignore) continue;
	bool in_core = tmp.stamp.Contains(i,j);

	//if( (!in_core) && (!fit_psf_tail)  ) continue; // if we fit tails we use the data outside the core, completely wrong, we need tails values everywhere
		
	nspots_in_pix++;
	


	double psfVal =  psf->PSFValueWithParamsXY(tmp.x,tmp.y, i, j, tmp.psf_all_params, gradPos_pointer, gradAllPar_pointer, in_core, compute_tail); // compute core part of psf only in core
	
	
	double flux = tmp.flux;
	
	if(!in_core) flux = tmp.frozen_flux; // to decorrelate tails

	
	if(fabs(flux*psfVal)>1.e20) {
	  SPECEX_ERROR("SEVERE BUG flux,psf,params " << flux << " " << psfVal << " " << tmp.psf_all_params);
	}

	
	signal += flux*psfVal;
	
	if (compute_ab) {
	  
	  if((fit_psf && in_core) || fit_psf_tail) {
	    size_t index = 0;
	    for(int p=0;p<npar_fixed_coord;p++) {
	      size_t m_size = psf_params->FitParPolXW[p]->coeff.size();
	      //blas::axpy(flux*gradPar[p],ublas::project(tmp.psf_monomials,ublas::range(index,index+m_size)),ublas::project(H,ublas::range(index,index+m_size))); // doesnt compile
	      ublas::noalias(ublas::project(H,ublas::range(index,index+m_size))) += (flux*gradAllPar(indices_of_fitpar_in_allpar[p]))*ublas::project(tmp.psf_monomials,ublas::range(index,index+m_size));
	      index += m_size;
	    }
	  }
	  if(fit_trace) {
	    ublas::noalias(ublas::project(H,ublas::range(tmp.trace_x_parameter_index,tmp.trace_x_parameter_index+tmp.trace_x_monomials.size()))) 
	      += (gradPos(0) * flux)*tmp.trace_x_monomials;
	    ublas::noalias(ublas::project(H,ublas::range(tmp.trace_y_parameter_index,tmp.trace_y_parameter_index+tmp.trace_y_monomials.size()))) 
	      += (gradPos(1) * flux)*tmp.trace_y_monomials;
	  }
	  
	  if(fit_flux && in_core) {

	    if(force_positive_flux)
	      H(tmp.flux_parameter_index) += tmp.flux*psfVal;
	    else
	      H(tmp.flux_parameter_index) += psfVal;

#ifdef FASTER_THAN_SYR
	    if(tmp.can_measure_flux)
	      other_indices.push_back(tmp.flux_parameter_index);
#endif
	  }
	  if(fit_position) {
	    H(tmp.x_parameter_index) += gradPos(0) * flux;
	    H(tmp.y_parameter_index) += gradPos(1) * flux;
#ifdef FASTER_THAN_SYR	    
	    other_indices.push_back(tmp.x_parameter_index);
	    other_indices.push_back(tmp.y_parameter_index);
#endif
	  }







	} // end of test compute_ab
      } // end of loop on spots
      
      
      double wscale=1;
      if(recompute_weight_in_fit) {
	w =  square(readnoise(i,j));
	if(signal>0) {
	   w += signal/psf->gain;
	   w += square(psf->psf_error*signal);	   
	}
	w = 1./w;

	// increase weight of side bands to avoid fiber to fiber degeneracy
	if(increase_weight_of_side_bands && in_side_band) {
	  w *= SIDE_BAND_WEIGHT_SCALE;
	  wscale *= SIDE_BAND_WEIGHT_SCALE;
	}
      }
      
      
      if(corefootprint_weight_boost>0 && corefootprint.n_rows()>0 && corefootprint(i,j)>0) {
	wscale *= corefootprint_weight_boost;
	w *= corefootprint_weight_boost;
      }
      
      res -= signal;
      
      chi2 += w*res*res;
      if(w>0) {
	npix_in_chi2 ++;
	sum_flux += signal;
      }
      
      /////////////////////////////////////////////////////
      // this is for debugging when developping the code
      //#define TESTING_H_VECTOR
#ifdef TESTING_H_VECTOR

      if(compute_ab && fit_psf && specex::dot(H,H)>0) {
	cout << "DEBUGGING TESTING_H_VECTOR" << endl;
	
	footprint_weight.data.clear();
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
	    double chi2_der_analytic = -2*w*resH(p); 
	  
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
	
	double bfact = w*res;
	if(recompute_weight_in_fit) {
	  if(signal>0)
	    bfact += (1./wscale)*0.5*square(w*res)*(1/psf->gain+2*square(psf->psf_error)*signal);
	}

	// doing A += w*h*h.transposed();  B += fact*h;

#ifdef FASTER_THAN_SYR
	if(do_faster_than_syr) {
	  if(index_of_spots_parameters==0) {
	    for(vector<int>::const_iterator i=other_indices.begin();i!=other_indices.end();i++) {
	      const double& hi = H(*i);
	      double whi = w*hi;
	      (*Ap)(*i,*i) += whi*hi;
	      for(vector<int>::const_iterator j=other_indices.begin();j!=i;j++) {
		(*Ap)(*i,*j) += whi*H(*j);
	      }
	    }
	  } else {
	    blas::syr(w,ublas::project(H,ublas::range(0,index_of_spots_parameters)),boost::numeric::bindings::lower(Ablock));
	    
	    for(vector<int>::const_iterator i=other_indices.begin();i!=other_indices.end();i++) {
	      const double& hi = H(*i);
	      double whi = w*hi;
	      (*Ap)(*i,*i) += whi*hi;
	      for(vector<int>::const_iterator j=other_indices.begin();j!=i;j++) (*Ap)(*i,*j) += whi*H(*j);
	      specex::axpy(whi,ublas::project(H,ublas::range(0,index_of_spots_parameters)),Arect[*i-index_of_spots_parameters]);
	    }
	  }
	}else{
	  specex::syr(w,H,*Ap);
	}
	specex::axpy(bfact,H,*Bp);
		
#else
	specex::syr(w,H,*Ap); //  this is still way too slow; this routine takes advantage of the zeros in h, its faster than sparse matrices
	specex::axpy(bfact,H,*Bp); 	
#endif
 
      } // end of test on compute_ab
      
    } // end of loop on pix coord. i
  } // end of loop on pix coord. j
  
#ifdef FASTER_THAN_SYR  
  if(do_faster_than_syr) {
    size_t n=Ablock.size1();
    for(size_t j=0;j<n;j++)
      for(size_t i=j;i<n;i++)
	(*Ap)(i,j) += Ablock(i,j);
    for(size_t s=0; s<Arect.size(); s++) {
      int i = index_of_spots_parameters + s;
      for(size_t j=0;j<n;j++)
	(*Ap)(i,j) += Arect[s](j);
    }
  }
#endif
  
  /*
  if(force_positive_flux) {
    
    double min_flux = 1.; // electrons
    double w = 100; //
    
    for(size_t s=0;s<spot_tmp_data.size();s++) { // loop on tmp spots data
      
      const specex::SpotTmpData &tmp = spot_tmp_data[s];
      
      if(tmp.flux>min_flux) continue;
      
      double res = min_flux-tmp.flux;
      chi2 += w*res*res;
      if(compute_ab && fit_flux) {
	(*Ap)(tmp.flux_parameter_index,tmp.flux_parameter_index) += w;
	(*Bp)(tmp.flux_parameter_index) += w*res;
      }
    }
  }
  */

  // psf priors , only once
  if(fit_psf && !(psf->Priors.empty()) && begin_j == stamp.begin_j) {
    
    int npar = psf->LocalNAllPar();

    for(size_t s=0;s<spot_tmp_data.size();s++) { // loop on tmp spots data
      
      const specex::SpotTmpData &tmp = spot_tmp_data[s];
      if(tmp.ignore) continue;
      
      int index=0;
      int fp=0; // fitted par. index
      for(int ap=0;ap<npar;ap++) {
	
	const specex::Pol_p AP=psf_params->AllParPolXW[ap];
	const specex::Pol_p FP=psf_params->FitParPolXW[fp];
	
	if(AP==FP) { // only apply prior to fit parameters  
	  
	  size_t c_size = AP->coeff.size();
	  std::map<int,Prior*>::const_iterator it = psf->Priors.find(ap);
       
	  if(it==psf->Priors.end()) {index += c_size; continue;} // no prior for this psf parameter
	

	  SPECEX_INFO("accounting for Prior at fit index " << fp << " and all index " << ap);
	  
	  const Prior* prior = it->second;
	  const double& par = tmp.psf_all_params(ap);
	
	  if (compute_ab) {
	    
	    
	    for (size_t c=0; c<c_size; c++, index++) {
	    
	      const double& monomial_val = tmp.psf_monomials(index);
	      (*Bp)(index)       += monomial_val * prior->hdChi2dx(par);
	      (*Ap)(index,index) += square(monomial_val) * prior->hd2Chi2dx2(par);
	    
	    }
	  }
	  chi2 += prior->Chi2(par);
	  fp++;
	  if(fp>=int(psf_params->FitParPolXW.size())) break;
	}
      }
    }

  }
  SPECEX_DEBUG("ComputeChi2AB chi2=" << chi2 << " npix=" << npix_in_chi2 << " sumflux=" << sum_flux);
  //SPECEX_INFO("ComputeChi2AB j range= " << input_begin_j << " " << input_end_j << " chi2= " << chi2);
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

  // definition of fitted region of image
  // ----------------------------------------------------
  stamp = compute_stamp(image,psf,spots,0,0,psf_params->bundle_id);
  
  if(spots.size()>1) {
    // compute psf footprint
    footprint_weight.resize(weight.Nx(),weight.Ny());
    footprint_weight.data.clear();
    
    *npix = 0;
    if(include_signal_in_weight) {
      
      SPECEX_INFO("WEIGHTS: Computing weights");
      SPECEX_INFO("WEIGHTS: gain = " << psf->gain);
      SPECEX_INFO("WEIGHTS: psf error = " << psf->psf_error);
      
      bool only_on_spots = !(fit_psf_tail || fit_continuum);
      bool only_psf_core = false;
      bool only_positive = true;
      
      if(only_on_spots) SPECEX_DEBUG("WEIGHTS: only on spots");
      if(only_psf_core) SPECEX_DEBUG("WEIGHTS: only psf core");
      if(only_positive) SPECEX_DEBUG("WEIGHTS: only positive");
      
      bool modified_tail_amplitude = false;
      harp::vector_double saved_tail_amplitude_coeff;
      if(fit_psf_tail || fit_continuum) {
	int index=psf->ParamIndex("TAILAMP");
	if(psf_params->AllParPolXW[index]->coeff(0)==0) {
	  saved_tail_amplitude_coeff=psf_params->AllParPolXW[index]->coeff;
	  psf_params->AllParPolXW[index]->coeff(0)=0.005;
	  modified_tail_amplitude = true;
	  SPECEX_DEBUG("WEIGHTS: setting tail amplitude = " << psf_params->AllParPolXW[index]->coeff(0) << " for weights");
	}
      }
      
      // generate error for a reason not understood
      parallelized_compute_model_image(footprint_weight,weight,psf,spots,only_on_spots,only_psf_core,only_positive,0,0,psf_params->bundle_id);
      
      //SPECEX_INFO("FOR DEBUG write model-image.fits");
      //write_new_fits_image("model-image.fits",footprint_weight);
      

      //compute_model_image(footprint_weight,weight,psf,spots,only_on_spots,only_psf_core,only_positive,-1,-1,0,0,psf_params->bundle_id);
      
      if(modified_tail_amplitude) {
	psf_params->AllParPolXW[psf->ParamIndex("TAILAMP")]->coeff = saved_tail_amplitude_coeff;
      }
      

      
      
      //if(fit_psf_tail || fit_continuum) SPECEX_INFO("debug, writing toto.fits and exit"); write_new_fits_image("toto.fits",footprint_weight); exit(12);
      

      // compute variance and weight
      
      SPECEX_DEBUG("WEIGHTS: Compute weights");
      
      for (int j=stamp.begin_j; j <stamp.end_j; ++j) {  // this stamp is a rectangle with all the spots in it
	for (int i=stamp.begin_i ; i < stamp.end_i; ++i) {

	  if( weight(i,j)==0) continue ; // pixels with ivar=0 still at 0
	  
	  double model_flux=footprint_weight(i,j);
	  if(model_flux==0) continue;
	  
	  double var = square(readnoise(i,j));
	  if(model_flux>0) { // else negative fluctuation
	    var += model_flux/psf->gain;
	    var += square(psf->psf_error*model_flux);
	  }
	  
	  footprint_weight(i,j) = 1./var;
	  (*npix)++;
	}
      }
      
      //SPECEX_INFO("FOR DEBUG write model-weight.fits and exit");
      //write_new_fits_image("model-weight.fits",footprint_weight);
      //exit(12);
      
    }else{
      SPECEX_DEBUG("WEIGHTS: inverse variance");
      if(fit_psf_tail || fit_continuum) {
	stamp = compute_stamp(image,psf,spots,0,0,psf_params->bundle_id);
	for(int j=stamp.begin_j;j<stamp.end_j;j++) {
	  int margin = min(MAX_X_MARGIN,psf->hSizeX); // 7 is half distance between center of ext. fibers of adjacent bundles
	  int begin_i = max(stamp.begin_i, int(floor(psf->GetTrace(psf_params->fiber_min).X_vs_Y.Value(double(j))+0.5))-margin);
	  int end_i   = min(stamp.end_i  , int(floor(psf->GetTrace(psf_params->fiber_max).X_vs_Y.Value(double(j))+0.5))+margin+1);
	  for(int i=begin_i;i<end_i;i++) {
	    footprint_weight(i,j)=weight(i,j);
	  }
	}
      }else{
	for(size_t s=0;s<spots.size();s++) {
	  specex::Spot_p& spot= spots[s];
	  Stamp spot_stamp(image);
	  SetStampLimitsFromPSF(spot_stamp,psf,spot->xc,spot->yc);
	  for(int j=spot_stamp.begin_j;j<spot_stamp.end_j;j++) {
	    int margin = min(MAX_X_MARGIN,psf->hSizeX); // 7 is half distance between center of ext. fibers of adjacent bundles
	    int begin_i = max(spot_stamp.begin_i, int(floor(psf->GetTrace(psf_params->fiber_min).X_vs_Y.Value(double(j))+0.5))-margin);
	    int end_i   = min(spot_stamp.end_i  , int(floor(psf->GetTrace(psf_params->fiber_max).X_vs_Y.Value(double(j))+0.5))+margin+1);
	    for(int i=begin_i;i<end_i;i++) {
	      footprint_weight(i,j)=weight(i,j);
	    }
	  }
	}
      }
    }
    bool zero_weight_for_core = ((fit_psf_tail || fit_continuum) && !fit_flux && !fit_psf);
      
      
    if(zero_weight_for_core) {
      
      SPECEX_INFO("WEIGHTS: Set weight to ZERO at the CORE of spots (to fit tails or continuum)");
      double saved_hsx = psf->hSizeX;
      double saved_hsy = psf->hSizeY;
      psf->hSizeX = min(12,psf->hSizeX);
      psf->hSizeY = min(12,psf->hSizeY);
      
      for(size_t s=0;s<spots.size();s++) {
	Stamp spot_stamp(image);
	psf->StampLimits(spots[s]->xc,spots[s]->yc,spot_stamp.begin_i,spot_stamp.end_i,spot_stamp.begin_j,spot_stamp.end_j);
	spot_stamp.begin_i = max(0,spot_stamp.begin_i);
	spot_stamp.end_i   = min(spot_stamp.Parent_n_cols(),spot_stamp.end_i);
	  spot_stamp.begin_j = max(0,spot_stamp.begin_j);
	  spot_stamp.end_j   = min(spot_stamp.Parent_n_rows(),spot_stamp.end_j);
	  //SPECEX_INFO("DEBUG " << spots[s]->xc << " " << spots[s]->yc << " " << spot_stamp.begin_i << " " << spot_stamp.end_i << " " << spot_stamp.begin_j << " " << spot_stamp.end_j);
	  for(int j=spot_stamp.begin_j;j<spot_stamp.end_j;j++) {
	    for(int i=spot_stamp.begin_i;i<spot_stamp.end_i;i++) {
	      footprint_weight(i,j) = 0;
	    }
	  }
      }
      psf->hSizeX = saved_hsx;
      psf->hSizeY = saved_hsy;
    }
    
    // mask out some regions with mis-understood lines
    mask.ApplyMaskToImage(footprint_weight,*psf,0);
    
    (*npix)=0;
    for (int j=stamp.begin_j; j <stamp.end_j; ++j) {  // this stamp is a rectangle with all the spots in it
      for (int i=stamp.begin_i ; i < stamp.end_i; ++i) {
	if(footprint_weight(i,j)>0) (*npix)++;
      }
    }
    
    SPECEX_DEBUG("WEIGHTS: number of pixels in fit = " << *npix);
    SPECEX_DEBUG("WEIGHTS: number of spots in fit  = " << spots.size());
    
    
    if(increase_weight_of_side_bands) {
      SPECEX_INFO("WEIGHTS: increase weight of side bands to avoid fiber to fiber degeneracy");
      
      int margin = min(MAX_X_MARGIN,psf->hSizeX); // 7 is half distance between center of ext. fibers of adjacent bundles
      
      
      int npix_side_band = 0;
      for(int j=stamp.begin_j;j<stamp.end_j;j++) {
	int i1 = int(floor(psf->GetTrace(psf_params->fiber_min).X_vs_Y.Value(double(j))+0.5));
	int begin1_i = max(stamp.begin_i,i1-margin);
	int end1_i   = i1+1;
	
	int i2 = int(floor(psf->GetTrace(psf_params->fiber_max).X_vs_Y.Value(double(j))+0.5));
	int begin2_i = i2;
	int end2_i   = min(stamp.end_i,i2+margin+1);
	
	for(int i=begin1_i;i<end1_i;i++) {
	  footprint_weight(i,j)*=SIDE_BAND_WEIGHT_SCALE;
	  npix_side_band++;
	}
	for(int i=begin2_i;i<end2_i;i++) {
	  footprint_weight(i,j)*=SIDE_BAND_WEIGHT_SCALE;
	  npix_side_band++;
	}
      }
      *npix += (SIDE_BAND_WEIGHT_SCALE-1)*npix_side_band; // artificial scaling to get meaningfull chi2
    } // end of increase of weight
    
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
    
  int npar_psf = 0;
  if(fit_psf) npar_psf = psf->BundleNFitPar(psf_params->bundle_id);
  npar_trace = 0;
  
  
  npar_fixed_coord = psf_params->FitParPolXW.size();
  npar_varying_coord = psf->BundleNFitPar(psf_params->bundle_id);

  nparTot  = NPar(spots.size());

  ComputeWeigthImage(spots,npix);
  
  InitTmpData(spots);
  int n_to_attach = 0;
  for(size_t s=0;s<spot_tmp_data.size();s++) {
    SpotTmpData& tmp = spot_tmp_data[s];
    if(!tmp.can_measure_flux) n_to_attach++;
  }
  if(fit_flux)
    nparTot -= n_to_attach;
  
  
  

  SPECEX_DEBUG("specex::PSF_Fitter::FitSeveralSpots npar_fixed_coord   = " << npar_fixed_coord);
  SPECEX_DEBUG("specex::PSF_Fitter::FitSeveralSpots npar_varying_coord = " << npar_varying_coord);
  
  
  
  // allocation of parameters, setting Params 
  // and recording some indices
  // ----------------------------------------------------
  Params.resize(nparTot);
  Params.clear();

  map<int,int> tmp_trace_x_parameter;
  map<int,int> tmp_trace_y_parameter;
  index_of_spots_parameters = 0;
  
  {
    int index=0;
    if(fit_psf || fit_psf_tail) {
      for(size_t p=0;p<psf_params->FitParPolXW.size();p++) {
	const harp::vector_double& coeff=psf_params->FitParPolXW[p]->coeff;
	size_t c_size = coeff.size();
	ublas::noalias(ublas::project(Params,ublas::range(index,index+c_size))) = coeff;
	index += c_size;
      }
    }
    
    if(fit_trace) {  
      for(std::map<int,specex::Trace>::const_iterator it=psf->FiberTraces.begin(); it!=psf->FiberTraces.end(); ++it) {
		
	if(it->first < psf_params->fiber_min || it->first > psf_params->fiber_max) continue;
	
	if(it->second.Off()) continue; // not fitted

	tmp_trace_x_parameter[it->first] = index;

	{
	  const harp::vector_double& coeff = it->second.X_vs_W.coeff;
	  size_t c_size = coeff.size();
	  ublas::noalias(ublas::project(Params,ublas::range(index,index+c_size))) = coeff;
	  index += c_size;
	  npar_trace += c_size;
	}
	
	tmp_trace_y_parameter[it->first] = index;
	
	{
	  const harp::vector_double& coeff = it->second.Y_vs_W.coeff;
	   size_t c_size = coeff.size();
	   ublas::noalias(ublas::project(Params,ublas::range(index,index+c_size))) = coeff;
	  index += c_size;
	  npar_trace += c_size;
	}
      }
    }

#ifdef CONTINUUM
    if(fit_continuum) {
      continuum_index = index;
      ublas::noalias(ublas::project(Params,ublas::range(continuum_index,continuum_index+psf_params->ContinuumPol.coeff.size()))) = psf_params->ContinuumPol.coeff;
      index += psf_params->ContinuumPol.coeff.size();
    }
    
#endif
    index_of_spots_parameters = index;
  }


  



  // ----------------------------------------------------
  
  //ComputeWeigthImage(spots,npix);
  
  
  
  


  //InitTmpData(spots);
  
  if(corefootprint_weight_boost>0) {
    corefootprint = image;
    corefootprint.data.clear();
    int core_hsize=2;
    for(size_t s=0;s<spot_tmp_data.size();s++) {
      SpotTmpData& tmp = spot_tmp_data[s];
      int iPix = int(floor(tmp.x+0.5));
      int jPix = int(floor(tmp.y+0.5));
      int BeginI = max(iPix-core_hsize,tmp.stamp.begin_i);
      int BeginJ = max(jPix-core_hsize,tmp.stamp.begin_j);
      int EndI   = min(iPix+core_hsize+1,tmp.stamp.end_i);
      int EndJ   = min(jPix+core_hsize+1,tmp.stamp.end_j); 
      for(int j= BeginJ; j<EndJ; j++) {
	for(int i= BeginI; i<EndI; i++) {
	  corefootprint(i,j)=1;
	}
      }
    }
  }
  
  // indexation and Params and monomials for tmp spots 
  {
  int index = index_of_spots_parameters;
  for(size_t s=0;s<spot_tmp_data.size();s++) {
    SpotTmpData& tmp = spot_tmp_data[s];

    if(fit_flux) {
      if(tmp.can_measure_flux) {
	tmp.flux_parameter_index = index;
	index++;  
      
      if(force_positive_flux) {
	if(tmp.flux<=1) tmp.flux=1;
	Params(tmp.flux_parameter_index) = log(tmp.flux);
      }else{
	Params(tmp.flux_parameter_index) = tmp.flux;
      }
      }
    }
    if(fit_position) {
      tmp.x_parameter_index = index; 
      Params(tmp.x_parameter_index) = tmp.x;
      index++;
      tmp.y_parameter_index = index; 
      Params(tmp.y_parameter_index) = tmp.y;
      index++;
    }
    if(fit_psf || fit_psf_tail) {
      tmp.psf_monomials.resize(npar_varying_coord);
      // specex::zero(tmp.psf_monomials);
      int index=0;
      for(int p=0;p<npar_fixed_coord;p++) {
	harp::vector_double legendre_monomials_for_this_psf_parameter = psf_params->FitParPolXW[p]->Monomials(tmp.x,tmp.wavelength);
	size_t m_size = legendre_monomials_for_this_psf_parameter.size();
	ublas::noalias(ublas::project(tmp.psf_monomials,ublas::range(index,index+m_size)))=legendre_monomials_for_this_psf_parameter;
	index += m_size;
	}
    }    
    // trace monomials
    if(fit_trace) {
      const specex::Trace& trace = psf->FiberTraces[tmp.fiber];
      tmp.trace_x_monomials = trace.X_vs_W.Monomials(tmp.wavelength);
      tmp.trace_y_monomials = trace.Y_vs_W.Monomials(tmp.wavelength);
      tmp.trace_x_parameter_index = tmp_trace_x_parameter[tmp.fiber];
      tmp.trace_y_parameter_index = tmp_trace_y_parameter[tmp.fiber];
    }
  }
  }
  
  // now deal with spots for which we cannot measure fluxes

  
  if(n_to_attach>0) {
    SPECEX_INFO("Need to attach " << n_to_attach << " spots");
    for(size_t s=0;s<spot_tmp_data.size(); s++) {
    
      SpotTmpData& tmp = spot_tmp_data[s];
      if(!tmp.can_measure_flux) {
	// find a nice neighbour
	SpotTmpData* neighbour = 0;
	int fiber_diff = 100000;
	for(size_t s2=0;s2<spot_tmp_data.size();s2++) {
	  SpotTmpData& tmp2 = spot_tmp_data[s2];
	  if(fabs(tmp2.wavelength - tmp.wavelength)>0.00000001) continue;
	  if(!tmp2.can_measure_flux) continue;
	  if(fabs(tmp2.x-tmp.x)>2*psf->hSizeX) continue;
	  if(fabs(tmp2.y-tmp.y)>2*psf->hSizeY) continue;
	  int tmp_fiber_diff = abs(tmp2.fiber-tmp.fiber);
	  if(fiber_diff > tmp_fiber_diff) {
	    neighbour  = &tmp2;
	    fiber_diff = tmp_fiber_diff;
	  }
	  if(fiber_diff==1) break; // it's ok
	}	
	if(!neighbour) continue;
	SPECEX_DEBUG("use neighbour flux");
	tmp.flux_parameter_index = neighbour->flux_parameter_index;
	tmp.flux = neighbour->flux;
	tmp.ignore = false;
      }
    }
  }

  
  SPECEX_INFO("specex::PSF_Fitter::FitSeveralSpots inc. signal in w=" << include_signal_in_weight << ", npix footprint = " << *npix);




  ////////////////////////////////////////////////////////////////////////// 
  number_of_image_bands = 1;
  if(parallelized) {
    char* OMP_NUM_THREADS = getenv("OMP_NUM_THREADS");
    if(OMP_NUM_THREADS) {
      number_of_image_bands = atoi(OMP_NUM_THREADS);
      SPECEX_DEBUG("Using " << number_of_image_bands << " image bands equal to value of OMP_NUM_THREADS");
    }
  }
 
  
  A_of_band.clear();
  B_of_band.clear();
  for(int i=0;i<number_of_image_bands;i++) {
    A_of_band.push_back(harp::matrix_double(nparTot, nparTot));
    B_of_band.push_back(harp::vector_double(nparTot));
  }


  //////////////////////////////////////////////////////////////////////////
  
  
  while(true) { // minimization loop 
      
    oldChi2 = *psfChi2;
    SPECEX_DEBUG("specex::PSF_Fitter::FitSeveralSpots iter=" << *niter << " old chi2=" << oldChi2);
    
    clock_t tstart = clock();
    SPECEX_DEBUG("specex::PSF_Fitter::FitSeveralSpots filling matrix (n="<< nparTot << ")...");
    if(parallelized) 
      *psfChi2 = ParallelizedComputeChi2AB(true);
    else
      *psfChi2 = ComputeChi2AB(true);
    clock_t tstop = clock();
    SPECEX_DEBUG("specex::PSF_Fitter::FitSeveralSpots time = " << float(tstop-tstart)/float(CLOCKS_PER_SEC) << "s");
    
    oldChi2 = *psfChi2;
    SPECEX_DEBUG("specex::PSF_Fitter::FitSeveralSpots chi2 = " << *psfChi2);
    SPECEX_DEBUG("specex::PSF_Fitter::FitSeveralSpots solving ...");
    
#ifdef SPARSE_H
#warning copy sparse mat and vect to dense for the moment 
    harp::matrix_double A = A_of_band[0];
    harp::vector_double B = B_of_band[0];
#else
    harp::matrix_double& A = A_of_band[0];
    harp::vector_double& B = B_of_band[0];
#endif 

    harp::matrix_double As=A;
    //harp::vector_double Bs=B;
    
    /*
    if(0 && fit_psf && npar_fixed_coord>2) {
      {
	ofstream os("A.xml");
	boost::archive::xml_oarchive xml_oa ( os );
	xml_oa << BOOST_SERIALIZATION_NVP(A);
	os.close();
      }
      {
	ofstream os("B.xml");
	boost::archive::xml_oarchive xml_oa ( os );
	xml_oa << BOOST_SERIALIZATION_NVP(B);
	os.close();
      }
      exit(12);
    }
    */

    
    //if(fit_psf && fit_flux) { specex::write_new_fits_image("A.fits",A); exit(12); }
    
    int status = cholesky_solve(A,B);

    SPECEX_DEBUG("specex::PSF_Fitter::FitSeveralSpots solving done");

    if (status != 0) {
      *psfChi2 = 1e30;
      for(size_t i=0;i<As.size1();i++) {
	if(As(i,i)<=0) {
	  SPECEX_DEBUG("DEBUG A(" <<i << "," << i << ")=" << As(i,i));
	  if(fit_flux) { // find it
	    for(size_t s=0;s<spot_tmp_data.size();s++) {
	      SpotTmpData& tmp = spot_tmp_data[s];
	      if(i==tmp.flux_parameter_index) {
		SPECEX_DEBUG(i << " it's spot x,y= " << tmp.x << "," << tmp.y << " fiber=" << tmp.fiber << " flux=" << tmp.flux);
	      }
	    }
	  }
	}
      }
      /*
	if(As.size1()>1) {
	specex::write_new_fits_image("A.fits",As);
	SPECEX_WARNING("wrote A.fits");
      }
      */
      if(fatal) {
	SPECEX_ERROR("cholesky_solve failed with status " << status);
      } else {
	if(spot_tmp_data.size()>1)
	  SPECEX_WARNING("cholesky_solve failed with status " << status);
	psf_params->fit_status = 1;
	psf_params->chi2 = *psfChi2;
	psf_params->nparams = Params.size();
	psf_params->ndata = *npix;
	psf_params->nspots_in_fit = spots.size();
	return false;
      }
    }
    
    
    
    bool linear = false;
    
    if( (fit_flux) && (!fit_position) && (!fit_trace) && (!fit_psf) ) {
      SPECEX_DEBUG("specex::PSF_Fitter::FitSeveralSpots linear because only fit flux");
      linear = true;
    }
    if( (fit_psf) && (!fit_position) && (!fit_trace) && (!fit_flux) && psf->IsLinear()) {
      SPECEX_DEBUG("specex::PSF_Fitter::FitSeveralSpots linear because only fit psf (linear wrt params)");
      linear = true;
    } 
    if( (fit_psf_tail || fit_continuum) && (!fit_position) && (!fit_trace) && (!fit_flux) && (!fit_psf)) {
      SPECEX_DEBUG("specex::PSF_Fitter::FitSeveralSpots linear because only fit tail and continuum");
      linear = true;
    } 

    
    bool use_brent = true;
    
    if(1) {
      if(linear) {
	SPECEX_DEBUG("specex::PSF_Fitter::FitSeveralSpots no brent because linear");
	use_brent = false;
      }else{
	if(!fit_trace) { // we use brent anyway for the fit of traces
	// check whether step indeed decreases chi2
	Params += B;
	SPECEX_DEBUG("specex::PSF_Fitter::FitSeveralSpots computing chi2 for step=1 ...");
	double chi2_1 = ParallelizedComputeChi2AB(false);
	Params -= B;
	if(chi2_1>oldChi2) {
	  use_brent = true;
	  SPECEX_DEBUG("specex::PSF_Fitter::FitSeveralSpots use brent because step 1 increases chi2");
	}else{
	  
	  use_brent = false;
	  *psfChi2 = chi2_1;
	  SPECEX_DEBUG("specex::PSF_Fitter::FitSeveralSpots don't brent because chi2 decreases");	  
	}
	}
      }
    }
    if(use_brent) {
      double brent_precision = 0.01;
      
      SPECEX_DEBUG("specex::PSF_Fitter::FitSeveralSpots starting brent with precision = " << brent_precision << " ...");
     
      
      
      if(fit_flux && force_positive_flux) {
	// force step for flux to be smaller than :
	double max_step = log(100.);
	for(size_t s=0;s<spot_tmp_data.size();s++) {
	  SpotTmpData& tmp = spot_tmp_data[s];
	  double step = fabs(B(tmp.flux_parameter_index));
	  if(step>max_step) B *= (max_step/step);
	}
      }
      if(fit_trace) {
	// don't want a step larger than N pix for all spots (can take some time to compute)
	double max_step=0.2; //pixel
	//SPECEX_INFO("trace fit:check max step lower than " << max_step << " pix");
	
	double scale=1;
	for(size_t s=0;s<spot_tmp_data.size();s++) { // loop on tmp spots data
	  const specex::SpotTmpData &tmp = spot_tmp_data[s];
	  
	  double dx = specex::dot(ublas::project(scale*B,ublas::range(tmp.trace_x_parameter_index,tmp.trace_x_parameter_index+tmp.trace_x_monomials.size())),tmp.trace_x_monomials);
	  double dy = specex::dot(ublas::project(scale*B,ublas::range(tmp.trace_y_parameter_index,tmp.trace_y_parameter_index+tmp.trace_y_monomials.size())),tmp.trace_y_monomials);
	  double dist=sqrt(dx*dx+dy*dy);
	  //SPECEX_INFO("trace fit: dx=" << dx << " " << dy << " for x,y,fiber,wavelength " << tmp.x << "," << tmp.y << "," << tmp.fiber << "," << tmp.wavelength);
	  if(dist>max_step)
	    scale *= (max_step/dist);
	}
	if(scale<1) {
	  SPECEX_INFO("trace fit: scaling down parameter step by " << scale);
	  B*=scale;
	}
      }
      
      // need to use brent here
      BrentBox bbox(*this,B,spots);
      
      double min_step=-0.05;
      double max_step=1.001;
      double prefered_step=1;
      int status=0;
      
      // check the chi2 decrement is not good enought with step=1
      double best_chi2=compute_chi2_for_a_given_step(1,&bbox);      
      double best_step = 1;
      
      if(fabs(best_chi2-*psfChi2)>brent_precision) { // really try brent now
	best_step = brent((AnalyticFunction*)(compute_chi2_for_a_given_step),
			  min_step,prefered_step,max_step,
			  brent_precision,&bbox,best_chi2,status,100);
      }else{
	SPECEX_DEBUG("brent not needed for required chi2 decrement dchi2=" << (-best_chi2+*psfChi2));
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
	    psf_params->fit_status = 2;
	    psf_params->chi2 = *psfChi2;
	    psf_params->nparams = Params.size();
	    psf_params->ndata = *npix;
	    psf_params->nspots_in_fit = spots.size(); 
	    return false;
	  }
	} else {
	  SPECEX_WARNING("specex::PSF_Fitter::FitSeveralSpots brent improves things so we continue");
	}
      }
      if(specex_is_verbose()) {
	double dchi2 = (-best_chi2+*psfChi2);
	if(dchi2>0 || fabs(dchi2)<brent_precision) {
	  SPECEX_INFO("specex::PSF_Fitter::FitSeveralSpots successful brent fit with step = " << best_step << " dchi2=" << dchi2
		      << " chi2pdf = " << best_chi2/(*npix-Params.size()) << " npar = " << Params.size()
		      );
	}else{
	  if(fatal) {
	    SPECEX_ERROR("problem with brent dchi2 = " << dchi2);
	  } else {
	    SPECEX_WARNING("problem with brent dchi2 = " << dchi2);
	    return false;
	  }
	  best_step = 0;
	  best_chi2 = *psfChi2;
	}
      }
      
      Params += best_step*B;
      *psfChi2 = best_chi2;
	
    } else { // didn't use brent
      Params += B;
      // *psfChi2 = ComputeChi2AB(false); // already computed above
      SPECEX_INFO("specex::PSF_Fitter::FitSeveralSpots dchi2=" << oldChi2-*psfChi2 << " chi2pdf = " << *psfChi2/(*npix-Params.size()) << " npar = " << Params.size());
      
    }
    
    // some sanity checks flux == nan is really bad
    for (unsigned k=0; k < nparTot; ++k) {
      if (std::isnan(Params(k))) {
	if(fatal) {
	  SPECEX_ERROR("specex::PSF_Fitter::FitSeveralSpots one parameter read nan");
	} else {
	  SPECEX_WARNING("specex::PSF_Fitter::FitSeveralSpots one parameter read nan");
	}
	*psfChi2 = 1e30;
	psf_params->fit_status = 3;
	psf_params->chi2 = *psfChi2;
	psf_params->nparams = Params.size();
	psf_params->ndata = *npix;
	psf_params->nspots_in_fit = spots.size();
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
      SPECEX_DEBUG("specex::PSF_Fitter::FitSeveralSpots quit loop because dchi2 small");
      break;
    }
    if( linear ) {
      SPECEX_DEBUG("specex::PSF_Fitter::FitSeveralSpots quit loop because linear");
      break;
    }
  } // end of minimization loop
    
  // We have to extract the weight matrix of the psf parameters (the first npar of params vector).
  //  This involves "marginalization" over the position and flux of the star. Compute 
  // full covariance matrix (from chi2 second derivatives), extract the psf params sub block,
  // and invert it back to get a weight matrix
  

  
  fitWeight = A_of_band[0];
  if (specex::cholesky_invert_after_decomposition(fitWeight) != 0) {
    SPECEX_ERROR("cholesky_invert_after_decomposition failed");
  }
  
  // fitWeight is now a covariance matrix
  // just to avoid confusion :
  harp::matrix_double& fitCovmat = fitWeight;
  
  SPECEX_DEBUG("specex::PSF_Fitter::FitSeveralSpots saving fitted params");
  
  int index=0;

  // copy fitted parameters in the right places:
  if (fit_psf || fit_psf_tail) {
    // save params to psf
    
    //SPECEX_INFO("TESTING : input params of first spot = " << psf->AllLocalParamsFW(spots[0]->fiber,spots[0]->wavelength,spots[0]->fiber_bundle));
    
    UpdateTmpData(false);
    
    //SPECEX_INFO("TESTING : fit params of first spot = " << spot_tmp_data[0].psf_all_params);
		

    for(size_t p=0;p<psf_params->FitParPolXW.size();p++) {
      harp::vector_double& coeff=psf_params->FitParPolXW[p]->coeff;
      for(size_t c=0;c<coeff.size();c++,index++)
	coeff(c)=Params(index);
    }
  
    //SPECEX_INFO("TESTING : output params of first spot = " << psf->AllLocalParamsFW(spots[0]->fiber,spots[0]->wavelength,spots[0]->fiber_bundle));
    
  }
  if (fit_trace) {
    for(std::map<int,specex::Trace>::iterator it=psf->FiberTraces.begin(); it!=psf->FiberTraces.end(); ++it) {

      if(it->first < psf_params->fiber_min || it->first > psf_params->fiber_max) continue;
      if(it->second.Off()) continue;
      
      harp::vector_double& PX = it->second.X_vs_W.coeff;
      for(size_t k=0;k<PX.size();k++,index++)
	PX(k)=Params(index);
      harp::vector_double& PY = it->second.Y_vs_W.coeff;
      for(size_t k=0;k<PY.size();k++,index++)
	PY(k)=Params(index);
    }
  }

#ifdef CONTINUUM
    if(fit_continuum) {
      ublas::noalias(psf_params->ContinuumPol.coeff) = ublas::project(Params,ublas::range(continuum_index,continuum_index+psf_params->ContinuumPol.coeff.size()));
    }
#endif


  // save results for spots
    SPECEX_DEBUG("specex::PSF_Fitter::FitSeveralSpots saving spots fluxes nspots=" << spots.size() << " ntmp=" << spot_tmp_data.size());
  for(size_t s=0;s<spots.size();s++) {
    specex::Spot_p& spot= spots[s];
    specex::SpotTmpData& tmp = spot_tmp_data[s]; // BUG IF WE ERASED DATA ???
    if(fit_flux) {
      if(force_positive_flux) {
	spot->flux = exp(min(max(Params(tmp.flux_parameter_index),-30.),+30.));
	spot->eflux = spot->flux*sqrt(fitCovmat(tmp.flux_parameter_index,tmp.flux_parameter_index));
      }else{
	spot->flux = Params(tmp.flux_parameter_index); 
	spot->eflux = sqrt(fitCovmat(tmp.flux_parameter_index,tmp.flux_parameter_index));
      }
    }
    if(fit_position) {
      spot->xc = Params(tmp.x_parameter_index);
      spot->yc = Params(tmp.y_parameter_index);
    }
    if(fit_trace) {
      spot->xc = psf->Xccd(spot->fiber,spot->wavelength);
      spot->yc = psf->Yccd(spot->fiber,spot->wavelength);

      //if(spot->xc != tmp.x || spot->yc != tmp.y) {
      //SPECEX_ERROR("spot coord. " << spot->xc << "," << spot->yc << " tmp " << tmp.x << "," << tmp.y);
      //}

    }
  }
    
  // spot->chi2=*psfChi2;
  
  bool ok=(*niter < maxiter);
  if(specex_is_verbose()) {
    
    if(ok) {
      if(specex_is_verbose()) {
	cout << "INFO specex::PSF_Fitter::FitSeveralSpots successful fit of ";
	if(fit_flux) cout << "+flux";
	if(fit_position) cout << "+pos";
	if(fit_psf) cout << "+psf";
	if(fit_trace) cout << "+trace";
	
#ifdef EXTERNAL_TAIL
	if(fit_psf_tail) cout << "+tail";
#endif
#ifdef CONTINUUM
	if(fit_continuum) cout << "+continuum";
#endif

	cout << " chi2= " << *psfChi2;
	cout << " niter=" << *niter;
	
	cout << endl;
	
#ifdef EXTERNAL_TAIL
	if(fit_psf_tail) {
	  
	  harp::vector_double spot_params = psf->AllLocalParamsXW_with_FitBundleParams(spot_tmp_data[spot_tmp_data.size()/2].x
										       ,spot_tmp_data[spot_tmp_data.size()/2].wavelength,psf_params->bundle_id,Params);
	  SPECEX_INFO("psf tail amplitudes, " << spot_params(psf->ParamIndex("TAILAMP")));
	}
#endif
	
#ifdef CONTINUUM
	if(fit_continuum)
	  SPECEX_INFO("continuum amplitude, " << Params(continuum_index) );
#endif

      }
    }else{
      SPECEX_WARNING("specex::PSF_Fitter::FitSeveralSpots failed because reached max number of iterations");
    }

  }
  
  psf_params->fit_status = 0;
  psf_params->chi2 = *psfChi2;
  psf_params->nparams = Params.size();
  psf_params->ndata = *npix;
  psf_params->nspots_in_fit = spots.size();
  return ok;
};
 
bool specex::PSF_Fitter::FitOneSpot(specex::Spot_p& spot, double *chi2, int *n_iterations) {

  specex::Spot saved_spot = *spot;

  vector<specex::Spot_p> spots; 
  spots.push_back(spot);
  
  
  // int npar_psf = psf->LocalNPar();
  // psf_global_params->clear();
  // for(int p=0;p<npar_psf;p++) {
  //   psf_global_params->push_back(specex::Legendre2DPol(0,0,image.Nx(),0,0,image.Ny()));
  //   (*psf_global_params)[p].coeff(0) = spot->PSFParams(p);
  // }
  
  parallelized=false;
  bool ok = FitSeveralSpots(spots,&(spot->chi2),n_iterations);
  parallelized=true;
  
  if(false) {
    
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

/*
void specex::PSF_Fitter::SetAllPSFParams(const harp::vector_double &ParamsToSet) {
  // also, define psf polynomes (should be in fitter)
  psf_params->AllParPolXW.clear();
  psf_params->FitParPolXW.clear();
  for (size_t k =0; k < ParamsToSet.size(); ++k) {
    psf_params->AllParPolXW.push_back(new specex::Legendre2DPol(0,0,image.Nx(),0,0,image.Ny()));
    psf_global_params->push_back(specex::Legendre2DPol(0,0,image.Nx(),0,0,image.Ny()));
    psf_global_params->back().coeff(0)=ParamsToSet(k);
  }
}
*/


bool specex::PSF_Fitter::FitTraces(vector<specex::Spot_p>& spots, int *n_fibers_fitted) {
  
  SPECEX_DEBUG("specex::PSF_Fitter::FitTraces starting");

  int nok_memory_slot;
  int *nok = &nok_memory_slot;
  if(n_fibers_fitted) nok = n_fibers_fitted;
  
  *nok = 0;
  for(map<int,specex::Trace>::iterator it=psf->FiberTraces.begin();
      it !=psf->FiberTraces.end(); ++it) {
    
    if(it->first < psf_params->fiber_min || it->first > psf_params->fiber_max) continue;
    if(it->second.Off()) continue;

    specex::Trace& trace=it->second;
    
    vector<specex::Spot_p> spots_of_fiber;
    for(size_t s=0;s<spots.size();s++) {
      specex::Spot_p& spot = spots[s];
      if(spot->fiber==trace.fiber)
	spots_of_fiber.push_back(spot);
    }
    bool ok = true;
    if(! trace.synchronized)
      ok = trace.Fit(spots_of_fiber);
    else
      SPECEX_DEBUG("No need to refit synchronized trace (wavemin,wavemax = " << trace.X_vs_W.xmin << "," << trace.X_vs_W.xmax << ")");
    
    if(ok) { 
      (*nok)++;
      
      for(size_t s=0;s<spots_of_fiber.size();s++) {
	specex::Spot_p& spot = spots_of_fiber[s];
	spot->xc=trace.X_vs_W.Value(spot->wavelength);
	spot->yc=trace.Y_vs_W.Value(spot->wavelength);
	//SPECEX_INFO("TRACE " << trace.fiber << " wave=" << spot->wavelength << " x=" << spot->xc << " y=" << spot->yc);
      }


    }
  }
  SPECEX_DEBUG("specex::PSF_Fitter::FitTraces ended nok = " << (*nok) << "/" << psf->FiberTraces.size());
 
  // detect failed traces via interpolation
  if(false && *nok>2) { // this may be problematic (need more check)

    map<int,specex::Trace>::iterator it_begin = psf->FiberTraces.begin();
    map<int,specex::Trace>::iterator it_end   = psf->FiberTraces.end();
    map<int,specex::Trace>::iterator it_last  = it_end; it_last--;
    
    for(map<int,specex::Trace>::iterator it=it_begin;it!=it_end;it++) {
      
      if(it->first < psf_params->fiber_min || it->first > psf_params->fiber_max) continue;
      if(it->second.Off()) continue;
      
      
      int fiber = it->first;
      int fiber1=-1;
      int fiber2=-1;      
      {
	vector<int> others;
	for(map<int,specex::Trace>::iterator it2=it_begin;it2!=it_end;it2++) {	
	  if(it2->first < psf_params->fiber_min || it2->first > psf_params->fiber_max) continue;
	  if(it2->second.Off()) continue;
	  if(it2==it) continue;
	  others.push_back(it2->first);
	}
	int mdist=1000;
	for(size_t f=0;f<others.size();f++) {
	  int dist=abs(others[f]-fiber);
	  if(dist<mdist) {
	    fiber1=others[f];
	    mdist=dist;
	  }
	}
	if(fiber1==-1) SPECEX_ERROR("could not find another valid fiber case #a, should not happen, for fiber " << fiber);
	mdist=1000;
	for(size_t f=0;f<others.size();f++) {
	  if(others[f]==fiber1) continue;
	  int dist=abs(others[f]-fiber);
	  if(dist<mdist) {
	    fiber2=others[f];
	    mdist=dist;
	  }
	}
	if(fiber2==-1) SPECEX_ERROR("could not find another valid fiber case #a, should not happen, for fiber " << fiber);
      }
            
      specex::Trace& trace=it->second;
      specex::Trace& trace1=psf->FiberTraces.find(fiber1)->second;
      specex::Trace& trace2=psf->FiberTraces.find(fiber2)->second;
      double f =  double(fiber);
      double f1 = double(fiber1);
      double f2 = double(fiber2);
      
      double wmin = min(min(trace.X_vs_W.xmin,trace1.X_vs_W.xmin),trace2.X_vs_W.xmin);
      double wmax = max(max(trace.X_vs_W.xmax,trace1.X_vs_W.xmax),trace2.X_vs_W.xmax);
      
      double dx=0;
      double dy=0;
      for(double w=wmin; w<wmax; w+=10) {
	double dxw=trace.X_vs_W.Value(w)-( (f2-f)*trace1.X_vs_W.Value(w) + (f-f1)*trace2.X_vs_W.Value(w) )/(f2-f1);
	double dyw=trace.Y_vs_W.Value(w)-( (f2-f)*trace1.Y_vs_W.Value(w) + (f-f1)*trace2.Y_vs_W.Value(w) )/(f2-f1);
	
	if(fabs(dxw)>fabs(dx)) dx=dxw;
	if(fabs(dyw)>fabs(dy)) dy=dxw;
	
      }
      SPECEX_INFO("checking fiber " << it->first << " using " << fiber1 << " and " << fiber2 << " dx=" << dx << " dy=" << dy);
      
      if(fabs(dx)>2 || fabs(dy)>2) {
	SPECEX_INFO("replacing trace of fiber " << it->first << " by interpolation");
	SPECEX_WARNING("replacing trace of fiber " << it->first << " by interpolation");
	trace.X_vs_W.deg  = min(trace1.X_vs_W.deg,trace2.X_vs_W.deg); trace.X_vs_W.coeff.resize(trace.X_vs_W.deg+1);
	trace.X_vs_W.coeff = ((f2-f)/(f2-f1))*trace1.X_vs_W.coeff + ((f-f1)/(f2-f1))*trace2.X_vs_W.coeff;
	trace.Y_vs_W.deg  = min(trace1.Y_vs_W.deg,trace2.Y_vs_W.deg); trace.Y_vs_W.coeff.resize(trace.Y_vs_W.deg+1);
	trace.Y_vs_W.coeff = ((f2-f)/(f2-f1))*trace1.Y_vs_W.coeff + ((f-f1)/(f2-f1))*trace2.Y_vs_W.coeff;
	trace.W_vs_Y = trace.X_vs_W.Invert(2); 
	trace.X_vs_Y = composed_pol(trace.X_vs_W,trace.W_vs_Y);
      }


    } // end of loop on fibers
  } // end of test

  
  return true;
}

bool specex::PSF_Fitter::FitIndividualSpotFluxes(std::vector<specex::Spot_p>& spots) {
  
  SPECEX_INFO("fitting independently the flux of each spot");

  // TURN OFF ALL MESSAGES HERE
  bool saved_debug = specex_is_debug();
  bool saved_info  = specex_is_verbose();
  specex_set_debug(false);
  specex_set_verbose(false);
  
  
  fit_flux                 = true;
  fit_position             = false;
  fit_psf                  = false;
  fit_trace                = false;
  fatal                    = false;

  bool saved_force_positive_flux = force_positive_flux;
  force_positive_flux = false;

#ifdef CONTINUUM
  fit_continuum  = false;
#endif
#ifdef EXTERNAL_TAIL
  fit_psf_tail   = false;
#endif
  
  int nok=0;
  int nfailed=0;
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
      nfailed++;
      
    }
    //if(int(s)%100==0 && s!=0) SPECEX_INFO("done " << s << "/" << spots.size() << " ...");
  }
  if(nfailed>0)
    SPECEX_WARNING("fit of flux of " << nfailed << " spots failed");
  
  //SPECEX_INFO("successful fit of each spot flux for " << nok << "/" << spots.size());

  force_positive_flux      = saved_force_positive_flux;

  // TURN BACK ALL MESSAGES TO REQUIRED VALUDE
  specex_set_debug(saved_debug);
  specex_set_verbose(saved_info);
  
  return true;
}

bool specex::PSF_Fitter::FitIndividualSpotPositions(std::vector<specex::Spot_p>& spots) {
  
  SPECEX_INFO("fitting independently the flux+position of each spot");
  include_signal_in_weight = false;
  fit_flux                 = true;
  fit_position             = true;
  fit_psf                  = false;
  fit_trace                = false;
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

void specex::PSF_Fitter::compare_spots_chi2_and_mask(std::vector<specex::Spot_p>& spots, const double& nsig) {
  for(size_t i=0;i<spots.size();i++) {
    specex::Spot_p spot= spots[i];
    // look at other spots with same wavelength
    double s = 0; // number of spots
    double sc = 0; // sum of spots chi2
    double sc2 = 0; // sum of spots chi2**2
    for(size_t j=0;j<spots.size();j++) {
      if(i==j) continue;
      specex::Spot_p ospot= spots[j];
      if(fabs(ospot->wavelength-spot->wavelength)>1.) continue;
      s += 1;
      sc += ospot->chi2;
      sc2 += (ospot->chi2)*(ospot->chi2);
    }
    if(s<2) continue;
    double mchi2=sc/s;
    double var=sc2/s-mchi2*mchi2;
    if(var<=0) continue;
    double rmschi2=sqrt(var);
    if(spot->chi2 > (mchi2 + nsig*rmschi2)) {
      SPECEX_WARNING("masking spot " << i << " at x=" << spot->xc << " y=" << spot->yc << " with large chi2 = " << spot->chi2 << ">" << mchi2 << "+" << nsig << "*" << rmschi2 << "=" << (mchi2 + nsig*rmschi2));
      spot->flux=0;
      spot->status=0;

      // now set weights to zero here
      {
	Stamp stamp(footprint_weight);
	SetStampLimitsFromPSF(stamp,psf,spot->xc,spot->yc);      
	for (int j=stamp.begin_j;j<stamp.end_j;j++) {
	  for (int i=stamp.begin_i;i<stamp.end_i;i++) {
	    footprint_weight(i,j)=0;
	  }
	}
      }
      {
	Stamp stamp(corefootprint);
	SetStampLimitsFromPSF(stamp,psf,spot->xc,spot->yc);      
	for (int j=stamp.begin_j;j<stamp.end_j;j++) {
	  for (int i=stamp.begin_i;i<stamp.end_i;i++) {
	    corefootprint(i,j)=0;
	  }
	}
      }
    }
  }
}

std::vector<specex::Spot_p> specex::PSF_Fitter::select_spots(std::vector<specex::Spot_p>& input_spots, double minimum_signal_to_noise, double min_wave_dist, double chi2_nsig) {
  
  // add a systematic test of chi2 of spots
  compare_spots_chi2_and_mask(input_spots,chi2_nsig);

  std::vector<specex::Spot_p> selected_spots;
  
  for(size_t s=0;s<input_spots.size();s++) {
    specex::Spot_p spot = input_spots[s];
    
    if(spot->eflux<=0 || spot->flux/spot->eflux<minimum_signal_to_noise) {
      spot->status=0;
      continue;
    }
    // check spot is in image 
    if( spot->yc<0 || spot->yc>=image.n_rows() || spot->xc<0 || spot->xc>=image.n_cols() ) {
      spot->status=0;
      continue;
    }
    // now loop on all spots to get distance
    if(min_wave_dist>0) {
      double dist=1000;
      for(size_t s2=0;s2<input_spots.size();s2++) {
	if(s==s2) continue;
	specex::Spot_p spot2 = input_spots[s2];
	if(spot2->fiber != spot->fiber) continue;
	  dist=min(dist,fabs(spot2->wavelength-spot->wavelength));
      }
      if(dist<min_wave_dist) {
	spot->status=0;
	continue;
      }
    }
    
    
    spot->status=1;
    selected_spots.push_back(spot);
    
  }
  
  SPECEX_INFO("selected " << selected_spots.size() << " spots out of " << input_spots.size() << " with S/N>" << minimum_signal_to_noise << " and min dist = " << min_wave_dist << " A");
  return selected_spots;
}


// version where we keep all spots in an array
// it is actually a bad idea because some can have very large flux and very large uncertainty
// which can pose problems for the PSF fit
/*
std::vector<specex::Spot_p> select_spots(std::vector<specex::Spot_p>& input_spots, double minimum_signal_to_noise, double min_wave_dist=0) {
  
  std::vector<specex::Spot_p> selected_spots;

  // need to either keep all or none of the spots at the same wavelength
  std::vector<specex::SpotArray> spot_arrays = find_spot_arrays(input_spots);

  for(size_t a=0;a<spot_arrays.size();a++) {
    // compute min S/N
    specex::SpotArray& spot_array = spot_arrays[a];
    double sum_snr2=0;
    int ns=0;
    for(size_t s=0;s<spot_array.size();s++) {
      specex::Spot_p spot = spot_array[s];      
      if(spot->eflux<=0) {
	continue;
      }
      sum_snr2 += specex::square(spot->flux/spot->eflux);
      ns++;
    }
    if(ns==0) continue;
    double snr=sqrt(sum_snr2/ns);

    SPECEX_INFO("select_spots : wave=" << spot_array.wavelength << " mean snr=" << snr);
    bool select_array = (snr>minimum_signal_to_noise);
    
    // now loop on all spots to get distance
    if(min_wave_dist>0) {
      for(size_t s=0;s<spot_array.size();s++) {
	specex::Spot_p spot = spot_array[s];    
	double dist=100000;
	for(size_t s2=0;s2<input_spots.size();s2++) {
	  specex::Spot_p spot2 = input_spots[s2];
	  if(spot2->fiber != spot->fiber) continue;
	  if(spot2->status==0) continue;
	  if(spot2->wavelength == spot->wavelength) continue;
	  dist=min(dist,fabs(spot2->wavelength-spot->wavelength));
	}
	if(dist<min_wave_dist) {
	  select_array=false;
	  break;
	}
      }
      
    }
    for(size_t s=0;s<spot_array.size();s++) {
      specex::Spot_p spot = spot_array[s];     
      if(select_array && spot->eflux>0) {
	spot->status=1;
	selected_spots.push_back(spot);
      }else{
	spot->status=0;
      }
    }
  } // end of loop on spot arrays
  
  SPECEX_INFO("selected " << selected_spots.size() << " spots out of " << input_spots.size() << " with S/N>" << minimum_signal_to_noise << " and min dist = " << min_wave_dist << " A");
  return selected_spots;
}
*/

bool specex::PSF_Fitter::FitEverything(std::vector<specex::Spot_p>& input_spots, bool init_psf) {

  if(input_spots.size()==0) {
    SPECEX_ERROR("in fit_several_spots spotarrays is empty");
  }

  
  
  
  
  
  double chi2=1e30;
  int npix = 0;
  int niter = 0;
  double min_snr_non_linear_terms = 10;
  double min_wave_dist_non_linear_terms = 4; // A 
  
  double min_snr_linear_terms = 3;
  double min_wave_dist_linear_terms = 0; 
  
  
  
  SPECEX_INFO("starting to fit PSF with " <<  input_spots.size() << " spots");
 
  if(init_psf) {
    
    SPECEX_INFO("init PSF");
    
    FitTraces(input_spots);

    int number_of_fibers_with_dead_columns = 0;
    
    if(1) {
      SPECEX_INFO("detecting dead columns in fiber traces ");
      for(map<int,specex::Trace>::iterator it=psf->FiberTraces.begin();
	  it !=psf->FiberTraces.end(); ++it) {
	
	if(it->first < psf_params->fiber_min || it->first > psf_params->fiber_max) continue;
	if(it->second.Off()) continue;
	
	specex::Trace& trace=it->second;
	int begin_j = max(0,int(floor(trace.Y_vs_W.Value(trace.Y_vs_W.xmin))));
	int end_j   = min(int(weight.n_rows()),int(floor(trace.Y_vs_W.Value(trace.Y_vs_W.xmax)))+1);
	int ndead=0;
	for(int j=begin_j;j<end_j;j++) {
	  int i_center = trace.X_vs_Y.Value(double(j));
	  int begin_i = max(0,i_center-3);
	  int end_i   = min(int(weight.n_cols()),i_center+4);
	  for(int i = begin_i ;i<end_i;i++) {
	    if(weight(i,j)==0) ndead++;
	  }
	} 
	if(ndead>0)
	  SPECEX_INFO("fiber " << it->first << " ndead=" << ndead);
	if(ndead>500) number_of_fibers_with_dead_columns++;
      }
      if(number_of_fibers_with_dead_columns>0)
	SPECEX_INFO("Number of fibers with dead columns = " << number_of_fibers_with_dead_columns);
    }
    
    // here we need to count number of wavelength and spots per wavelength
    {
      
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
      // can happend when testing
      if(max_x==min_x) max_x = min_x+1;
      if(max_wave==min_wave) max_wave = min_wave+1;
      
      // if the traces are already synchronized, as is the case for DESI sims,
      // we use this predefined range of coordinates
      
      bool modified=false;
      for(map<int,specex::Trace>::iterator it=psf->FiberTraces.begin();
	  it !=psf->FiberTraces.end(); ++it) {
	const Trace& trace = it->second;
	if(trace.synchronized) {
	  min_wave = min( min_wave , trace.X_vs_W.xmin);
	  max_wave = max( max_wave , trace.X_vs_W.xmax);
	  min_x = min( min_x , trace.X_vs_W.Value(trace.X_vs_W.xmin));
	  max_x = max( max_x , trace.X_vs_W.Value(trace.X_vs_W.xmax));
	  modified=true;
	}
      }
      if(modified)
	SPECEX_INFO("Use range for traces wmin,wmax=" << min_wave << "," << max_wave << " xmin,xmax=" << min_x << "," << max_x);

      
      // get max number of spots per fiber
      map<int,int> number_of_spots_per_fiber;
      for(size_t s=0;s<input_spots.size(); ++s) {
	const specex::Spot_p spot = input_spots[s];
	number_of_spots_per_fiber[spot->fiber]++;
      }
      int number_of_fibers = number_of_spots_per_fiber.size();
      int max_number_of_spots_per_fiber = 0;
      for(map<int,int>::const_iterator it=number_of_spots_per_fiber.begin();
	  it != number_of_spots_per_fiber.end(); ++it) {
	if(it->second > max_number_of_spots_per_fiber ) max_number_of_spots_per_fiber = it->second;
      }
      
      SPECEX_INFO("Number of fibers " << number_of_fibers);
      SPECEX_INFO("Max number of spots per fiber " << max_number_of_spots_per_fiber);
     
      
      if(polynomial_degree_along_x>number_of_fibers-1) {
	polynomial_degree_along_x=number_of_fibers-1;
	SPECEX_INFO("Reducing polynomial degree along x to " << polynomial_degree_along_x << " because of number of fibers");
      }
      if(polynomial_degree_along_wave>max_number_of_spots_per_fiber-1) {
	polynomial_degree_along_wave=max_number_of_spots_per_fiber-1;
	SPECEX_INFO("Reducing polynomial degree along wavelength to " << polynomial_degree_along_wave << " because of number of spots");
      }
      
      SPECEX_INFO("Setting PSF polynomial degrees " << polynomial_degree_along_x << " " << polynomial_degree_along_wave);
      
      psf_params->AllParPolXW.clear();
      
      int npar = psf->LocalNAllPar();
      
      harp::vector_double default_params = psf->DefaultParams();
      std::vector<string> param_names = psf->DefaultParamNames();
      SPECEX_INFO("Default PSF params = " <<  default_params);
      
      if(int(default_params.size()) != npar) SPECEX_ERROR("Fatal inconsistency in number of parameters of PSF");
      if(int(param_names.size()) != npar) SPECEX_ERROR("Fatal inconsistency in number of parameters of PSF");
      
      int npar_tot = 0;
      for(int p=0;p<npar;p++) {
	const string& name = param_names[p];
	
	int degx=polynomial_degree_along_x;
	int degw=polynomial_degree_along_wave;
	
	// don't want to play with second gaussian parameters if exist
	if(name=="GHSIGX2" || name=="GHSIGY2" || name=="GHSCAL2" || name=="GHNSIG") {
	  degx=0;
	  degw=0;
	}
	
	// this is to readjust fiber positions
	if(name=="GH-1-0" || name=="GH-0-1" || name=="GH-1-1" ) {
	  degx=number_of_fibers-1;
	  if(number_of_fibers_with_dead_columns>0) {
	    degx=max(0,degx-number_of_fibers_with_dead_columns);
	    // remove one more degree for precaution (not well understool why it's needed)
	    degx=max(0,degx-1);
	  }
	}
	
	// minimal variation of abberation if exist
	if(name.find("GH2")!=name.npos) {
	  degx=0;
	  degw=min(2,int(polynomial_degree_along_wave));
	}
	
	// minimal variation of tail parameters is exist
	if(name=="TAILAMP") {
	  degx=0;
	  degw=min(1,int(polynomial_degree_along_wave));
	}
	if(name=="TAILCORE" || name=="TAILXSCA" || name=="TAILYSCA" || name=="TAILINDE" ) {
	  degx=0;
	  degw=0;
	}
	
	if(reduce_psf_params_variation){
	  // loop on high order GH coefficients if exist
	  // limit to legendre deg 4 higher orders
	  char label[1000];
	  for(int gh_i=0; gh_i<12; gh_i++) { // hard coded big number
	    for(int gh_j=0; gh_j<12; gh_j++) { // hard coded big number
	      if(gh_i+gh_j<=2) continue; // we leave upto second order as it	      
	      sprintf(label,"GH-%d-%d",gh_i,gh_j);
	      if(name==label) {
		degw=min(4,degw);
	      }
	    }
	  }
	}


	specex::Pol_p pol(new specex::Pol(degx,min_x,max_x,degw,min_wave,max_wave));
	pol->name = name;
	
	
	// insert by hand
	if(reduce_psf_params_variation) {
	  for(int i=0;i<=degx;i++) { // x coordinate
	    for(int j=0;j<=degw;j++) { // wave coordinate	      
	      if(i==0) {
		pol->Add(i,j); // full wavelength resolution
	      }else if(i==1) {
		if(j<2) pol->Add(i,j); // only first x * wavelength cross-term
	      }else{
		if(j==0) pol->Add(i,j); // only  x terms
	      }	    
	    }
	  }
	}else{
	  pol->Fill();
	}
	
	

	pol->coeff(0) = default_params(p);
	SPECEX_DEBUG("Init P" << p << " " << pol->name << " =" << pol->coeff(0) << ", degw=" << degw << " degx=" << degx << " ncoef=" << pol->coeff.size());

	npar_tot += pol->coeff.size();

	psf_params->AllParPolXW.push_back(pol);
	
      }
      SPECEX_INFO("Total number of PSF params = " << npar_tot);
      
      //exit(12); // debug

#ifdef CONTINUUM
      //SPECEX_WARNING("I set deg=0 to cont and tail");
      psf_params->ContinuumPol.deg = 2;
      psf_params->ContinuumPol.coeff.resize(psf_params->ContinuumPol.deg+1);
      psf_params->ContinuumPol.coeff.clear();
      psf_params->ContinuumPol.xmin = min_wave;
      psf_params->ContinuumPol.xmax = max_wave;
#endif
    }
    
    for(map<string,Prior*>::const_iterator it=priors.begin(); it!=priors.end(); ++it) {
      SPECEX_INFO("Setting Gaussian prior on param " << it->first);
      psf->SetPrior(it->first,it->second);
    }
  }
  
  fatal = true;
  include_signal_in_weight = false;
  //include_signal_in_weight = true;
  chi2_precision = 50.;
  bool ok = true;
  
  SPECEX_INFO("Starting FitSeveralSpots FLUX");
  
  int saved_psf_hsizex = psf->hSizeX;
  int saved_psf_hsizey = psf->hSizeY;
  psf->hSizeX=min(3,psf->hSizeX);
  psf->hSizeY=min(3,psf->hSizeY);
  
  include_signal_in_weight = false;
  //include_signal_in_weight = true;
  ok = FitIndividualSpotFluxes(input_spots);
  
  std::vector<specex::Spot_p> selected_spots = select_spots(input_spots,min_snr_non_linear_terms,min_wave_dist_non_linear_terms);
  if(write_tmp_results) {
    write_spots_xml(selected_spots,"spots-before-gaussian-fit.xml");
    write_psf_xml(psf,"psf-before-gaussian-fit.xml");
  }

  
  if(scheduled_fit_of_traces) {
    
    // reduce trace degree if not enough spots to fit
    // --------------------------------------------  
    for(int fiber=psf_params->fiber_min; fiber<=psf_params->fiber_max; fiber++) 
      {
	int nspots=0;
	for(size_t s=0;s<selected_spots.size();s++) {
	  Spot_p spot= selected_spots[s];
	  if(spot->fiber==fiber) nspots++;
	}
	
	
	
	specex::Trace& trace = psf->FiberTraces.find(fiber)->second;
	if(trace.Off()) continue;
	
	SPECEX_INFO("Fiber " << fiber << " nspots= " << nspots 
		    << " xdeg=" << trace.X_vs_W.deg
		    << " ydeg=" << trace.Y_vs_W.deg);

	if(trace.X_vs_W.deg > (nspots-1) || trace.Y_vs_W.deg > (nspots-1) 
	   || trace.W_vs_Y.deg > (nspots-1) || trace.X_vs_Y.deg > (nspots-1) ) {
	  SPECEX_WARNING("Reducing degree of trace of fiber "<< fiber << " to match number of spots = " << nspots);
	  trace.resize(nspots);
	}

	// check current offset wrt initial position
	double sum=0;
	double sdx=0;
	double sdy=0;
	for(size_t s=0;s<selected_spots.size();s++) {
	  Spot_p spot= selected_spots[s];
	  if(spot->fiber==fiber) {
	    sum++;
	    sdx+= (spot->xc - spot->initial_xc);
	    sdy+= (spot->yc - spot->initial_yc);	  
	  }
	}
	if(sum>0) {
	  double mean_dx=sdx/sum;
	  double mean_dy=sdy/sum;
	  SPECEX_INFO("Fiber " << fiber << " dx=" << mean_dx << " dy=" << mean_dy);
	}else{
	  SPECEX_WARNING("No selected spot for fiber " << fiber);
	  trace.mask=3;
	}

      }
    

    SPECEX_INFO("Starting FitSeveralSpots TRACE");

    fit_flux       = false;
    fit_position   = false;
    fit_psf        = false;
    fit_trace      = true;
    ok = FitSeveralSpots(selected_spots,&chi2,&npix,&niter);
    if(!ok) SPECEX_ERROR("FitSeveralSpots failed for TRACE");
    
    
     
    { // testing traces
      vector<int> fibers_with_large_offsets;
      map<int,bool> fiber_is_ok;      
      for(int fiber = psf_params->fiber_min; fiber<=psf_params->fiber_max;fiber++) {
	int nbad=0;
	double delta=2; //pix
	double max_dx=0;
	double max_dy=0;
	for(size_t s=0;s<input_spots.size();s++) {
	  specex::Spot_p spot= input_spots[s];
	  if(spot->fiber != fiber) continue;
	  max_dx=max(max_dx,fabs(psf->Xccd(spot->fiber,spot->wavelength)-spot->initial_xc));
	  max_dy=max(max_dy,fabs(psf->Yccd(spot->fiber,spot->wavelength)-spot->initial_yc));	  
	  if(fabs(psf->Xccd(spot->fiber,spot->wavelength)-spot->initial_xc)>delta || fabs(psf->Yccd(spot->fiber,spot->wavelength)-spot->initial_yc)>delta) nbad++;
	}
	if(nbad>2) {
	  SPECEX_WARNING("Large x y offset for fiber " << fiber << " max dx,dy=" << max_dx << "," << max_dy);
	  fibers_with_large_offsets.push_back(fiber);
	  fiber_is_ok[fiber]=false;
	}else{
	  fiber_is_ok[fiber]=true;
	}
	
      }
      
      if(fibers_with_large_offsets.size()>=3) {
	SPECEX_WARNING("There are more than 2 fibers with large offsets, this is not due to dead columns, so we continue");
      }
      
      if(false && fibers_with_large_offsets.size()>0 && fibers_with_large_offsets.size()<3) {
	// try to fix this : use interpolation of other fibers : this assume fiber slit heads allow it, as for BOSS
	for(size_t f=0;f<fibers_with_large_offsets.size();f++) {
	  int bad_fiber   = fibers_with_large_offsets[f];
	  int fiber1  = -1;
	  int fiber2  = -1;
	  
	  if(bad_fiber==psf_params->fiber_min) {
	    for(int fiber=bad_fiber+1;fiber<=psf_params->fiber_max;fiber++)
	      if(fiber_is_ok[fiber]) { fiber1=fiber; break;}
	  }else{
	    for(int fiber=bad_fiber-1;fiber>=psf_params->fiber_min;fiber--)
	      if(fiber_is_ok[fiber]) { fiber1=fiber; break;}
	  }
	  if(bad_fiber==psf_params->fiber_max) {
	    for(int fiber=min(bad_fiber,fiber1)-1;fiber>=psf_params->fiber_min;fiber--)
	      if(fiber_is_ok[fiber]) { fiber2=fiber; break;}
	  }else{
	    for(int fiber=max(bad_fiber,fiber1)+1;fiber<=psf_params->fiber_max;fiber++)
	      if(fiber_is_ok[fiber]) { fiber2=fiber; break;}
	  }
	  if(fiber1==-1 || fiber2==-1) {
	    SPECEX_ERROR("can't fix this, abort");
	  }
	  SPECEX_INFO("trying to fix bad fiber trace " << bad_fiber << " with " << fiber1 << " and " << fiber2);
	  specex::Trace &bad_trace = psf->FiberTraces.find(bad_fiber)->second;
	  const specex::Trace &trace1 = psf->FiberTraces.find(fiber1)->second;
	  const specex::Trace &trace2 = psf->FiberTraces.find(fiber2)->second;
	  bad_trace.X_vs_W.coeff = (float(fiber2-bad_fiber)/(fiber2-fiber1))*trace1.X_vs_W.coeff + (float(bad_fiber-fiber1)/(fiber2-fiber1))*trace2.X_vs_W.coeff;
	  bad_trace.Y_vs_W.coeff = (float(fiber2-bad_fiber)/(fiber2-fiber1))*trace1.Y_vs_W.coeff + (float(bad_fiber-fiber1)/(fiber2-fiber1))*trace2.Y_vs_W.coeff;
	  
	}
	scheduled_fit_of_traces = false;
      }
      
    
      for(size_t s=0;s<input_spots.size();s++) {
	specex::Spot_p& spot= input_spots[s];
	spot->xc = psf->Xccd(spot->fiber,spot->wavelength);
	spot->yc = psf->Yccd(spot->fiber,spot->wavelength);
      }
    }
    
    if(write_tmp_results)
      write_spots_xml(input_spots,"spots-after-trace-fit.xml");
  }
  
  if(scheduled_fit_of_traces) {
    
    for(size_t s=0;s<input_spots.size();s++) {
      specex::Spot_p& spot= input_spots[s];
      spot->xc = psf->Xccd(spot->fiber,spot->wavelength);
      spot->yc = psf->Yccd(spot->fiber,spot->wavelength);
    }
    include_signal_in_weight = false;
    //include_signal_in_weight = true;
    ok = FitIndividualSpotFluxes(input_spots);
    selected_spots = select_spots(input_spots,min_snr_non_linear_terms,min_wave_dist_non_linear_terms);
  



    SPECEX_INFO("Starting FitSeveralSpots FLUX+TRACE ");
    SPECEX_INFO("=======================================");
    fit_flux       = true;
    fit_position   = false;
    fit_psf        = false;
    fit_trace      = true;
    ok = FitSeveralSpots(selected_spots,&chi2,&npix,&niter);
    if(!ok) SPECEX_ERROR("FitSeveralSpots failed for FLUX+TRACE");
  
    for(size_t s=0;s<input_spots.size();s++) {
      specex::Spot_p& spot= input_spots[s];
      spot->xc = psf->Xccd(spot->fiber,spot->wavelength);
      spot->yc = psf->Yccd(spot->fiber,spot->wavelength);
    }
    include_signal_in_weight = false;
    //include_signal_in_weight = true;
    ok = FitIndividualSpotFluxes(input_spots);
    selected_spots = select_spots(input_spots,min_snr_non_linear_terms,min_wave_dist_non_linear_terms);
    
    /*
    if(scheduled_fit_of_sigmas) {
      SPECEX_INFO("Starting FitSeveralSpots PSF+FLUX+TRACE only gaussian terms");
      SPECEX_INFO("========================================================");
      fit_flux       = true;
      fit_position   = false;
      fit_psf        = true;
      fit_trace      = true;
      ok = FitSeveralSpots(selected_spots,&chi2,&npix,&niter);
      if(!ok) SPECEX_ERROR("FitSeveralSpots failed for PSF+FLUX+TRACE only gaussian terms");
      
    }
    */
  }
  
  if(scheduled_fit_of_sigmas) {
    {
      SPECEX_INFO("Choose the parameters that participate to the fit : only gaussian terms");
      psf_params->FitParPolXW.clear();
      int npar = psf->LocalNAllPar();
      for(int p=0;p<npar;p++) {
	const string& name = psf->ParamName(p);
	bool ok = false;
	ok |= (name=="GHSIGX" || name=="GHSIGY" || name=="RADIUS"  || name=="SIGMA" );
	//ok |= (name=="GHSIGX2");
	//ok |= (name=="GHSIGY2");
	//ok |= (name=="GHSCAL2");
	if(ok)
	  psf_params->FitParPolXW.push_back(psf_params->AllParPolXW[p]);
      }
    }
    chi2=1e30;

    for(int i=0;i<5;i++) {
      SPECEX_INFO("Starting FitSeveralSpots PSF gaussian terms then FLUX (loop="<<i<<")");
      
      fit_flux       = false;
      fit_position   = false;
      fit_psf        = true;
      fit_trace      = false;
      float previous_chi2 = chi2;
      ok = FitSeveralSpots(selected_spots,&chi2,&npix,&niter);
      if(!ok) SPECEX_ERROR("FitSeveralSpots failed for PSF");
      
      ok = FitIndividualSpotFluxes(input_spots);
      selected_spots = select_spots(input_spots,min_snr_non_linear_terms,min_wave_dist_non_linear_terms);
      if(fabs(previous_chi2 - chi2)<chi2_precision) break;
    }
    chi2_precision = 0.1;
    force_positive_flux = true;
    increase_weight_of_side_bands = false;
  
    SPECEX_INFO("Starting FitSeveralSpots PSF+FLUX only gaussian terms");

    fit_flux       = true;
    fit_position   = false;
    fit_psf        = true;
    fit_trace      = false;
    ok = FitSeveralSpots(selected_spots,&chi2,&npix,&niter);
    if(!ok) SPECEX_ERROR("FitSeveralSpots failed for PSF+FLUX");
    
    
    /* dump central values for all fibers */

    for(int fiber=psf_params->fiber_min; fiber<=psf_params->fiber_max; fiber++) {
      specex::Trace& trace = psf->FiberTraces.find(fiber)->second;
      double w1 = trace.X_vs_W.xmin;
      double w2 = (trace.X_vs_W.xmin+trace.X_vs_W.xmax)/2;
      double w3 = trace.X_vs_W.xmax;
    
      SPECEX_INFO("fiber " << fiber 
		  << " wave " << w1 << " :" << ublas::project(psf->AllLocalParamsFW(fiber,w1,psf_params->bundle_id),ublas::range(0,2)) 
		  << " wave " << w2 << " :" << ublas::project(psf->AllLocalParamsFW(fiber,w2,psf_params->bundle_id),ublas::range(0,2))
		<< " wave " << w3 << " :" << ublas::project(psf->AllLocalParamsFW(fiber,w3,psf_params->bundle_id),ublas::range(0,2)));
      
    }
    if(write_tmp_results) {
      write_spots_xml(selected_spots,"spots-after-gaussian-fit.xml");
      write_psf_xml(psf,"psf-after-gaussian-fit.xml");
    }
  }
    
  if(psf->HasParam("GHNSIG")) {
    double inner_core_radius_n_sigma = 3.5;
    psf_params->AllParPolXW[psf->ParamIndex("GHNSIG")]->coeff(0)=inner_core_radius_n_sigma;
    SPECEX_INFO("Setting PSF inner core radius to " << inner_core_radius_n_sigma);
  }

  // refit individual fluxes
  for(size_t s=0;s<input_spots.size();s++) {
    specex::Spot_p& spot= input_spots[s];
    spot->xc = psf->Xccd(spot->fiber,spot->wavelength);
    spot->yc = psf->Yccd(spot->fiber,spot->wavelength);
  }
  include_signal_in_weight = false;
  //include_signal_in_weight = true;
  
  ok = FitIndividualSpotFluxes(input_spots);
  
  
  selected_spots = select_spots(input_spots,min_snr_linear_terms,min_wave_dist_linear_terms);
  psf->hSizeX=saved_psf_hsizex;
  psf->hSizeY=saved_psf_hsizey;

#ifdef CONTINUUM
#ifdef EXTERNAL_TAIL  

  if(scheduled_fit_of_psf_tail || scheduled_fit_of_continuum) {

    {
      psf_params->FitParPolXW.clear();
      for(int p=0;p<psf->LocalNAllPar();p++) {
	const string& name = psf->ParamName(p);
	ok = (name=="TAILAMP");
	if(ok)
	  psf_params->FitParPolXW.push_back(psf_params->AllParPolXW[p]);
      }
    }
    
    SPECEX_INFO("Starting FitSeveralSpots TAIL&CONTINUUM");
    SPECEX_INFO("=======================================");
    fit_flux       = false;
    fit_position   = false;
    fit_psf        = false;
    fit_trace      = false;
    fit_psf_tail   = scheduled_fit_of_psf_tail;
    fit_continuum  = scheduled_fit_of_continuum;
    ok = FitSeveralSpots(selected_spots,&chi2,&npix,&niter);
    fit_psf_tail   = false; // don't fit this anymore
    fit_continuum   = false; // don't fit this anymore
    {
      psf_params->FitParPolXW.clear();
      for(int p=0;p<psf->LocalNAllPar();p++) {
	const string& name = psf->ParamName(p);
	ok = true;
	ok &= (name!="GHSIGX" && name!="GHSIGY" && name!="GHNSIG");
	ok &= (name!="GHSIGX2" && name!="GHSIGY2" && name!="GHSCAL2");
	//ok &= (name!="RADIUS" && name!="SIGMA");
	ok &= (name!="RADIUS");
	ok &= (name!="TAILAMP" && name!="TAILCORE" && name!="TAILXSCA" && name!="TAILYSCA" && name!="TAILINDE");
	if(ok)
	  psf_params->FitParPolXW.push_back(psf_params->AllParPolXW[p]);
      }
    }
  }
#endif
#endif
  
  {
    psf_params->FitParPolXW.clear();
    for(int p=0;p<psf->LocalNAllPar();p++) {
      const string& name = psf->ParamName(p);
      ok = true;
      ok &= (name!="GHSIGX" && name!="GHSIGY" && name!="GHNSIG");
      ok &= (name!="GHSIGX2" && name!="GHSIGY2" && name!="GHSCAL2");
      ok &= (name!="RADIUS" && name!="SIGMA");
      ok &= (name!="TAILAMP" && name!="TAILCORE" && name!="TAILXSCA" && name!="TAILYSCA" && name!="TAILINDE");
      if(ok)
	psf_params->FitParPolXW.push_back(psf_params->AllParPolXW[p]);
    }
  }
  
  fit_flux       = false;
  fit_position   = false;
  fit_trace      = false;
  include_signal_in_weight = false;
  //include_signal_in_weight = true;
  
  int count=1;

  /*
  fit_flux = true; fit_psf = false;
  SPECEX_INFO("Starting FitSeveralSpots FLUX #" << count);
  ok = FitSeveralSpots(selected_spots,&chi2,&npix,&niter);
  selected_spots = select_spots(selected_spots,min_snr_linear_terms,min_wave_dist_linear_terms);
  
  if(write_tmp_results) {
    char filename[1000];
    sprintf(filename,"spots-after-flux-%d.xml",count);
    write_spots_xml(selected_spots,filename);
  }
  */
  if(scheduled_fit_of_psf) {
    fit_flux = false; fit_psf = true;
    SPECEX_INFO("Starting FitSeveralSpots PSF #" << count);
    ok = FitSeveralSpots(selected_spots,&chi2,&npix,&niter);
    fit_flux = true; fit_psf = true;
    SPECEX_INFO("Starting FitSeveralSpots PSF+FLUX #" << count);
    ok = FitSeveralSpots(selected_spots,&chi2,&npix,&niter);
    if(!ok) SPECEX_ERROR("FitSeveralSpots failed for PSF+FLUX");
    if(write_tmp_results) {
      char filename[1000];
      sprintf(filename,"spots-after-psf+flux-%d.xml",count);
      write_spots_xml(selected_spots,filename);
    }
    
    
    // if we have continuum and tail, additional loop
#ifdef CONTINUUM
#ifdef EXTERNAL_TAIL  
    
    if(scheduled_fit_of_psf_tail || scheduled_fit_of_continuum) {
      psf_params->FitParPolXW.clear();
      for(int p=0;p<psf->LocalNAllPar();p++) {
	const string& name = psf->ParamName(p);
	ok = (name=="TAILAMP");
	if(ok)
	  psf_params->FitParPolXW.push_back(psf_params->AllParPolXW[p]);
      }
  
      
      SPECEX_INFO("Starting FitSeveralSpots TAIL&CONTINUUM");
      fit_flux       = false;
      fit_position   = false;
      fit_psf        = false;
      fit_trace      = false;
      fit_psf_tail   = scheduled_fit_of_psf_tail;
      fit_continuum  = scheduled_fit_of_continuum;
      ok = FitSeveralSpots(selected_spots,&chi2,&npix,&niter);
      
      
      fit_psf_tail   = false; // don't fit this anymore
      fit_continuum   = false; // don't fit this anymore
      {
	psf_params->FitParPolXW.clear();
	for(int p=0;p<psf->LocalNAllPar();p++) {
	  const string& name = psf->ParamName(p);
	  ok = true;
	  ok &= (name!="GHSIGX" && name!="GHSIGY" && name!="GHNSIG");
	  ok &= (name!="GHSIGX2" && name!="GHSIGY2" && name!="GHSCAL2");
	  //ok &= (name!="RADIUS" && name!="SIGMA");
	  ok &= (name!="RADIUS");
	  ok &= (name!="TAILAMP" && name!="TAILCORE" && name!="TAILXSCA" && name!="TAILYSCA" && name!="TAILINDE");
	  if(ok)
	    psf_params->FitParPolXW.push_back(psf_params->AllParPolXW[p]);
	}
      }
      
      fit_flux = true; fit_psf = false;
      count++;
      
      /*
	SPECEX_INFO("Starting FitSeveralSpots FLUX #" << count);
	ok = FitSeveralSpots(selected_spots,&chi2,&npix,&niter);
	selected_spots = select_spots(selected_spots,min_snr_linear_terms,min_wave_dist_linear_terms);
	
	if(write_tmp_results) {
	char filename[1000];
	sprintf(filename,"spots-after-flux-%d.xml",count);
	write_spots_xml(selected_spots,filename);
	}
      */
      
      fit_flux = false; fit_psf = true;
      SPECEX_INFO("Starting FitSeveralSpots PSF #" << count);
      ok = FitSeveralSpots(selected_spots,&chi2,&npix,&niter);
      
    
      fit_flux = true; fit_psf = true;
      SPECEX_INFO("Starting FitSeveralSpots PSF+FLUX #" << count);
      ok = FitSeveralSpots(selected_spots,&chi2,&npix,&niter);
      
      if(!ok) SPECEX_ERROR("FitSeveralSpots failed for PSF+FLUX");
      if(write_tmp_results) {
	char filename[1000];
	sprintf(filename,"spots-after-psf+flux-%d.xml",count);
	write_spots_xml(selected_spots,filename);
      }
      

    } // end of test of fit of tail or continuum
  
#endif
#endif
  
    if(scheduled_fit_with_weight_model)
      { 
	// use model instead of data to weight the objects
	// to avoid biases on the PSF
	include_signal_in_weight = true;
	fit_flux = true; fit_psf = false;
	count++;
	SPECEX_INFO("Starting FitSeveralSpots FLUX(w) #" << count);
	ok = FitSeveralSpots(selected_spots,&chi2,&npix,&niter);
	selected_spots = select_spots(selected_spots,min_snr_linear_terms,min_wave_dist_linear_terms);
	
	if(write_tmp_results) {
	  char filename[1000];
	  sprintf(filename,"spots-after-flux-%d.xml",count);
	  write_spots_xml(selected_spots,filename);
	}
	
	fit_flux = false; fit_psf = true;
	SPECEX_INFO("Starting FitSeveralSpots PSF(w) #" << count);
	ok = FitSeveralSpots(selected_spots,&chi2,&npix,&niter);
	
	if(write_tmp_results) {
	  char filename[1000];
	  sprintf(filename,"spots-after-psf-%d.xml",count);
	  write_spots_xml(selected_spots,filename);
	}
      }
  } // end of test on scheduled_fit_of_psf
  
  SPECEX_DEBUG("Compute in-core chi2");
  
  increase_weight_of_side_bands = false;
  
  fit_flux       = false;
  fit_position   = false;
  fit_psf        = false;
  fit_trace      = false;
#ifdef CONTINUUM
  fit_continuum  = false;
#endif
#ifdef EXTERNAL_TAIL
  fit_psf_tail   = false;
#endif  

  int saved_hsizex = psf->hSizeX;
  int saved_hsizey = psf->hSizeY;
  psf->hSizeX=2;
  psf->hSizeY=2;
  if(scheduled_fit_with_weight_model)
    include_signal_in_weight = true;
  ComputeWeigthImage(selected_spots,&npix);  
  psf_params->ndata_in_core = npix;
  InitTmpData(selected_spots);

  for(size_t s=0;s<spot_tmp_data.size();s++) 
    if(spot_tmp_data[s].flux<0) spot_tmp_data[s].flux=0;
  
  psf_params->chi2_in_core = ParallelizedComputeChi2AB(false);
  psf->hSizeX=saved_hsizex;
  psf->hSizeY=saved_hsizey;

 
  SPECEX_DEBUG("Compute final chi2");
  
  fit_flux       = false;
  fit_position   = false;
  fit_psf        = false;
  fit_trace      = false;
#ifdef CONTINUUM
  fit_continuum  = false;
#endif
#ifdef EXTERNAL_TAIL
  fit_psf_tail   = false;
#endif  
  if(scheduled_fit_with_weight_model)
    include_signal_in_weight = true;
  ComputeWeigthImage(selected_spots,&npix); 
  psf_params->ndata = npix; 
  InitTmpData(selected_spots);
  
  for(size_t s=0;s<spot_tmp_data.size();s++) 
    if(spot_tmp_data[s].flux<0) spot_tmp_data[s].flux=0;
  
  psf_params->chi2 = ParallelizedComputeChi2AB(false);
  return ok;
}
