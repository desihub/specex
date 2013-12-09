#include <assert.h>

#include "harp.hpp"

#include "specex_base_analytic_psf.h"
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

void specex::PSF_Fitter::SetStampLimitsFromPSF(specex::Stamp& stamp, const specex::PSF& PSF, const double &X, const double &Y) {
  psf.StampLimits(X,Y,stamp.begin_i,stamp.end_i,stamp.begin_j,stamp.end_j);
  // now check image bounds
  stamp.begin_i = max(0,stamp.begin_i);
  stamp.end_i   = min(stamp.Parent_n_cols(),stamp.end_i);
  stamp.begin_j = max(0,stamp.begin_j);
  stamp.end_j   = min(stamp.Parent_n_rows(),stamp.end_j);
}
 
void specex::PSF_Fitter::SetStampLimitsFromPSF(specex::Stamp& stamp, const specex::PSF& psf, const double &xc_min, const double &xc_max, const double &yc_min, const double &yc_max) {
  int k,p;
  psf.StampLimits(xc_min,yc_min,stamp.begin_i,k,stamp.begin_j,p);
  psf.StampLimits(xc_max,yc_max,k,stamp.end_i,p,stamp.end_j);
  
  // now check image bounds
  stamp.begin_i = max(0,stamp.begin_i);
  stamp.end_i   = min(stamp.Parent_n_cols(),stamp.end_i);
  stamp.begin_j = max(0,stamp.begin_j);
  stamp.end_j   = min(stamp.Parent_n_rows(),stamp.end_j);
}

struct BrentBox { // handler to data and params needed to compute chi2 in brent routine
  specex::PSF_Fitter& fitter;
  harp::vector_double& delta_P;
  vector<specex::Spot*>& spots;
  BrentBox(specex::PSF_Fitter& i_fitter,  harp::vector_double& i_delta_P, vector<specex::Spot*>& i_spots) :
    fitter(i_fitter),
    delta_P(i_delta_P),
    spots(i_spots)
  {}
};



double compute_chi2_for_a_given_step(const double &current_step, BrentBox* bbox) {
  bbox->fitter.Params += current_step*bbox->delta_P;
  double chi2 = bbox->fitter.ComputeChi2AB(bbox->spots,false);
  bbox->fitter.Params -= current_step*bbox->delta_P; // go back after testing
  if(bbox->fitter.verbose) {
    SPECEX_INFO("brent step=" << current_step << " chi2=" << chi2);
  }
  return chi2;
}



int specex::PSF_Fitter::NPar(int nspots) const {
  int npar = 0; // fluxes
  if(fit_psf) npar += psf.VaryingCoordNPar();
  if(fit_trace) npar += psf.TracesNPar();
  if(fit_flux) npar += nspots;
  if(fit_position) npar += 2*nspots;
  return npar;
}


double specex::PSF_Fitter::ComputeChi2AB(vector<specex::Spot*>& spots, bool compute_ab) 
{
  // The idea here is that the layout of A and B does not depend on
  // what is to be fitted. So, we manipulate "big" matrices (~ 6x6) even if 
  // only fitting pos and flux. 

  if(fit_trace && fit_position) {
    SPECEX_ERROR("specex::PSF_Fitter::ComputeChi2AB cannot fit traces and spot positions at the same time");
  }

  vector<harp::vector_double> psf_Params_of_spots;
  vector<harp::vector_double> psf_Monomials_of_spots; // still compute it if pos is varying
  
  {
    int index = 0;
    if(fit_psf)   index += npar_varying_coord;
    if(fit_trace) index += npar_traces;
    
    for(int s=0;s<nspots;s++) {
      const specex::Spot& spot=*spots[s];
      if(fit_flux) 
	spots_flux(s) = Params(index++);

      double& xs =spots_x(s);
      double& ys =spots_y(s);
      
      if(fit_trace) {
	specex::Trace &trace = psf.FiberTraces[spot.fiber];
	{
	  const harp::vector_double& M=TraceXvsW_Monomials_of_spots[s];
	  size_t index = XvsW_index_of_fiber[spot.fiber];
	  xs = specex::dot(ublas::project(Params,ublas::range(index,index+M.size())),M);
	  //const double *p=&Params.Data()[XvsW_index_of_fiber[spot.fiber]];
	  //const double *m=M.Data();
	  //xs=0;
	  //for(int k=0;k<M.size();k++,p++,m++)
	  //xs += (*p)*(*m);
	}
	{
	  const harp::vector_double& M = TraceYvsW_Monomials_of_spots[s];	  
	  size_t index = YvsW_index_of_fiber[spot.fiber];
	  ys = specex::dot(ublas::project(Params,ublas::range(index,index+M.size())),M);
	  //const double *p=&Params.Data()[YvsW_index_of_fiber[spot.fiber]];
	  //const double *m=M.Data();	  
	  //ys=0;
	  //for(int k=0;k<M.size();k++,p++,m++)
	  //ys += (*p)*(*m);
	}
	
	// if(s==0) cout << "DEBUG in ComputeChi2AB coords of first spot = " << x(s) << " " << y(s) << endl;
	
	
      }else if(fit_position) {
	xs = Params(index++);
	ys = Params(index++);
      }
      
      if(fit_psf)
	psf_Params_of_spots.push_back(psf.FixedCoordParams(xs,ys,Params));
      else
	psf_Params_of_spots.push_back(spot.PSFParams);
      
      if(compute_ab && fit_psf) {
	
	harp::vector_double gM(npar_varying_coord); 
	
	//double *mglob = M.NonConstData();
	int index=0;
	
	for(int p=0;p<npar_fixed_coord;p++) {

	  harp::vector_double pM = psf.Params[p].Monomials(xs,ys);
	  size_t pM_size = pM.size();
	  ublas::project(gM,ublas::range(index,index+pM_size))=pM;
	  index += pM_size;
	  
	  //const double *mloc = pM.Data();
	  //for(size_t c=0;c<pM.size();c++,mglob++,mloc++)
	  //*mglob=*mloc;
	}
	
	psf_Monomials_of_spots.push_back(gM);
      }
    }
  }
  
  
  int npar_per_spot = 0;
  if(fit_flux) npar_per_spot++;
  if(fit_position) npar_per_spot+=2;

  
  

  harp::vector_double H;
  if(compute_ab) {
    H.resize(nparTot);
  }
  
  double chi2 = 0;
  
  harp::matrix_double &Anc  = const_cast<specex::PSF_Fitter*>(this)->A;
  harp::vector_double &Bnc = const_cast<specex::PSF_Fitter*>(this)->B;
  
  if(compute_ab) {
    Anc *= 0;
    Bnc *= 0;
  }
  bool use_footprint = (nspots>1);
  use_footprint &= (footprint_weight.Nx()>0);
  
  // cout << "DEBUG stamp " <<  stamp << endl;
  // cout << "DEBUG nspots " << nspots  << "  use_footprint = " << use_footprint << endl;
  
  for (int j=stamp.begin_j; j <stamp.end_j; ++j) {
    // if(verbose)  if (j%100==0) cout << "INFO specex::PSF_Fitter::ComputeChi2AB doing j=" << j << " in [" << stamp.begin_j << "," << stamp.end_j << "]" << endl;
    
    for (int i=stamp.begin_i ; i < stamp.end_i; ++i) {
      
     
      
      double w = 1;
      if(use_footprint) { 
	w=footprint_weight(i,j); 
      } else{
	w=weight(i,j);
	//cout << "DEBUG16 w = " << w << " " << i << " " << j << endl;
      }
      if (w<=0) continue;
      
      double res = double(image(i,j));
      
      if(compute_ab)
	H *= 0;
      
      int nspots_in_pix = 0;
      for(int s=0;s<nspots;s++) {
	
	specex::Spot &spot = *(spots[s]);
	
	if( ! spot_stamps[s].Contains(i,j)) {
	  continue;
	}
	
	
	nspots_in_pix++;
	
	const double& x_s=spots_x(s);
	const double& y_s=spots_y(s);
	const double& flux_s=spots_flux(s);
	
	double psfVal = 0;
	
	if(compute_ab) {
	  psfVal = psf.PSFValueWithParams(x_s,y_s, i, j, psf_Params_of_spots[s], 
					  (fit_position || fit_trace) ? &gradPos : NULL, 
					  fit_psf ? &gradPar : NULL);
	}else{
	  psfVal = psf.PSFValueWithParams(x_s,y_s, i, j, psf_Params_of_spots[s],NULL,NULL);
	}
	
	// cout << "DEBUG i j psf " << i << " " << j << " " << psfVal << endl;

	if(psfVal==PSF_NAN_VALUE) return 1.e20;

	res-= flux_s*psfVal;
	
	if (compute_ab) {

	  int index=0;

	  if(fit_psf) {

	    const harp::vector_double& sM = psf_Monomials_of_spots[s];
	    for(int p=0;p<npar_fixed_coord;p++) {
	      double flux_times_dpsfdpar = flux_s*gradPar[p]; // this depends on xyij coordinate;
	      size_t coef_size = psf.Params[p].coeff.size();
	      ublas::project(H,ublas::range(index,index+coef_size)) += flux_times_dpsfdpar*ublas::project(sM,ublas::range(index,index+coef_size));
	      index += coef_size;
	    }
	  }
	  if(fit_trace) {
	    {
	      double dvdx = gradPos(0) * spots_flux(s);
	      int tindex = XvsW_index_of_fiber[spot.fiber];
	      const harp::vector_double& M =  TraceXvsW_Monomials_of_spots[s];
	      ublas::project(H,ublas::range(tindex,tindex+M.size())) += dvdx*M;
	    }
	    {
	      double dvdy = gradPos(1) * spots_flux(s);
	      int tindex = YvsW_index_of_fiber[spot.fiber];
	      const harp::vector_double& M =  TraceYvsW_Monomials_of_spots[s];
	      ublas::project(H,ublas::range(tindex,tindex+M.size())) += dvdy*M;
	    }
	    index += npar_traces;
	  }
	  
	  index+=npar_per_spot*s; // we are fitting here spot number s
	  
	  if(fit_flux) 
	    H[index++] += psfVal;
	  
	  if(fit_position) {
	    H[index++] += gradPos(0) * flux_s;
	    H[index++] += gradPos(1) * flux_s;
	  }
	} // end of test compute_ab
      } // end of loop on spots
      
      
      chi2 += w*res*res;
      
      if (compute_ab) { 	
	specex::syr(w,H,A); //  A += w*h*h.transposed();
	specex::axpy(w*res,H,B); // B += w*res*h;
      }
      
    } // end of loop on pix coord. i
  } // end of loop on pix coord. j
  
  //cout << "DEBUG12 chi2 = " << chi2 << endl;

  // psf priors
  if(fit_psf && !(psf.Priors.empty())) {

    for(int s=0;s<nspots;s++) {
      const harp::vector_double& SpotParams=psf_Params_of_spots[s];
      
      
      int index=0;
      for(int p=0;p<npar_fixed_coord;p++) {
	
	const double& par = SpotParams(p);
	
	int ncoef=psf.Params[p].coeff.size();
	std::map<int,Prior*>::const_iterator it = psf.Priors.find(p);
	
	if(it==psf.Priors.end()) {index += ncoef; continue;} // no prior for this psf parameter
	const Prior* prior = it->second;
	
	if (compute_ab) {
	  
	  const harp::vector_double& sM = psf_Monomials_of_spots[s];
	  
	  for (int c=0; c<ncoef; c++, index++) {
	    
	    const double& monomial_val = sM(index);
	    Bnc(index)       += monomial_val * prior->hdChi2dx(par);
	    Anc(index,index) += square(monomial_val) * prior->hd2Chi2dx2(par);
	    
	    //cout << "spot=" << s << " p=" << p << " c=" << c << " index=" << index 
	    //     << " monomial=" << monomial_val << " dchi2=" << monomial_val * prior->hdChi2dx(par) << endl;
	    
	  }
	}
	//cout << "spot=" << s << " p=" << p << " chi2=" << prior->Chi2(par) << endl;
	chi2 += prior->Chi2(par);
      }
    
    //cout << "INFO adding prior to param(" << k << ")=" << pk << " prior chi2 = " << it->second->Chi2(pk) << endl;
      
    }

  }
  
  return chi2;
}

   
static double sign(const double& a, const double& b) {
  if(b>0) return fabs(a);
  return -fabs(a);
}




bool specex::PSF_Fitter::FitSeveralSpots(vector<specex::Spot*>& spots, double *chi2, int *n_pixels, int *n_iterations) {
    
  double chi2_memory_slot = 0;
  double *psfChi2 = &chi2_memory_slot;
  if(chi2) psfChi2 = chi2;
  int niter_memory_slot = 0;
  int *niter = &niter_memory_slot;
  if(n_iterations) niter = n_iterations;
  int npixels_memory_slot = 0;
  int *npix = &npixels_memory_slot;
  if(n_pixels) npix = n_pixels;
  
  
  double xc_min=1e20;
  double xc_max=-1e20;
  double yc_min=1e20;
  double yc_max=-1e20;
  for(size_t s=0;s<spots.size();++s) {
    const specex::Spot &spot = *(spots[s]);
    if(spot.xc<xc_min) xc_min=spot.xc;
    if(spot.xc>xc_max) xc_max=spot.xc;
    if(spot.yc<yc_min) yc_min=spot.yc;
    if(spot.yc>yc_max) yc_max=spot.yc;
  }
  SetStampLimitsFromPSF(stamp,psf,xc_min,xc_max,yc_min,yc_max);
  
  // const specex::Spot &spot = *(spots[0]); cout << "DEBUG13 " << spot.xc << " " << spot.yc << " " << stamp.begin_i << " " << stamp.end_i << " " << stamp.begin_j << " " << stamp.end_j << endl;

  if(spots.size()>1) {
    // remove boundaries of stamp to avoid contamination of spots not included in the fit
    stamp.begin_i += (psf.hSizeX-min(psf.hSizeX,5)); // leave 5 pix
    stamp.end_i   -= (psf.hSizeX-min(psf.hSizeX,5)); // leave 5 pix
    //stamp.begin_j += max(0,psf.hSizeY-2);
    //stamp.end_j   -= max(0,psf.hSizeY-2);
  }
  
  
  


  int npar_psf = 0;
  if(fit_psf) npar_psf = psf.VaryingCoordNPar();
  int npar_trace = 0;
  if(fit_trace) npar_trace = psf.TracesNPar();
  
  nparTot  = NPar(spots.size());
  
  if(verbose && specex_verbose()) {
    cout << "INFO specex::PSF_Fitter::FitSeveralSpots starting fitting ";
    if(fit_flux) cout << "+flux";
    if(fit_psf) cout << "+psf";
    if(fit_trace) cout << "+trace";
    if(fit_position) cout << "+position";
    if(fit_psf) cout << " npar_psf=" << npar_psf;
    if(fit_trace) cout << " npar_trace=" << npar_trace;
    cout << endl;
  }
  
  
  // allocation
  B.resize(nparTot);
  A.resize(nparTot, nparTot);
  Params.resize(nparTot);
  
  // setting parameters
  {
    int index=0;
    if(fit_psf) {
      for(size_t p=0;p<psf.Params.size();p++) {
	const harp::vector_double& coeff=psf.Params[p].coeff;
	for (size_t k=0; k < coeff.size(); ++k, ++index) 
	  Params(index) = coeff(k);
      }
    }
    
    if(fit_trace) {
      
      for(std::map<int,specex::Trace>::const_iterator it=psf.FiberTraces.begin(); it!=psf.FiberTraces.end(); ++it) {
	const harp::vector_double& PX = it->second.X_vs_lW.coeff;
	for(size_t k=0;k<PX.size();k++,index++)
	  Params(index)=PX(k);
	const harp::vector_double& PY = it->second.Y_vs_lW.coeff;
	for(size_t k=0;k<PY.size();k++,index++)
	  Params(index)=PY(k);
      }
    }

    for(size_t s=0;s<spots.size();s++) {
      const specex::Spot& spot= *(spots[s]);
      if(fit_flux) Params(index++) = spot.flux;
      if(fit_position) {
	Params(index++) = spot.xc;
	Params(index++) = spot.yc;
      }
    }
  }
  
  
  *niter=0;
  int maxiter = 100;
  
  double oldChi2=1e30;
  if(*psfChi2<=0) *psfChi2 = 1e20;
  
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
	specex::Spot &spot= *(spots[s]);
	if(spot.flux<=0) continue;
	harp::vector_double psfParams = psf.FixedCoordParams(spot.xc,spot.yc);
	
	Stamp spot_stamp(image);
	SetStampLimitsFromPSF(spot_stamp,psf,spot.xc,spot.yc);
	for(int j=spot_stamp.begin_j;j<spot_stamp.end_j;j++) {
	  for(int i=spot_stamp.begin_i;i<spot_stamp.end_i;i++) {
	    if(weight(i,j)==0) continue; // no need to compute anything
	    footprint_weight(i,j) += spot.flux*psf.PSFValueWithParams(spot.xc,spot.yc,i,j,psfParams,0,0,0);
	  }
	}
      }
      cout << "for debugging, write image" << endl;
      specex::write_new_fits_image("debugging-model-for-weights.fits",footprint_weight);
      
      for(int j=0;j<footprint_weight.Ny();j++) {
	for(int i=0;i<footprint_weight.Nx();i++) {
	  double model_flux=footprint_weight(i,j);
	  if(model_flux!=0) { // has been computed, pixel participates to footprint
	    double var = readout_noise*readout_noise;
	    if(model_flux>0) { // else negative fluctuation
	      var += model_flux/gain;
	      var += square(flatfield_error*model_flux);
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
	specex::Spot &spot= *(spots[s]);
	Stamp spot_stamp(image);
	SetStampLimitsFromPSF(spot_stamp,psf,spot.xc,spot.yc);
	for(int j=spot_stamp.begin_j;j<spot_stamp.end_j;j++)
	  for(int i=spot_stamp.begin_i;i<spot_stamp.end_i;i++) {
	    footprint_weight(i,j)=weight(i,j);
	    (*npix)++;
	  }
      }
    }
    // mask out some regions with mis-understood lines
    mask.ApplyMaskToImage(footprint_weight,psf,0);
    
    
    if(include_signal_in_weight) 
    {
      cout << "for debugging, write image" << endl;
      write_new_fits_image("debugging-weights.fits",footprint_weight);
      //exit(12); // debug
    }
    
  } // end of test on spot size

  // stuff that was in computechi2ab
  nspots = spots.size();
  npar_fixed_coord = psf.FixedCoordNPar();
  npar_varying_coord = psf.VaryingCoordNPar();
  if(fit_psf) gradPar.resize(npar_fixed_coord); // will remain zero if (!fit_psf)
  if(fit_position || fit_trace) gradPos.resize(2);
  
  spots_flux.resize(nspots);
  spots_x.resize(nspots);
  spots_y.resize(nspots);
  {
    for(int s=0;s<nspots;s++) {
      const specex::Spot& spot=*spots[s];
      spots_flux(s) = spot.flux;
      spots_x(s) = spot.xc;
      spots_y(s) = spot.yc;      
    }
  }
  
  {
    npar_traces = 0;
    if(fit_trace) {
      int index = 0;
      if(fit_psf) index += npar_varying_coord;
      
      // compute first index of trace for each fiber 
      
      XvsW_index_of_fiber.clear();
      YvsW_index_of_fiber.clear();
      int n=0;
      for(std::map<int,specex::Trace>::const_iterator it=psf.FiberTraces.begin(); it!=psf.FiberTraces.end(); ++it) {
	XvsW_index_of_fiber[it->first]=index;
	n = it->second.X_vs_lW.coeff.size();
	index += n;
	npar_traces += n;
	
	YvsW_index_of_fiber[it->first]=index;
	n = it->second.Y_vs_lW.coeff.size();
	index += n;
	npar_traces += n;
      }
      
      // compute monomials of traces for spots
      TraceXvsW_Monomials_of_spots.clear();
      TraceYvsW_Monomials_of_spots.clear();
      for(int s=0;s<nspots;s++) {
	const specex::Spot& spot=*spots[s];   
	specex::Trace &trace = psf.FiberTraces[spot.fiber];
	TraceXvsW_Monomials_of_spots.push_back(trace.X_vs_lW.Monomials(spot.log10_wavelength));
	TraceYvsW_Monomials_of_spots.push_back(trace.Y_vs_lW.Monomials(spot.log10_wavelength));
      }
    }
  }
  
  // compute stamps (don't change them while minimizing)
  spot_stamps.clear();
  for(int s=0;s<nspots;s++) {
    Stamp spot_stamp(image);
    SetStampLimitsFromPSF(spot_stamp,psf,spots_x(s),spots_y(s));
    spot_stamp = spot_stamp.Intersection(stamp); 
    spot_stamps.push_back(spot_stamp);
  }
  


  if(verbose)
    SPECEX_INFO("specex::PSF_Fitter::FitSeveralSpots inc. signal in w=" << include_signal_in_weight << " npix footprint = " << *npix);
  
  
  *niter = 0;
    
  while(true) { // minimization loop 
      
    oldChi2 = *psfChi2;
    if(verbose) SPECEX_INFO("specex::PSF_Fitter::FitSeveralSpots iter=" << *niter << " old chi2=" << oldChi2);
    
    if(verbose) SPECEX_INFO("specex::PSF_Fitter::FitSeveralSpots filling matrix ...");
    *psfChi2 = ComputeChi2AB(spots,true);
    
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
    if( (fit_flux) && (!fit_position) && (!fit_trace) && (!fit_psf) ) {
      if( verbose) SPECEX_INFO("specex::PSF_Fitter::FitSeveralSpots linear because only fit flux");
      linear = true;
    }
    if( (fit_psf) && (!fit_position) && (!fit_trace) && (!fit_flux) && psf.IsLinear()) {
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
	double chi2_0 = ComputeChi2AB(spots,false);
	Params += B;
	double chi2_1 = ComputeChi2AB(spots,false);
	Params -= B;
	if(chi2_1>chi2_0) {
	  use_brent = true;
	  if( verbose) SPECEX_INFO("specex::PSF_Fitter::FitSeveralSpots use brent because step1 increases chi2");
	}else{
	   use_brent = false;
	  if( verbose) SPECEX_INFO("specex::PSF_Fitter::FitSeveralSpots don't brent because chi2 decreases");	  
	}
      }
    }
    if(use_brent) {
      if(verbose) SPECEX_INFO("specex::PSF_Fitter::FitSeveralSpots starting brent ...");
     
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
			  chi2_precision/10.,&bbox,best_chi2,status,100);
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
	if(dchi2>0) {
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
      *psfChi2 = ComputeChi2AB(spots,false);
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
    
    for(size_t p=0;p<psf.Params.size();p++) {
      harp::vector_double& coeff=psf.Params[p].coeff;
      for(size_t c=0;c<coeff.size();c++,index++)
	coeff(c)=Params(index);
    }
    
    
  }
  if (fit_trace) {
    for(std::map<int,specex::Trace>::iterator it=psf.FiberTraces.begin(); it!=psf.FiberTraces.end(); ++it) {
      harp::vector_double& PX = it->second.X_vs_lW.coeff;
      for(size_t k=0;k<PX.size();k++,index++)
	PX(k)=Params(index);
      harp::vector_double& PY = it->second.Y_vs_lW.coeff;
      for(size_t k=0;k<PY.size();k++,index++)
	PY(k)=Params(index);
    }
  }
  
  // save fluxes
  if(verbose) SPECEX_INFO("specex::PSF_Fitter::FitSeveralSpots saving spots fluxes");
  for(size_t s=0;s<spots.size();s++) {
    specex::Spot &spot= *(spots[s]);
    if(fit_flux) {
      spot.flux = Params(index); 
      spot.eflux = sqrt(fitCovmat(index,index));
      index++;
    }
    if(fit_position) {
      int fluxindex = index-1;
      spot.xc = Params(index++);
      spot.yc = Params(index++);
      if(fit_position && fit_flux) {
	if (spot.fxy_CovMat.size1() != 3) spot.fxy_CovMat.resize(3,3);
	for (unsigned k=0; k<3; ++k) 
	  for (unsigned l=0; l<3 ;++l)
	    spot.fxy_CovMat(k,l) = fitCovmat(fluxindex+k, fluxindex+l);
      }
    }
  }
  
  // also save psf in spots;
  if (fit_psf) {
    if(verbose) SPECEX_INFO("specex::PSF_Fitter::FitSeveralSpots saving spots psf params");
    for(size_t s=0;s<spots.size();s++) {
      specex::Spot &spot= *(spots[s]);
      spot.PSFParams = psf.FixedCoordParams(spot.xc,spot.yc);
    } 
  }
  // also save coordinates of spots
  if (fit_trace) {
    if(verbose) SPECEX_INFO("specex::PSF_Fitter::FitSeveralSpots saving spots coord from fitted traces");
    for(size_t s=0;s<spots.size();s++) {
      specex::Spot &spot= *(spots[s]);
      spot.xc = psf.FiberTraces[spot.fiber].X_vs_lW.Value(spot.log10_wavelength);
      spot.yc = psf.FiberTraces[spot.fiber].Y_vs_lW.Value(spot.log10_wavelength);
    } 
  }
   
  //
  // spot.chi2=*psfChi2;
  
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
 
bool specex::PSF_Fitter::FitOneSpot(specex::Spot& spot, double *chi2, int *n_iterations) {

  specex::Spot saved_spot = spot;

  vector<specex::Spot*> spots; 
  spots.push_back(&spot);
  
  int npar_psf = psf.FixedCoordNPar();
  psf.Params.clear();
  for(int p=0;p<npar_psf;p++) {
    psf.Params.push_back(specex::Legendre2DPol(0,0,image.Nx(),0,0,image.Ny()));
    psf.Params[p].coeff(0) = spot.PSFParams(p);
  }
  
  
  bool ok = FitSeveralSpots(spots,&(spot.chi2),n_iterations);
  
  if(verbose) {
    
    if(ok) {
      cout << "INFO specex::PSF_Fitter::FitOneSpot successful fit of flux";
      if(fit_position) cout << "+pos";
      if(fit_psf) cout << "+psf";
      cout << " chi2= " << spot.chi2;
      cout << " dchi2= " <<  saved_spot.chi2-spot.chi2;
      if(saved_spot.flux>0)
	 cout << " dflux/flux=" << spot.flux/saved_spot.flux-1;
      cout << " dx=" << spot.xc-saved_spot.xc;
      cout << " dy=" << spot.yc-saved_spot.yc;
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
  psf.Params.clear();
  for (size_t k =0; k < ParamsToSet.size(); ++k) {
    psf.Params.push_back(specex::Legendre2DPol(0,0,image.Nx(),0,0,image.Ny()));
    psf.Params.back().coeff(0)=ParamsToSet(k);
  }
}

bool specex::PSF_Fitter::InterpolateSpotPSFs(vector<specex::Spot*>& spots, double *chi2_val, int *n_iterations) {
  
  double chi2_memory_slot;
  double *chi2 = &chi2_memory_slot;
  if(chi2_val) chi2 = chi2_val;
  
  int niter_memory_slot;
  int *niter = &niter_memory_slot;
  if(n_iterations) niter = n_iterations;
  
  int number_of_wavelength =  find_spot_arrays(spots).size();
  int polynomial_degree_along_x_for_psf_tail = 0;
  int polynomial_degree_along_y_for_psf_tail = min(polynomial_degree_along_y,2);
  


  // allocate psf polynomials
  int npar = psf.FixedCoordNPar();

  SPECEX_INFO("Number of PSF params at fixed coord = " << npar);

  psf.Params.clear();

  int first_index_for_tail = npar;
  if(psf.analyticPSF->HasParam("YTailNorm")) first_index_for_tail = min(first_index_for_tail,psf.analyticPSF->ParamIndex("YTailNorm"));
  if(psf.analyticPSF->HasParam("TailNorm")) first_index_for_tail = min(first_index_for_tail,psf.analyticPSF->ParamIndex("TailNorm"));
  
  for(int p=0;p<npar;p++) {
    if(p<first_index_for_tail)
      psf.Params.push_back(specex::Legendre2DPol(polynomial_degree_along_x,0,image.Nx(),polynomial_degree_along_y,0,image.Ny()));
    else
      psf.Params.push_back(specex::Legendre2DPol(polynomial_degree_along_x_for_psf_tail,0,image.Nx(),polynomial_degree_along_y_for_psf_tail,0,image.Ny()));
  }

  
  
  int nparpsf = psf.Params.size();
  int npartot = 0;
  for(int p=0;p<nparpsf;p++) {
    const specex::Legendre2DPol &pol = psf.Params[p];
    int nparpol = (pol.xdeg+1)*(pol.ydeg+1);
    npartot += nparpol;
  }

  if( int(spots.size())<npartot/nparpsf ) {
    cout << "WARNING specex::PSF_Fitter::InterpolateSpotPSFs cannot fit because not enough spots" << endl;
    return false;
  }

  harp::matrix_double A(npartot,npartot); A*=0;
  harp::vector_double B(npartot); B*=0;
  
  // dpsfdcoef = dpsf/dpar * dpar/dcoef

  //harp::matrix_double W(nparpsf,nparpsf); W*= 0;
  //for(int i=0;i<nparpsf;i++)
  //W(i,i)=1;
  
  harp::matrix_double H(nparpsf,npartot);
  
  for(size_t s=0;s<spots.size();s++) {
    const specex::Spot& spot = *(spots[s]);
    
    // filling the matrix of derivatives
    H *= 0;
    int index=0;
    for(int p=0;p<nparpsf;p++) {
      const specex::Legendre2DPol &pol = psf.Params[p];
      harp::vector_double m=pol.Monomials(spot.xc,spot.yc);
      for(size_t k=0;k<m.size();k++,index++)
	H(p,index)=m[k];
    }
    
#ifdef USE_PARAMS_UNCERTAINTIES
    { 
      const harp::matrix_double &W=spot.PSFParams_WeightMat;
      if(W.size1()!=nparpsf)  SPECEX_ERROR("PSFParams weight mat not set");
      
      harp::matrix_double WH(W.size1(),H.size2());
      
      //specex::symm('L',1,W,H,1,WH); // cannot compile, WH = W * H
      specex::gemm(1,W,H,0,WH); //  WH = W * H
      specex::gemm(1,ublas::trans(H),WH,1,A);  // A += Ht*spot.PSFParams_WeightMat*H;
      specex::gemv(1,ublas::trans(WH),spot.PSFParams,1,B); // B += Ht*W*Mat(spot.PSFParams);
    }
#else
    specex::syrk(1,ublas::trans(H),1,A); // faster than  specex::gemm(1,ublas::trans(H),H,1,A); // = A += HHt
    specex::gemv(1,ublas::trans(H),spot.PSFParams,1,B);    
#endif
    
    //if(int(s)%100==0 && s>0) cout << "processed " << s << " spots" << endl;
  }

  
  
  // actually solve
  //harp::matrix_double As=A;
  int status = cholesky_solve(A,B);
  if(status) {
    //write_new_fits_image("A.fits",As);
    SPECEX_ERROR("cholesky_solve failed in specex::PSF_Fitter::InterpolateSpotPSFs with status = " << status);
  }
  // push result into PSF
  {
    int index=0;
    for(int p=0;p<nparpsf;p++) {
      specex::Legendre2DPol &pol = psf.Params[p];
      for(size_t k=0;k<pol.coeff.size();k++,index++)
	pol.coeff(k)=B(index);
    }
  }
  
  // add global psf parameters to spots
  for(size_t s=0;s<spots.size();s++) {
    specex::Spot& spot = *(spots[s]);
    spot.GlobalPSFParams = psf.FixedCoordParams(spot.xc,spot.yc);
  }
  
  cout << "INFO specex::PSF_Fitter::InterpolateSpotPSFs successful" << endl;
  
  return true;
}

bool specex::PSF_Fitter::FitTraces(vector<specex::Spot*>& spots, int *n_fibers_fitted) {
  
  SPECEX_INFO("specex::PSF_Fitter::FitTraces starting");

  int nok_memory_slot;
  int *nok = &nok_memory_slot;
  if(n_fibers_fitted) nok = n_fibers_fitted;
  
  *nok = 0;
  for(map<int,specex::Trace>::iterator it=psf.FiberTraces.begin();
      it !=psf.FiberTraces.end(); ++it) {
    
    specex::Trace& trace=it->second;
    
    vector<specex::Spot*> spots_of_fiber;
    for(size_t s=0;s<spots.size();s++) {
      specex::Spot* spot = spots[s];
      if(spot->fiber==trace.fiber)
	spots_of_fiber.push_back(spot);
    }
    bool ok = trace.Fit(spots_of_fiber);
    if(ok) { 
      (*nok)++;
      
      for(size_t s=0;s<spots_of_fiber.size();s++) {
	specex::Spot& spot = *(spots_of_fiber[s]);
	spot.xc=trace.X_vs_lW.Value(spot.log10_wavelength);
	spot.yc=trace.Y_vs_lW.Value(spot.log10_wavelength);

      }


    }
  }
  SPECEX_INFO("specex::PSF_Fitter::FitTraces ended nok = " << (*nok) << "/" << psf.FiberTraces.size());
 
  return true;
}

bool specex::PSF_Fitter::FitIndividualSpotFluxes(std::vector<specex::Spot*>& spots) {
  
  SPECEX_INFO("fitting independently the flux of each spot");
  include_signal_in_weight = false;
  fit_flux                 = true;
  fit_position             = false;
  fit_psf                  = false;
  fit_trace                = false;
  verbose                  = false; // debug
  fatal                    = false; // debug
  psf.hSizeX = 4; // 
  psf.hSizeY = 4; //
  
  int nok=0;
  for(size_t s=0;s<spots.size();s++) {
    specex::Spot& spot = *(spots[s]);    
    spot.eflux = 0;
    spot.flux = 0;	
    spot.status=-1;
    
    bool ok = FitOneSpot(spot);
    
    if(ok) {
      nok++;
      spot.status=1;
    }else{
      spot.status=0;
      SPECEX_WARNING("fit of flux of spot at x = " << spot.xc << " " << spot.yc << " failed");
    }
    if(int(s)%100==0 && s!=0) SPECEX_INFO("done " << s << "/" << spots.size() << " ...");
  }
  SPECEX_INFO("successful fit of each spot flux for " << nok << "/" << spots.size());
  return true;
}
bool specex::PSF_Fitter::FitIndividualSpotPositions(std::vector<specex::Spot*>& spots) {
  
  SPECEX_INFO("fitting independently the flux+position of each spot");
  include_signal_in_weight = false;
  fit_flux                 = true;
  fit_position             = true;
  fit_psf                  = false;
  fit_trace                = false;
  verbose                  = false;
  fatal                    = false;
  psf.hSizeX = 4; // 
  psf.hSizeY = 4; //
  int nok=0;
  for(size_t s=0;s<spots.size();s++) {
    specex::Spot& spot = *(spots[s]);    
    spot.initial_xc = spot.xc;
    spot.initial_yc = spot.yc;
    spot.eflux = 0;
    spot.flux = 0;	
    spot.status=-1;
    
    bool ok = FitOneSpot(spot);
    
    if(ok) {
      nok++;
      spot.status=1;
    }else{
      spot.status=0;
    }
    if(int(s)%100==0 && s!=0) SPECEX_INFO("done " << s << "/" << spots.size() << " ...");
  }
  SPECEX_INFO("successful fit of each spot flux+pos for " << nok << "/" << spots.size());
  return true;
}

bool specex::PSF_Fitter::FitEverything(std::vector<specex::Spot*>& spots, bool init_psf) {

  if(spots.size()==0) {
    SPECEX_ERROR("in fit_several_spots spotarrays is empty");
  }

  if(init_psf) {

    harp::vector_double startParams; 
    // set default parameters
    psf.analyticPSF->InitParams(1.0,startParams);
    int nparpsf = psf.analyticPSF->NPar();
    for(size_t s=0;s<spots.size();s++) {
      specex::Spot& spot = *(spots[s]);
      spot.PSFParams = startParams;
      spot.PSFname = psf.analyticPSF->Name();
      spot.PSFParams_WeightMat.resize(nparpsf,nparpsf); spot.PSFParams_WeightMat*=0;
      spot.PSFParams_CovMat.resize(nparpsf,nparpsf); spot.PSFParams_CovMat *= 0;
      for(int i=0;i<nparpsf;i++) {
	spot.PSFParams_WeightMat(i,i)=1;
	spot.PSFParams_CovMat(i,i)=1;
      }
    }
  }
  
  
  FitIndividualSpotFluxes(spots);
  
  verbose = true;
  
  double chi2=1e30;
  int npix = 0;
  int niter = 0;
    
  SPECEX_INFO("starting to fit PSF with " <<  spots.size() << " spots");
 
  if(init_psf) {
    
    SPECEX_INFO("init PSF");
    
    FitTraces(spots);

    psf.hSizeX = 4; //
    psf.hSizeY = 4; // 
    
    // here we need to count number of wavelength and spots per wavelength
    {
      std::vector<SpotArray> spot_arrays = find_spot_arrays(spots);
      int nspots_per_wave=0;
      for(size_t a=0;a<spot_arrays.size();a++) {
	int asize = spot_arrays[a].size();
	if(asize>nspots_per_wave) 
	  nspots_per_wave=asize;
      }
      polynomial_degree_along_x = min(nspots_per_wave-1,1);
      polynomial_degree_along_y = min(int(spot_arrays.size())-1,3);
      
      SPECEX_INFO("Setting PSF polynomial degrees " << polynomial_degree_along_x << " " << polynomial_degree_along_y);
    }

    
    bool ok = InterpolateSpotPSFs(spots,&chi2,&niter);
    if(ok) {
      SPECEX_INFO("InterpolateSpotPSFs successful");
      for(size_t s=0;s<spots.size();s++)
	spots[s]->status=2;
    }
    else {
      SPECEX_ERROR("InterpolateSpotPSFs failed");
    }
    
  }
  
  fatal = true;
  include_signal_in_weight = false;
  chi2_precision = 10;
  bool ok = true;
  
  if(psf.analyticPSF->Name()=="GAUSSHERMITE") {
    psf.hSizeX = 4; // no fit of tails to start
    psf.hSizeY = 4; // no fit of tails to start
    psf.ClearPriors();
    // from a not yet perfect fit of laser data in r1 CCD at lambda=8702 , expnum=00120726
    cout << "INFO Setting priors to PSF tails" << endl;
    cout << "DEBUG tail indices :" << endl;
    cout << "DEBUG TailNorm : " << psf.analyticPSF->ParamIndex("TailNorm") << endl;
    cout << "DEBUG TailPower : " << psf.analyticPSF->ParamIndex("TailPower") << endl;
    cout << "DEBUG TailXposScale : " << psf.analyticPSF->ParamIndex("TailXposScale") << endl;
    cout << "DEBUG TailXnegScale : " << psf.analyticPSF->ParamIndex("TailXnegScale") << endl;
    cout << "DEBUG TailYnegScale : " << psf.analyticPSF->ParamIndex("TailYnegScale") << endl;
    cout << "DEBUG Npar : " << psf.analyticPSF->NPar() << endl;
    
    double epsilon = 1.e-12;
    if(psf.analyticPSF->HasParam("YTailNorm")) psf.SetPrior("YTailNorm",new GaussianPrior(0,epsilon)); // 
    if(psf.analyticPSF->HasParam("TailNorm")) psf.SetPrior("TailNorm",new GaussianPrior(0,epsilon)); // 
    if(psf.analyticPSF->HasParam("TailPower")) psf.SetPrior("TailPower",new GaussianPrior(2.,epsilon)); // 
    if(psf.analyticPSF->HasParam("TailXposScale")) psf.SetPrior("TailXposScale",new GaussianPrior(1,epsilon));
    if(psf.analyticPSF->HasParam("TailXnegScale")) psf.SetPrior("TailXnegScale",new GaussianPrior(1,epsilon));
    if(psf.analyticPSF->HasParam("TailYnegScale")) psf.SetPrior("TailYnegScale",new GaussianPrior(1,epsilon));
    
    //exit(12);
  }

  
  
  SPECEX_INFO("Starting FitSeveralSpots FLUX");
  SPECEX_INFO("=======================================");
  fit_flux       = true;
  fit_position   = false;
  fit_psf        = false;
  fit_trace      = false;
  ok = FitSeveralSpots(spots,&chi2,&npix,&niter);
  if(!ok) SPECEX_ERROR("FitSeveralSpots failed for FLUX");
  
  // write_spots_list(spots,"spots-tmp-FLUX.list",PSF);
  
  SPECEX_INFO("Starting FitSeveralSpots TRACE");
  SPECEX_INFO("=======================================");
  fit_flux       = false;
  fit_position   = false;
  fit_psf        = false;
  fit_trace      = true;
  ok = FitSeveralSpots(spots,&chi2,&npix,&niter);
  if(!ok) SPECEX_ERROR("FitSeveralSpots failed for TRACE");
  
  SPECEX_INFO("Starting FitSeveralSpots TRACE+FLUX");
  SPECEX_INFO("=======================================");
  fit_flux       = true;
  fit_position   = false;
  fit_psf        = false;
  fit_trace      = true;
  ok = FitSeveralSpots(spots,&chi2,&npix,&niter);
  if(!ok) SPECEX_ERROR("FitSeveralSpots failed for TRACE+FLUX");
  
  // write_spots_list(spots,"spots-tmp-TRACE+FLUX.list",PSF);

  
  SPECEX_INFO("Starting FitSeveralSpots PSF");
  SPECEX_INFO("=======================================");
  fit_flux       = false;
  fit_position   = false;
  fit_psf        = true;
  fit_trace      = false;
  ok = FitSeveralSpots(spots,&chi2,&npix,&niter);
  if(!ok) SPECEX_ERROR("FitSeveralSpots failed for PSF");
  
  // write_spots_list(spots,"spots-tmp-psf.list",PSF);
  
  SPECEX_INFO("Starting FitSeveralSpots PSF+FLUX");
  SPECEX_INFO("=======================================");
  fit_flux       = true;
  fit_position   = false;
  fit_psf        = true;
  fit_trace      = false;
  ok = FitSeveralSpots(spots,&chi2,&npix,&niter);
  if(!ok) SPECEX_ERROR("FitSeveralSpots failed for PSF+FLUX");

  // write_spots_list(spots,"spots-tmp-PSF+FLUX.list",PSF);
  // psf.write("psf-tmp-PSF+FLUX.dat");

  return ok;

#ifdef LA_SUITE

  cout << "==== Starting FitSeveralSpots TRACE+FLUX+PSF ==== " << endl;
  fit_flux       = true;
  fit_position   = false;
  fit_psf        = true;
  fit_trace      = true;
  chi2_precision = 10;
  ok = FitSeveralSpots(spots,&chi2,&npix,&niter);
  if(!ok) {
    cout << "FitSeveralSpots failed for TRACE" << endl;
    return false;
  }
  
  write_spots_list(spots,"spots-TRACE+FLUX+psf.list",PSF);
  write_spots_data(spots,"spots-TRACE+FLUX+psf.dat");
  psf.write("psf-TRACE+FLUX+psf.dat");
  
  if(fit_psf_tails) {
    if(psf.analyticPSF->Name()=="GAUSSHERMITE") {
      cout << "==== Starting FitSeveralSpots PSF TAILS ==== " << endl;
      psf.hSizeX = 50; //  fit  tails !
      psf.hSizeY = 500; //  fit  tails !
      psf.ClearPriors();
      
      
      // from a not yet perfect fit of laser data in r1 CCD at lambda=8702 , expnum=00120726
      double epsilon=1.e-12;
      //if(psf.analyticPSF->HasParam("TailNorm")) psf.SetPrior("TailNorm",new GaussianPrior(1.8,epsilon)); // slope
      if(psf.analyticPSF->HasParam("TailPower")) psf.SetPrior("TailPower",new GaussianPrior(2.,epsilon)); // slope
      if(psf.analyticPSF->HasParam("TailXposScale")) psf.SetPrior("TailXposScale",new GaussianPrior(1.2,epsilon));
      if(psf.analyticPSF->HasParam("TailXnegScale")) psf.SetPrior("TailXnegScale",new GaussianPrior(1.2,epsilon));
      if(psf.analyticPSF->HasParam("TailYnegScale")) psf.SetPrior("TailYnegScale",new GaussianPrior(0.9,epsilon));
      
      include_signal_in_weight = false; // to deweight the core for the determination of tails
      fit_flux       = false;
      fit_position   = false;
      fit_psf        = true;
      fit_trace      = false;
      chi2_precision = 10;
      ok = FitSeveralSpots(spots,&chi2,&npix,&niter);
      if(!ok) {
	cout << "FitSeveralSpots failed for PSF+FLUX" << endl;
	return false;
      }
      write_spots_data(spots,"spots-PSF+TAILS.dat");
      write_spots_list(spots,"spots-PSF+TAILS.list",PSF);
      psf.write("psf-PSF+TAILS.dat");
      
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
      ok = FitSeveralSpots(spots,&chi2,&npix,&niter);
      if(!ok) {
	cout << "FitSeveralSpots failed for PSF+TAILS+FLUX+TRACE" << endl;
	return false;
      }
    }
    write_spots_data(spots,"spots-PSF+TAILS+FLUX+TRACE.dat");
    write_spots_list(spots,"spots-PSF+TAILS+FLUX+TRACE.list",PSF);
    psf.write("psf-PSF+TAILS+FLUX+TRACE.dat");
  }// end of test fit_psf_tails

#endif
  
  return ok;
}
