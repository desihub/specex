

#include <assert.h>
#include <iomanip>

#include "specex_linalg.h"
#include "specex_spot.h"


using namespace std;


void specex::Spot::write_list_header(ostream& os)  const {

  int index_psf_y_tail_norm = -1;
  int index_psf_sym_tail_norm = -1;
  
  /*
  if(PSF) {
    if(PSF->analyticPSF->HasParam("TailNorm")) index_psf_sym_tail_norm = PSF->analyticPSF->ParamIndex("TailNorm");
    if(PSF->analyticPSF->HasParam("YTailNorm")) index_psf_y_tail_norm = PSF->analyticPSF->ParamIndex("YTailNorm");
    
  }
  */
  
  os << setprecision(10);
  os << "# lwave : log10(wavelength) with wavelength in A" << endl;
  os << "# fiber : " << endl;
  os << "# fiber_bundle : " << endl;
  os << "# x : " << endl;
  os << "# y : " << endl;
  os << "# flux : " << endl;
  os << "# ex : error on x" << endl;
  os << "# ey : error on y" << endl;
  os << "# eflux: error on flux" << endl;
  for(size_t i=0;i<PSFParams.size(); i++) {
    os << "# psf" << i << ": ";
    if(i==index_psf_sym_tail_norm) os << "PSF symmetric tail amplitude";
    if(i==index_psf_y_tail_norm) os << "PSF y tail amplitude";
    
    os << endl;
  }
  
  
  os << "# chi2 : " << endl;
  os << "# chi2pdf : " << endl;
  os << "# status : " << endl;
  os << "# ix : initial x coordinate" << endl;
  os << "# iy : initial y coordinate " << endl;
  os << "# iflux : initial flux " << endl;
  os << "#end" << endl;
}

void specex::Spot::write_list_entry(ostream& os)  const {
  os << log10_wavelength << " ";
  os << fiber << " ";
  os << fiber_bundle << " ";
  os << xc << " ";
  os << yc << " ";
  os << flux << " ";
  if(fxy_CovMat.size1()==3 && fxy_CovMat.size2()) {
    os << sqrt(fxy_CovMat(1,1)) << " ";
    os << sqrt(fxy_CovMat(2,2)) << " ";
    os << sqrt(fxy_CovMat(0,0)) << " ";
  }else{
    os << "0 0 0 ";
  }
  for(size_t i=0;i<PSFParams.size(); i++)
    os << PSFParams(i) << " ";
  os << chi2 << " ";
  os << 0 << " ";
  os << status << " ";
  os << initial_xc << " ";
  os << initial_yc << " ";
  os << initial_flux << " ";
  os << endl;
  
}

/*
void specex::Spot::write_full_entry(ostream& os)  const {
  os << "BeginSpot" << endl;
  os << log10_wavelength << " " << fiber << " " << fiber_bundle << endl;
  os << xc << " " << yc << " " << flux << endl;
  os << initial_xc << " " << initial_yc << " " << initial_flux << endl;
  os << PSFname << endl;
  PSFParams.writeASCII(os);
  os << fxy_CovMat << endl;
  os << PSFParams_WeightMat << endl;
  os << chi2 << " " << status << endl;
  os << "EndSpot" << endl;
}


bool specex::Spot::read_full_entry(istream& is) {
  string label;
  if(! (is >> label)) return false;
  if(! (label=="BeginSpot")) return false;
  if(! (is >> log10_wavelength >> fiber >> fiber_bundle)) return false;
  if(! (is >> xc >> yc >> flux)) return false;
  if(! (is >> initial_xc >> initial_yc >> initial_flux)) return false;
  if(! (is >> PSFname)) return false;
  if( PSFParams.readASCII(is)!=0 ) return false;
  if( fxy_CovMat.readASCII(is)!=0 ) return false;
  if( PSFParams_WeightMat.readASCII(is)!=0 ) return false;
  if(! (is >> chi2 >> status)) return false;
  if(! (is >> label)) return false;
  if(! (label=="EndSpot")) return false;
  
  // sanity checks
  if(fxy_CovMat.SizeX()!=0) {
    if(fxy_CovMat.SizeX() != 3 || fxy_CovMat.SizeY() != 3) {
      cout << "ERROR specex::Spot::read_full_entry xyflux matrix wrong size" << endl;
      return false;
    }
  }
  if(fxy_CovMat.SizeX()>0)
    eflux = sqrt(fxy_CovMat(0,0));
  else
    eflux = 0;
  
  if(PSFParams_WeightMat.SizeX()!=0) {
    if(PSFParams.size() != PSFParams_WeightMat.SizeX() || PSFParams.size() != PSFParams_WeightMat.SizeY()) {
      cout << "ERROR specex::Spot::read_full_entry vector and matrix don't match" << endl;
      return false;
    }
    
    // PSFParams_CovMat.allocate(,0); // we probably don't need it now
    PSFParams_CovMat = PSFParams_WeightMat;
    PSFParams_CovMat.CholeskyInvert("L"); // should never fail
    PSFParams_CovMat.Symmetrize("L");
  }
  
  return true;
  }
*/
