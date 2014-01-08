

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
  os << "# wave : in A" << endl;
  os << "# fiber : " << endl;
  os << "# fiber_bundle : " << endl;
  os << "# x : " << endl;
  os << "# y : " << endl;
  os << "# flux : " << endl;
  os << "# eflux: error on flux" << endl;
  os << "# chi2 : " << endl;
  os << "# chi2pdf : " << endl;
  os << "# status : " << endl;
  os << "# ix : initial x coordinate" << endl;
  os << "# iy : initial y coordinate " << endl;
  os << "# iflux : initial flux " << endl;
  os << "#end" << endl;
}

void specex::Spot::write_list_entry(ostream& os)  const {
  os << wavelength << " ";
  os << fiber << " ";
  os << fiber_bundle << " ";
  os << xc << " ";
  os << yc << " ";
  os << flux << " ";
  os << eflux << " ";
  os << chi2 << " ";
  os << 0 << " ";
  os << status << " ";
  os << initial_xc << " ";
  os << initial_yc << " ";
  os << initial_flux << " ";
  os << endl;
  
}

