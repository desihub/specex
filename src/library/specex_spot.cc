

#include <assert.h>
#include <iomanip>

#include "specex_linalg.h"
#include "specex_spot.h"


using namespace std;


void specex::Spot::write_list_header(ostream& os)  const {

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
}

