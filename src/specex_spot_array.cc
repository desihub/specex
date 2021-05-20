#include <vector>
#include <string>
#include <cmath>
#include <fstream>

#include "specex_spot_array.h"
#include "specex_message.h"

using namespace std;

void specex::write_spots_list(vector<Spot_p> spots, const string& filename) {
  ofstream osl(filename.c_str());  
  for(size_t s=0;s<spots.size();s++) {
    Spot_p& spot = spots[s];
    if(s==0) {spot->write_list_header(osl); osl << "#end" << endl;  }
    spot->write_list_entry(osl);
    osl << endl;
  }
  osl.close();
  
  SPECEX_INFO("wrote spots in " << filename);
  
}

void specex::write_spots_list(vector<Spot_p> spots, const PSF_p psf, const string& filename) {
  ofstream osl(filename.c_str());
  
  std::vector<std::string> pnames = psf->DefaultParamNames();
  
  for(size_t s=0;s<spots.size();s++) {
    Spot_p& spot = spots[s];
    
    if(s==0) {
      spot->write_list_header(osl); 
      for(size_t p=0;p<pnames.size();p++)
	osl << "# " << pnames[p] << " :" << endl;
      osl << "#end" << endl;  
    }
    
    spot->write_list_entry(osl);

    unbls::vector_double params = psf->AllLocalParamsFW(spot->fiber,spot->wavelength,spot->fiber_bundle);
    for(size_t p=0;p<params.size();p++)
      osl << " " << params[p];
    
    osl << endl;
    
  }
  osl.close();
  
  SPECEX_INFO("wrote spots in " << filename);
  
}
