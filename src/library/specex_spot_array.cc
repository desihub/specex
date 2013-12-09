#include <vector>
#include <string>
#include <cmath>

#include "specex_spot_array.h"
#include "specex_message.h"

using namespace std;

vector<specex::Spot*> specex::array_of_pointer(vector<specex::Spot> &spots) {
  
  vector<specex::Spot*> spot_pointers;
  for(size_t s=0;s<spots.size();s++) {
    specex::Spot& spot = spots[s];
    spot_pointers.push_back(&spot);
  }
  return spot_pointers;
}


vector<specex::SpotArray> specex::find_spot_arrays( vector<specex::Spot*> &spots, bool check_status) {
  
  vector<specex::SpotArray> spotarrays;
    
  {
    for(size_t s=0;s<spots.size();s++) {
      specex::Spot& spot = *(spots[s]);
      if(check_status && spot.status==0) continue;
      specex::SpotArray* spotarray=0;
      for(size_t a=0;a<spotarrays.size();a++) {
	specex::SpotArray& testarray=spotarrays[a];
	if(spot.fiber_bundle == testarray.fiber_bundle &&
	   fabs(pow(10,spot.log10_wavelength) - pow(10,testarray.log10_wavelength))<1) {
	  spotarray = &testarray;
	  break;
	}
      }
      if(!spotarray) {
	specex::SpotArray newarray;
	newarray.fiber_bundle = spot.fiber_bundle;
	newarray.log10_wavelength = spot.log10_wavelength;
	spotarrays.push_back(newarray);
	spotarray=&(spotarrays[spotarrays.size()-1]);
      }
      spotarray->push_back(&spot);
    }
  }
  return spotarrays;
}

vector<specex::SpotArray> specex::find_spot_arrays( vector<specex::Spot> &spots, bool check_status) {
  vector<specex::Spot*> spot_pointers = specex::array_of_pointer(spots);
  return find_spot_arrays(spot_pointers,check_status);
}

vector<specex::SpotArray> specex::isolated_spot_arrays(  vector<specex::SpotArray>& other, const double& delta_lambda) {
  vector<specex::SpotArray> iarray;
  for(size_t a1=0;a1<other.size();a1++) {
    const double& w1=pow(10,other[a1].log10_wavelength);
    bool ok=true;
    for(size_t a2=0;a2<other.size();a2++) {
      if(a1!=a2) {
	const double& w2=pow(10,other[a2].log10_wavelength);
	if(fabs(w1-w2)<delta_lambda) {
	  ok = false;
	  break;
	}
      }
    }
    if(ok) iarray.push_back(other[a1]);
  }
  return iarray;
}

vector<specex::SpotArray> specex::valid_spot_arrays(  vector<specex::SpotArray>& other, int required_status, const double& min_valid_fraction) {
  vector<specex::SpotArray> varrays;
  for(size_t a=0;a<other.size();a++) {
     specex::SpotArray& other_array = other[a];
    specex::SpotArray array;
    int ntot=other.size();
    
    for(size_t s=0;s<other_array.size();s++) {
      specex::Spot& spot=*(other_array[s]);
      if(spot.status>=required_status)
	array.push_back(&spot);
    }
    float frac=float(array.size()/other_array.size());
    if(frac>=min_valid_fraction)
      varrays.push_back(array);
  }
  return varrays;
}

specex::SpotArray  specex::merge_spot_arrays(  vector<specex::SpotArray>& other) {
  specex::SpotArray array;
  for(size_t a=0;a<other.size();a++) {
    specex::SpotArray& other_array = other[a];
    for(size_t s=0;s<other_array.size();s++) {
      array.push_back(other_array[s]);
    }
  }
  return array;
}

void specex::write_spots_list(vector<Spot*> spots, const string& filename) {
  ofstream osl(filename.c_str());  
  for(size_t s=0;s<spots.size();s++) {
    Spot& spot = *(spots[s]);
    
    if(spot.GlobalPSFParams.size()==0) {
      spot.GlobalPSFParams.resize(spot.PSFParams.size());
      spot.GlobalPSFParams *= 0;
    }
    if(s==0) spot.write_list_header(osl);
    spot.write_list_entry(osl);
  }
  osl.close();
  
  SPECEX_INFO("wrote spots in " << filename);
  
}

void specex::write_spots_list(std::vector<Spot> spots, const std::string& filename) {
  write_spots_list(array_of_pointer(spots),filename);
}
