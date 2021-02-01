#include <string>
#include <iostream>
#include <specex_pyimage.h>

using namespace std;

specex::PyImage::PyImage(
			 py::array ds_pix,
			 py::array ds_ivar,
			 py::array ds_mask,
			 py::array ds_readnoise,
			 std::map<std::string,std::string> hdr_mt
			 ){
  
  auto pix_prop       = ds_pix.request();
  auto ivar_prop      = ds_ivar.request();
  auto mask_prop      = ds_mask.request();
  auto readnoise_prop = ds_readnoise.request();

  double *pix_vals       = (double*) pix_prop.ptr;
  double *ivar_vals      = (double*) ivar_prop.ptr;
  int    *mask_vals      = (   int*) mask_prop.ptr;
  double *readnoise_vals = (double*) readnoise_prop.ptr;

  // assume all images equal in size to ds_pix
  size_t nx = pix_prop.shape[0];
  size_t ny = pix_prop.shape[1];

  // this is to match the existing ordering as used in specex with
  // harp::image i.e. specex::image_data objects
  image.resize(ny,nx);
  weight.resize(ny,nx);
  mask.resize(ny,nx);
  rdnoise.resize(ny,nx);
  
  for (size_t i = 0; i < nx*ny; i++){
    image.data[i]   = pix_vals[i];
    weight.data[i]  = ivar_vals[i];
    mask.data[i]    = mask_vals[i];
    rdnoise.data[i] = readnoise_vals[i];
  }

  header = hdr_mt;
  
}

std::vector<double> specex::PyImage::get_data(std::string tag){

  auto v = this->image.data;
  if(tag == "weight") { v = this->weight.data;}
  if(tag == "mask")   { v = this->mask.data;}
  if(tag == "rdnoise"){ v = this->rdnoise.data;}

  if(tag == "image")
  cout << tag << " array sizes " << image.n_cols() << " " <<
    image.n_rows() << " " << image.data.size() << endl;

  std::vector<double> vi(v.size());
  for (unsigned i = 0; i < v.size(); ++i) vi[i]=v[i];
  return vi;

}

std::map<std::string,std::string> specex::PyImage::get_header(){
  return this->header;
}
