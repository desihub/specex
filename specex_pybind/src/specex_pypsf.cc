#include <string>
#include <iostream>
#include <specex_pypsf.h>
using namespace std;

specex::image_data specex::PyPSF::get_trace(std::string axis){

  ///specex::PSF_p psf = this->psf;
  
  if(axis == "x"){
    return this->psf->pydata.coeff2d_x;
  }
  else{
    return this->psf->pydata.coeff2d_y;
  }
}

void specex::PyPSF::get_table(std::vector<std::string>         &table_string,
			      std::vector<std::vector<int>>    &table_int,
			      std::vector<std::vector<double>> &table_double){

  int status = 0;
  int nrows  = this->psf->pydata.table.data.size();
  int ncols  = this->psf->pydata.table.columns.size();

  cout << "nrows in get_tablestring = " << nrows << endl;
  cout << "ncols in get_tablestring = " << ncols << endl;
  
}

void specex::PyPSF::SetParamsOfBundle(){

  for(std::map<int,specex::PSF_Params>::const_iterator bundle_it = this->psf->ParamsOfBundles.begin();
      bundle_it != this->psf->ParamsOfBundles.end(); ++bundle_it) {
      
    const specex::PSF_Params & params_of_bundle = bundle_it->second;
      
    int ndf = params_of_bundle.ndata - params_of_bundle.nparams;
    double chi2pdf = 0;
    
    if(ndf>0) chi2pdf = params_of_bundle.chi2/ndf;

    bundle_id.push_back(params_of_bundle.bundle_id);
    bundle_ndata.push_back(params_of_bundle.ndata);
    bundle_nparams.push_back(params_of_bundle.nparams);
    bundle_chi2pdf.push_back(chi2pdf);
    
  }
}
