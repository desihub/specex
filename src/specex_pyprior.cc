#include <string>
#include <iostream>
#include <specex_message.h>
#include <specex_pyprior.h>

using namespace std;

int specex::PyPrior::set_priors(specex::PyOptions opts){

  map<string,Prior*> priors;
  {
    int np=int(opts.argurment_priors.size());
    if( ! (np%3==0) ) {
      cerr << "error in parsing, priors must be of the form 'name value error' (np=" << np << ")" << endl;
      cerr << opts.desc << endl;
      return EXIT_FAILURE;
    }
    try {
      for(int i=0;i<np/3;i++) {
	string pname=opts.argurment_priors[3*i];
	double val=atof(opts.argurment_priors[3*i+1].c_str());
	double err=atof(opts.argurment_priors[3*i+2].c_str());
	cout << "priors[" << i << "]= " << pname << " " << val << " " << err << endl;
	priors[pname]= new GaussianPrior(val,err);
      }
    }catch(std::exception) {
      cerr << "error in parsing arguments of priors" << endl;
      cerr << "priors must be of the form 'name value error'" << endl;
      cerr << opts.desc << endl;
      return EXIT_FAILURE;
    } 
  }

  return EXIT_SUCCESS;
  
}

