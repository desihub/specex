
#include <iostream>
#include <cstdio>

#include <boost/program_options.hpp>

#include <harp.hpp>

#include <specex_message.h>
#include <specex_image_data.h>
#include <specex_fits.h>
#include <specex_stamp.h>

#include <harp_plugin_specex.h>


using namespace std;
namespace popts = boost::program_options;

class Patch {
  public :
  size_t x_offset;
  size_t y_offset;
  harp::matrix_double h;
  specex::Stamp stamp;
  Patch() : x_offset(0),y_offset(0) {};
};

typedef std::map<size_t,Patch> PatchMap;
typedef std::map<size_t,Patch>::iterator PatchIterator;
typedef std::map<size_t,Patch>::const_iterator PatchConstIterator;


void extract(const harp::psf* psf, const specex::image_data& image, const specex::image_data& weight, double& rflux, double& invar, size_t spec_index, size_t lambda_index, PatchMap& patches) {
  
  SPECEX_INFO("Extracting fiber " << spec_index << " lambda " << lambda_index);
  
  // defining data stamp
  
  size_t bx,by,nx,ny;
  psf->extent(spec_index,lambda_index,bx,by,nx,ny);
  
  specex::Stamp stamp(image);
  stamp.begin_i=max(stamp.begin_i,int(bx));
  stamp.begin_j=max(stamp.begin_j,int(by));
  stamp.end_i=min(stamp.end_i,int(bx+nx));
  stamp.end_j=min(stamp.end_j,int(by+ny));
  
  //cout << "stamp begin " << stamp.begin_i << " " << stamp.begin_j << endl;
  //cout << "stamp end " << stamp.end_i << " " << stamp.end_j << endl;
  
  
  size_t nspec   = psf->n_spec();
  size_t nlambda = psf->n_lambda();
  

  // loading patches

  // book keeping
  std::map<size_t,size_t> used_patches;
  
  int spec_margin = 2; // 20
  int lambda_margin = 2; // 1000
  
  size_t param_index = 0;
  //  load neighbouring patches( compute only if not in store)
  for(size_t other_spec = size_t(max(0,int(spec_index)-spec_margin)); other_spec<min(nspec,spec_index+spec_margin+1); other_spec++) {
   for(size_t other_lambda = size_t(max(0,int(lambda_index)-lambda_margin)); other_lambda<min(nlambda,lambda_index+lambda_margin+1); other_lambda++) {
     
     
     psf->extent(other_spec,other_lambda,bx,by,nx,ny);
     
     specex::Stamp other_stamp; other_stamp.SetParent(stamp); // stamp in the coordinate system of the target stamp
     int dx = int(bx)-stamp.begin_i;
     int dy = int(by)-stamp.begin_j;
     //cout << "begin: " << dx << " " << dy << endl;
     other_stamp.begin_i=max(other_stamp.begin_i,dx);
     other_stamp.begin_j=max(other_stamp.begin_j,dy);
     //cout << "begin: " << other_stamp.begin_i << " " << other_stamp.begin_j << endl;
     dx = int(bx+nx)-stamp.begin_i;
     dy = int(by+ny)-stamp.begin_j;
     //cout << "end: " << dx << " " << dy << endl;
     other_stamp.end_i=min(other_stamp.end_i,dx);
     other_stamp.end_j=min(other_stamp.end_j,dy);
     //cout << "end: " << other_stamp.end_i << " " << other_stamp.end_j << endl;
     
     if(other_stamp.n_rows()>0 && other_stamp.n_cols()>0) {
       
       //SPECEX_INFO(other_spec << "," << other_lambda << " is a neighbour of " << spec_index << "," << lambda_index << " size= " << other_stamp.n_cols() << "x" << other_stamp.n_rows());
       
       int index = other_lambda+nlambda*other_spec;

       used_patches[index] = param_index; param_index++;

       if( patches.find(index) == patches.end()) {

	 //SPECEX_INFO("Loading patch at " << index%nlambda << " "<<   index/nlambda);

	 Patch patch;
	 psf->response(other_spec,other_lambda,patch.x_offset,patch.y_offset,patch.h);
	 patches[index] = patch;
       }
       
       patches.find(index)->second.stamp = other_stamp; // save stamp for filling A
       
     }
   }
  }
  

  // clear cache of unused patches to save memory
  PatchIterator it=patches.begin();
  while(it!=patches.end()) {
    
    if( used_patches.find(it->first) == used_patches.end()) {
      // SPECEX_INFO("Erasing patch at " << it->first%nlambda << " "<<   it->first/nlambda);
      
      patches.erase(it);
    }
    it++;
  }
  
  
  
  // fill A matrix (keep notations of spectroperf.)
  size_t nparams     = used_patches.size();
  size_t stamp_ncols = stamp.n_cols();
  size_t stamp_nrows = stamp.n_rows();
  size_t npix        = stamp_ncols*stamp_nrows;
  
  
  harp::vector_double P(npix); // the image
  harp::vector_double Ninv(npix); // the noise diagonal matrix
  
  int pixindex=0;
  for(int j=stamp.begin_j;j<stamp.end_j;j++) {
    for(int i=stamp.begin_i;i<stamp.end_i;i++,pixindex++) {
      P(pixindex)=image(i,j);
      Ninv(pixindex)=weight(i,j);
    }
  }
  
  harp::matrix_double A(npix,nparams); // P=A*f
  A.clear();
  
  for(PatchIterator it=patches.begin(); it!=patches.end();it++) {
    
    size_t param_index = used_patches.find(it->first)->second;
    
    // tricky part : putting the patch at the right place in matrix A
    
    const specex::Stamp & other_stamp = it->second.stamp;
    
    for(int j=other_stamp.begin_j;j<other_stamp.end_j;j++) {
      
      for(int i=other_stamp.begin_i;i<other_stamp.end_i;i++) {
	
	// this is the row-major index of the pixels where i and j are pixel coordinates in the stamp
	int pix_index   = i+j*stamp_ncols;

	// i and j are pixel coordinates in the stamp, so we have to add the coordinate of the bottom left corner in the CCD image
	int ccd_x   = i+stamp.begin_i; 
	int ccd_y   = j+stamp.begin_j;

	// the patch coordinates begins at x_offset y_offset
	
	//int x_offset = it->second.x_offset; // for gdb
	//int y_offset = it->second.y_offset; // for gdb
	//int patchnx = it->second.h.size2(); // for gdb
	//int patchny = it->second.h.size1(); // for gdb
	
	
	int patch_x = ccd_x-it->second.x_offset;
	int patch_y = ccd_y-it->second.y_offset;

	// the patch is flipped in harp
	A(pix_index,param_index) = it->second.h(patch_y,patch_x); // h is flipped in harp
	
      }
    }

  }  
  
  // now brute force algebra
  
  harp::matrix_double AtNinv(nparams,npix); // At*Ni
  for ( size_t i = 0; i < nparams; ++i ) {
    for ( size_t j = 0; j < npix; ++j ) {
      AtNinv(i,j) = A(j,i)*Ninv(j);
    }
  }
  
  harp::matrix_double Cinv(nparams,nparams);
  blas::gemm(1,AtNinv,A,1,Cinv); // Cinv=At*Ni*A
  
  // eigen decomposition
  harp::vector_double Dinv;
  harp::matrix_double W;
  
  // this routine is scary
  // harp::eigen_decompose(Cinv,Dinv,W); // Cinv = Wt*Dinv*W

  cout << "Dinv = " << Dinv << endl;

  {
    // debug : is Cinv=Wt*Dinv*W or Cinv=W*Dinv*Wt ?
    harp::matrix_double Dinvmat(nparams,nparams);
    Dinvmat.clear();
    for(size_t i=0;i<nparams;i++) Dinvmat(i,i)=Dinv(i);
    {
      harp::matrix_double DinvW(nparams,nparams);
      blas::gemm(1,Dinvmat,W,1,DinvW);
      harp::matrix_double WtDinvW(nparams,nparams);
      blas::gemm(1,boost::numeric::bindings::trans(W),DinvW,1,WtDinvW);
      harp::matrix_double tmp = Cinv - WtDinvW; 
      specex::write_new_fits_image("zero_WtDinvW.fits",tmp);
    }
    {
      harp::matrix_double DinvWt(nparams,nparams);
      blas::gemm(1,Dinvmat,boost::numeric::bindings::trans(W),1,DinvWt);
      harp::matrix_double WDinvWt(nparams,nparams);
      blas::gemm(1,W,DinvWt,1,WDinvWt);
      harp::matrix_double tmp = Cinv - WDinvWt; 
      specex::write_new_fits_image("zero_WDinvWt.fits",Cinv); // closer to zero but noisy !
    }
    exit(12);
  }

  // computing S, Sij=1/sum_j Qij , where Q=Wt*Dinv^1/2*W
  // Q is computed in harp::norm with the call harp::eigen_compose ( EIG_SQRT, Dinv, W, temp ) 
  // then harp::column_norm ( temp, S ) is called
  harp::vector_double S;
  harp::norm(Dinv,W,S);
  
  harp::matrix_double& RC = Cinv; // we don't need Cinv anymore, we keep the memory
  harp::eigen_compose(harp::EIG_INVSQRT,Dinv,W,RC ); // now it is Wt D^(1/2) W 
  harp::apply_norm(S,RC); // now it is S Wt D^(1/2) W 
  
  harp::vector_double AtNinvP(nparams);
  blas::gemv(1,AtNinv,P,1,AtNinvP); // AtNinvP = AtNinv P
  
  harp::vector_double Rf(nparams); // here we are, the reconvolved fluxes Rf 
  blas::gemv(1,RC,AtNinvP,1,Rf); // Rf = S Wt D^(1/2) W At Ninv P
  
  // collecting results
  int index_of_target_param = used_patches.find(lambda_index+nlambda*spec_index)->second;
  
  rflux=Rf(index_of_target_param);
  invar=S(index_of_target_param); invar*=invar; // Ctilde=S^-2 and we want here Ctilde^-1=S^2
  
  // we discard all the rest

}

int main ( int argc, char *argv[] ) {

  
  string jsonpar = "";
  string outfilename = "toto.fits";
  
  // Declare options 
  popts::options_description desc ( "Allowed Options" );
  
  desc.add_options()
  ( "help,h", "display usage information" )
  ( "verbose,v", "supress information printing" )
  ( "debug,d", "write out intermediate data products for debugging" )
  ( "out", popts::value<string>( &outfilename ), "root path of output files (default = \"harp_\")" )
  ( "par", popts::value < string > ( &jsonpar ), "JSON parameter file" )
  ;

  popts::variables_map vm;

  popts::store(popts::command_line_parser( argc, argv ).options(desc).run(), vm);
  
  popts::notify(vm);

  if ( ( argc < 2 ) || vm.count( "help" ) || ( ! vm.count( "par" ) ) ) {
    cerr << endl;
    cerr << desc << endl;
    return 0;
  }
  
  specex_set_verbose(vm.count("verbose")>0);
  
  SPECEX_INFO("Reading " << jsonpar);
  
  boost::property_tree::ptree par;
  boost::property_tree::json_parser::read_json ( jsonpar, par );
  
  SPECEX_INFO("Loading PSF");
  
  boost::property_tree::ptree psf_props = par.get_child ( "psf" );
  harp::psf* psf = harp::specex_psf_create(psf_props);

  
  SPECEX_INFO("Loading image");
  
  boost::property_tree::ptree image_props = par.get_child ( "image" );
  string path = image_props.get("path","");
  specex::image_data image,weight;
  specex::read_fits_images(path,image,weight);
  
  SPECEX_INFO("Starting extraction");

  size_t nspec   = psf->n_spec();
  size_t nlambda = psf->n_lambda();
  
  specex::image_data rflux(nspec,nlambda);
  specex::image_data invar(nspec,nlambda);

  PatchMap patches;
  
  for(size_t s=0;s<nspec;s++) {
    for(size_t l=0;l<nlambda;l++) {
      extract(psf,image,weight,rflux(s,l),invar(s,l),s,l,patches);
      SPECEX_INFO("spec " << s << " lambda " << l << " flux " <<  rflux(s,l) << " error " << 1./sqrt(invar(s,l)));
      if(l>3) {cout << "EXIT FOR DEBUG" << endl; exit(12);}
    }
  }
  
  

  return 0;
}
