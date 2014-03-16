
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







//////////////////////////////////////////
size_t saved_npix;
size_t saved_nparams;
harp::vector_double P;
harp::vector_double Ninv;
harp::matrix_double A;
harp::matrix_double AtNinv;
harp::matrix_double Cinv;
harp::vector_double Dinv;
harp::vector_double Dquartinv;
harp::vector_double Dquart;
harp::matrix_double WDquartinv;
harp::matrix_double WDquart;
harp::matrix_double WDhalfinvWt;
harp::matrix_double WDhalfWt;
harp::matrix_double W;
harp::matrix_double temp;
harp::vector_double S;
harp::vector_double AtNinvP;
harp::vector_double Rf;
boost::numeric::ublas::vector < int > support;
//////////////////////////////////////////
#define MINEIGENVAL 1.e-40

void specex_column_norm (const harp::matrix_double & mat,  harp::vector_double & S ) {
  S.clear();
  for ( size_t i = 0; i < mat.size2(); ++i ) {
    for ( size_t j = 0; j < mat.size1(); ++j ) {
      S[ j ] += mat( j, i );
    }
  }
  // Invert
  for ( size_t i = 0; i < S.size(); ++i ) {
    S[i] = 1.0 / S[i];
  }
  return;
}
void specex_apply_norm_left(const harp::vector_double& v, harp::matrix_double & mat) {
  for ( size_t i = 0; i < mat.size2(); ++i ) {
    for ( size_t j = 0; j < mat.size1(); ++j ) {
      mat( j, i ) *= v[j];
    }
  }
  return;
}

void specex_apply_norm_right(const harp::vector_double& v, harp::matrix_double & mat) {
  for ( size_t i = 0; i < mat.size2(); ++i ) {
    for ( size_t j = 0; j < mat.size1(); ++j ) {
      mat( i, j ) *= v[j];
    }
  }
  return;
}
   
  
void extract(const harp::psf* psf, const specex::image_data& image, const specex::image_data& weight, double& rflux, double& invar, size_t spec_index, size_t lambda_index, PatchMap& patches) {
  
  double tstart;
  double tstop;
  double tstart1;
  double tstop1;
  
  //SPECEX_INFO("Extracting fiber " << spec_index << " lambda " << lambda_index);
  
  // defining data stamp
  size_t bx,by,nx,ny;
  psf->extent(spec_index,lambda_index,bx,by,nx,ny);
  
  specex::Stamp stamp(image);
  stamp.begin_i=max(stamp.begin_i,int(bx));
  stamp.begin_j=max(stamp.begin_j,int(by));
  stamp.end_i=min(stamp.end_i,int(bx+nx));
  stamp.end_j=min(stamp.end_j,int(by+ny));
  size_t nspec   = psf->n_spec();
  size_t nlambda = psf->n_lambda();
  
  // loading patches

  // book keeping
  std::map<size_t,size_t> used_patches;
  
  int spec_margin = 20; // 20
  int lambda_margin = 1000; // 1000
  
  tstart = harp::wtime();
  
  size_t param_index = 0;
  //  load neighbouring patches( compute only if not in store)
  for(size_t other_spec = size_t(max(0,int(spec_index)-spec_margin)); other_spec<min(nspec,spec_index+spec_margin+1); other_spec++) {
   for(size_t other_lambda = size_t(max(0,int(lambda_index)-lambda_margin)); other_lambda<min(nlambda,lambda_index+lambda_margin+1); other_lambda++) {
     
     
     psf->extent(other_spec,other_lambda,bx,by,nx,ny);
     
     specex::Stamp other_stamp; other_stamp.SetParent(stamp); // stamp in the coordinate system of the target stamp
     int dx = int(bx)-stamp.begin_i;
     int dy = int(by)-stamp.begin_j;
     other_stamp.begin_i=max(other_stamp.begin_i,dx);
     other_stamp.begin_j=max(other_stamp.begin_j,dy);
     dx = int(bx+nx)-stamp.begin_i;
     dy = int(by+ny)-stamp.begin_j;
     other_stamp.end_i=min(other_stamp.end_i,dx);
     other_stamp.end_j=min(other_stamp.end_j,dy);
     
     if(other_stamp.n_rows()>0 && other_stamp.n_cols()>0) {
       
       //SPECEX_INFO(other_spec << "," << other_lambda << " is a neighbour of " << spec_index << "," << lambda_index << " size= " << other_stamp.n_cols() << "x" << other_stamp.n_rows());
       
       int index = other_lambda+nlambda*other_spec;

       used_patches[index] = param_index; param_index++;

       if( patches.find(index) == patches.end()) {

	 //SPECEX_INFO("Loading patch at lambda=" << other_lambda << " spec=" << other_spec);

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
      //SPECEX_INFO("Erasing patch at " << it->first%nlambda << " "<<   it->first/nlambda);
      patches.erase(it);
      it=patches.begin(); // better start from begining
    }else{
      it++;
    }
  }
  tstop = harp::wtime();
  cout << "load/unload patches : " << tstop-tstart << " s" << endl;
  
  
  
  // fill A matrix (keep notations of spectroperf.)
  size_t nparams     = used_patches.size();
  size_t stamp_ncols = stamp.n_cols();
  size_t stamp_nrows = stamp.n_rows();
  size_t npix        = stamp_ncols*stamp_nrows;
  
    
  if(npix > saved_npix) {
    cout << "resizing npix" << endl;
    P.resize(npix);
    Ninv.resize(npix);    
  }
  if(nparams > saved_nparams) {
    cout << "resizing nparams" << endl;
    Cinv.resize(nparams,nparams);
    Dinv.resize(nparams);
    Dquartinv.resize(nparams);
    Dquart.resize(nparams);
    W.resize (nparams,nparams); 
    WDquartinv.resize(nparams,nparams);
    WDquart.resize(nparams,nparams);
    WDhalfinvWt.resize(nparams,nparams);
    WDhalfWt.resize(nparams,nparams);
    temp.resize(nparams,nparams);
    S.resize(nparams);
    AtNinvP.resize(nparams);
    Rf.resize(nparams);
    support.resize(2*nparams);
  }
  if((npix > saved_npix) || (nparams > saved_nparams)) {
    A.resize(npix,nparams);
    AtNinv.resize(nparams,npix);
  }
  
  if(npix>saved_npix) saved_npix = npix;
  if(nparams>saved_nparams) saved_nparams = nparams;
  
  P.clear(); // the image
  Ninv.clear();// the noise diagonal matrix
  
  int pixindex=0;
  for(int j=stamp.begin_j;j<stamp.end_j;j++) {
    for(int i=stamp.begin_i;i<stamp.end_i;i++,pixindex++) {
      P(pixindex)=image(i,j);
      Ninv(pixindex)=weight(i,j);
    }
  }
  
  A.clear(); // P=A*f
  
  
  tstart = harp::wtime();
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
	int patch_x = ccd_x-it->second.x_offset;
	int patch_y = ccd_y-it->second.y_offset;

	// the patch is flipped in harp
	A(pix_index,param_index) = it->second.h(patch_y,patch_x); // h is flipped in harp
	
      }
    }

  }  
  tstop = harp::wtime();
  cout << "filling A : " << tstop-tstart << " s" << endl;
  // now brute force algebra
  
  tstart = harp::wtime();
  AtNinv.clear(); // At*Ni
  for ( size_t i = 0; i < nparams; ++i ) {
    for ( size_t j = 0; j < npix; ++j ) {
      AtNinv(i,j) = A(j,i)*Ninv(j);
    }
  }
  
  Cinv.clear();
  blas::gemm(1,AtNinv,A,1,Cinv); // Cinv=At*Ni*A
  
  //specex::write_new_fits_image("Cinv.fits",Cinv);

  // eigen decomposition  
  // harp::eigen_decompose(Cinv,Dinv,W); // Cinv = W*Dinv*Wt
  tstart1=harp::wtime();
  {
    ublas::noalias(temp)=Cinv;
    int nfound;
    boost::numeric::bindings::lapack::syevr ( 'V', 'A', boost::numeric::bindings::lower ( temp ), 0.0, 0.0, 0, 0, 0.0, nfound, Dinv, W, support );
  }
  tstop1=harp::wtime();
  cout << "eigen decompose   : " << tstop1-tstart1 << " s " << endl;
  

  for(size_t i=0;i<Dinv.size();i++) {
    if(Dinv(i)<MINEIGENVAL) { 
      //cout << "fixing eigenval=" << Dinv(i) << endl; 
      Dinv(i)=MINEIGENVAL; 
    }
    Dquartinv(i)=sqrt(sqrt(Dinv(i)));
    Dquart(i)=1/Dquartinv(i);
  }
  noalias(WDquartinv)=W;
  specex_apply_norm_right(Dquartinv,WDquartinv); // WDquartinv
  harp::matrix_double& WDquart=W; // just to remember what it is
  specex_apply_norm_right(Dquart,WDquart); // WDquart
  WDhalfWt.clear(); WDhalfinvWt.clear();
  blas::gemm(1,WDquart,boost::numeric::bindings::trans(WDquart),1,WDhalfWt);
  blas::gemm(1,WDquartinv,boost::numeric::bindings::trans(WDquartinv),1,WDhalfinvWt);
  //blas::syrk(1,boost::numeric::bindings::trans(WDquart),1,WDhalfWt);
  //blas::syrk(1,boost::numeric::bindings::trans(WDquartinv),1,WDhalfinvWt);
  
  // computing S, Sij=1/sum_j Qij , where Q=W*Dinv^1/2*Wt
  specex_column_norm ( WDhalfinvWt, S );
  
  specex_apply_norm_left(S,WDhalfWt); // now it is S W D^(1/2) Wt 
  harp::matrix_double& SWDhalfWt = WDhalfWt; // just to remember
  
  AtNinvP.clear();
  blas::gemv(1,AtNinv,P,1,AtNinvP); // AtNinvP = AtNinv P
  
  Rf.clear(); // here we are, the reconvolved fluxes Rf 
  blas::gemv(1,SWDhalfWt,AtNinvP,1,Rf); // Rf = S Wt D^(1/2) W At Ninv P
  
  // collecting results
  int index_of_target_param = used_patches.find(lambda_index+nlambda*spec_index)->second;
  
  rflux=Rf(index_of_target_param);
  invar=S(index_of_target_param); invar*=invar; // Ctilde=S^-2 and we want here Ctilde^-1=S^2
  tstop = harp::wtime();
  
  cout << "algebra   : " << tstop-tstart << " s " << endl;
  // we discard all the rest

}

int main ( int argc, char *argv[] ) {

  saved_npix=0;
  saved_nparams=0;

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
  
  specex::image_data rflux(nlambda,nspec);
  specex::image_data invar(nlambda,nspec);

  PatchMap patches;
  
  for(size_t s=0;s<nspec;s++) {
    if(s%2==0) {
      for(int l=0;l<nlambda;l++) {
	extract(psf,image,weight,rflux(l,s),invar(l,s),s,l,patches);
	SPECEX_INFO("spec " << s << " lambda " << l << " flux " <<  rflux(l,s) << " error " << 1./sqrt(invar(l,s)));	
      }
    } else { // go the other way to minimize number of patches to reload
      for(int l=nlambda-1;l>=0;l--) {
	extract(psf,image,weight,rflux(l,s),invar(l,s),s,l,patches);
	SPECEX_INFO("spec " << s << " lambda " << l << " flux " <<  rflux(l,s) << " error " << 1./sqrt(invar(l,s)));	
      }
    }
  }

  { // writing output

    harp::vector_double lambda = psf->lambda();
    specex::image_data wave(nlambda,1);
    for(size_t i=0;i<nlambda;i++) wave.data(i)=lambda(i);
    
    fitsfile * fp;  
    harp::fits::create ( fp, outfilename );
    harp::fits::img_append < double > ( fp, rflux.n_rows(), rflux.n_cols() );
    harp::fits::img_write ( fp, rflux.data ,false);
    harp::fits::key_write(fp,"WHAT","RFLUX","");
    harp::fits::img_append < double > ( fp, invar.n_rows(), invar.n_cols() );
    harp::fits::img_write ( fp, invar.data ,false);
    harp::fits::key_write(fp,"EXTNAME","IVAR","");
    harp::fits::img_append < double > ( fp, 1, wave.n_cols() );
    harp::fits::img_write ( fp, wave.data ,false);
    harp::fits::key_write(fp,"EXTNAME","WAVELENGTH","");
    harp::fits::close ( fp );
  }
  return 0;
}
