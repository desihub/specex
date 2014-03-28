#include <iostream>
#include "specex_linalg.h"
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/basic_text_oarchive.hpp>
#include <time.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp> 
#include <boost/random.hpp>



using namespace std;

int main() {
  
  harp::vector_double x(4);
  x(1)=1;
  x[2]=2;

  if(0) {
  boost::archive::xml_oarchive xml_oa ( cout );
  xml_oa << BOOST_SERIALIZATION_NVP(x);
  
  boost::archive::text_oarchive text_oa ( cout );
  text_oa << BOOST_SERIALIZATION_NVP(x);
  
  boost::archive::text_oarchive basic_text_oa ( cout );
  basic_text_oa << BOOST_SERIALIZATION_NVP(x);
  }


  {


    //ublas::mapped_vector<double> sv (4000);
    //ublas::compressed_vector<double> sv (4000);
     cout << "==================" << endl;
  cout << " SPARSE VECTORS   " << endl;
  cout << "==================" << endl;
  
  int n = 8000;
  //ublas::coordinate_vector<double> sv (n);
  //ublas::mapped_vector<double> sv (n);
  ublas::compressed_vector<double> sv (n);
  //ublas::vector<double> sv (n);
  harp::vector_double v (n);
  harp::vector_double v2 (n);
  double w = 0.2;
  v.clear();
  v2.clear();
  sv.clear();
  
  for (unsigned i = 0; i < sv.size (); i ++) {
    v2(i)=i;
  }
  for (unsigned i = 0; i < sv.size (); i += 10) {
    sv(i) = i;
    v(i)=i;
  }
  
  {
    harp::matrix_double a(10,10);
    a.clear();
    harp::vector_double h = ublas::project(v2,ublas::range(0,10));
    cout << h << endl;
    specex::syr(1,ublas::project(h,ublas::range(0,5)),a);
    cout << a << endl;
    a.clear();
    specex::syr(1,ublas::project(h,ublas::range(5,10)),a);
    cout << a << endl;
    
    exit(12);
  }	  


    // see http://www.boost.org/doc/libs/1_38_0/libs/numeric/ublas/doc/blas.htm
  
  harp::matrix_double m(n,n);
  
  int N=500;
  {
    m.clear();
    clock_t tstart = clock();
    for(int i=0;i<N;i++)
      specex::syr(w,v,m);
    clock_t tstop = clock();
    cout << "#1 n clocks = " << tstop-tstart << " " << float(tstop-tstart)/float(CLOCKS_PER_SEC) << endl;
  }
  
  //ublas::coordinate_matrix<double> sm (n,n);
  ublas::compressed_matrix<double> sm (n,n);
  //ublas::triangular_matrix<double,  ublas::lower> sm (n,n); // crashed
  //ublas::mapped_matrix<double> sm (n,n);
  //ublas::matrix<double> sm (n,n);
  {
    
    clock_t  tstart = clock();
    //ublas::outer_prod(sv,sv);
    clock_t tstop = clock();
    cout << "#1 n clocks = " << tstop-tstart << " " << float(tstop-tstart)/float(CLOCKS_PER_SEC) << endl;
  }
  {
    
    clock_t  tstart = clock();
    //ublas::outer_prod(sv,sv);
    
    //ublas::triangular_adaptor<ublas::compressed_matrix<double>,ublas::lower> tsm(sm);
    for(int i=0;i<N;i++)
      ublas::noalias(sm) += w*ublas::outer_prod(sv,sv);
    //ublas::noalias(m) += ublas::sparse_prod(sv,sv);
      
    //ublas::blas_3::srk(m,1,1,v);
    
    clock_t tstop = clock();
    cout << "#2 n clocks = " << tstop-tstart << " " << float(tstop-tstart)/float(CLOCKS_PER_SEC) << endl;
  }
  harp::matrix_double dm = m - sm;
  for(unsigned i = 0; i < 100;i++) {
    cout << "m dm ("<< i << ") =" << m(i,i) << " " << dm(i,i) << endl;
  }
  

  if(0) {
    harp::matrix_double m(n,n); m.clear();
    clock_t tstart = clock();
    for(int i=0;i<N;i++)
      specex::syr(w,v2,m);
    clock_t tstop = clock();
    cout << "#3 n clocks = " << tstop-tstart << " " << float(tstop-tstart)/float(CLOCKS_PER_SEC) << endl;
  }
  
    return 0;

}
  cout << "==================" << endl;
  cout << " LINEAR SYSTEMS   " << endl;
  cout << "==================" << endl;
  
  typedef boost::ecuyer1988 base_generator_type;
  base_generator_type generator(42u);
  boost::uniform_real < double > dist ( 0.0, 1.0 );
  boost::variate_generator < base_generator_type&, boost::uniform_real < double > > uni ( generator, dist );

  int ndata = 12;
  int nparams = 10;
  vector<harp::vector_double> h;
  for( int i=0; i<ndata; i++) {
    harp::vector_double hi(nparams);
    for(int p=0; p<nparams ;p++) {
      hi(p)=uni(); // stupid number
    }
    //cout << "data " << i << " h=" << hi << endl;
    h.push_back(hi);
  }
  
  harp::vector_double true_params(nparams);
  for(int p=0; p<nparams ;p++) {
    true_params(p)=p;
  }
  
  harp::matrix_double A(nparams,nparams);
  harp::vector_double B(nparams);
  A.clear();
  B.clear();
  
  harp::vector_double data(ndata);
  for( int i=0; i<ndata; i++) {
    data[i]=specex::dot(true_params,h[i]);
    double w=uni()+1; // stupid number
    specex::syr(w,h[i],A);
    specex::axpy(w*data(i),h[i],B);
  }
  
  //cout << "A=" << A << endl;
  //cout << "B=" << B << endl;
  

  int status = specex::cholesky_solve(A,B);
  cout << "cholesky_solve status = " << status << endl;
  harp::vector_double diff = B-true_params;
  cout << "parameter residuals = " << diff << endl;
  harp::vector_double res(ndata);
  for( int i=0; i<ndata; i++) {
    res(i) = data(i)-specex::dot(B,h[i]);
  }
  cout << "data residuals = " << res << endl;
  
  

  cout << "==================" << endl;
  cout << " PROJECTIONS   " << endl;
  cout << "==================" << endl;
  
  harp::vector_double one(16);
  harp::vector_double zero(16);
  zero.clear();
  
  for(size_t i=0;i<one.size();i++) one[i]=i+1;
  for(int j=0;j<4;j++)
    ublas::project(zero,ublas::range(4*j,4*(j+1))) += pow(10,j)* ublas::project(one,ublas::range(4*j,4*(j+1)));
  cout << zero << endl;

  return 0;
  harp::vector_double y(12);
  for(int i=7;i<7+4;i++) y[i]=1;

  try {
    //cout << specex::dot(x,y) << endl; // this throws an exception
    
    
    ublas::range r(0,1);
    cout << specex::dot(x,ublas::project(y,ublas::range(0,x.size()))) << endl;
    cout << specex::dot(x,ublas::project(y,ublas::range(7,7+x.size()))) << endl;


    y.clear();
    size_t x_size = x.size();
    for(size_t k=0;(k+1)*x_size<=y.size();k++) {
      ublas::project(y,ublas::range(k*x_size,(k+1)*x_size)) = (k+1)*x;
    }

    cout << y << endl;
    cout << ublas::inner_prod (y,y) << endl;
    cout << specex::dot(y,y) << endl;
    cout << ublas::outer_prod (y,y) << endl;



    cout << "==================" << endl;
    cout << " SPARSE VECTORS   " << endl;
    cout << "==================" << endl;

    //ublas::mapped_vector<double> sv (8000);
    ublas::compressed_vector<double> sv (8000);
    //ublas::coordinate_vector<double> sv (8000);
    harp::vector_double v (8000);
    
   

    if(0) {
      
      cout << sv << endl;
      // for coordinate vector
      //for(ublas::coordinate_vector<double>::iterator it=sv.begin();it!=sv.end();++it) {
      //cout << *it << " ";
      //}
      //cout << endl;
    }
    
    double res1=0;
    {
      clock_t tstart = clock();
      ublas::vector<double> v2(v.size());
      ublas::matrix<double> m = ublas::outer_prod(v,v); boost::numeric::bindings::blas::gemv(1,m,v,0,v2);
      //boost::numeric::bindings::blas::gemv(1,(ublas::matrix<double>&)ublas::outer_prod(v,v),v,0,v2);
      //cout << v2 << endl;
      clock_t tstop = clock();
      res1=ublas::inner_prod(v,v2);
      cout << res1 << endl;
      cout << "n clocks = " << tstop-tstart << " " << float(tstop-tstart)/float(CLOCKS_PER_SEC) << endl;
    }
    if(1){ // ok 
      clock_t tstart = clock();
      ublas::vector<double> v2(sv.size());
      ublas::matrix<double> m = ublas::outer_prod(sv,sv); boost::numeric::bindings::blas::gemv(1,m,v,0,v2);
      double res2 = ublas::inner_prod(sv,v2);
      clock_t tstop = clock();
      cout << res2 << " " << res2/res1-1 << endl;
      cout << "n clocks = " << tstop-tstart << " " << float(tstop-tstart)/float(CLOCKS_PER_SEC) << endl;
    }
    if(1){ // ok 
      clock_t tstart = clock();
      double res2 = ublas::inner_prod(sv,ublas::prod(ublas::outer_prod(sv,sv),sv));
      clock_t tstop = clock();
      cout << res2 << " " << res2/res1-1 << endl;
      cout << "n clocks = " << tstop-tstart << " " << float(tstop-tstart)/float(CLOCKS_PER_SEC) << endl;
    }
    
    // same with "pre-existing" matrix
    { 
      ublas::matrix<double, ublas::column_major> m(v.size(),v.size());
      clock_t tstart = clock();
      //m += ublas::outer_prod(sv,sv); // terribly slow
      // ublas::add(m,ublas::outer_prod(sv,sv)); // doesn't compile 
      boost::numeric::bindings::blas::syr(1.,v,boost::numeric::bindings::lower(m)); // compiles, and fast, but fills only part of the matrix !!
      //boost::numeric::bindings::blas::syr(1.,sv,boost::numeric::bindings::lower(m)); // doesn't compile
      
      //specex::syr(1.,v,m);
      double res2 = ublas::inner_prod(sv,ublas::prod(m,sv));
      clock_t tstop = clock();
      cout << res2 << " " << res2/res1-1 << endl;
      cout << "n clocks = " << tstop-tstart << " " << float(tstop-tstart)/float(CLOCKS_PER_SEC) << endl;
    }


  }catch(harp::exception e) {
    cout << "ohoh " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
