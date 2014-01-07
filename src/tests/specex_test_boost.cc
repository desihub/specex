#include <iostream>
#include "specex_linalg.h"
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/basic_text_oarchive.hpp>
#include <time.h>
#include <boost/numeric/ublas/matrix.hpp>
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
  A *= 0;
  B *= 0;
  
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
  zero *= 0;
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


    y *= 0;
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

    //ublas::mapped_vector<double> sv (4000);
    //ublas::compressed_vector<double> sv (4000);
    ublas::coordinate_vector<double> sv (4000);
    ublas::vector<double> v (4000);
    
    for (unsigned i = 0; i < sv.size (); i += 20) {
        sv (i) = i;
	v(i)=i;
    }
    
    
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
