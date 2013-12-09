#include <specex_hermite.h>
#include <harp.hpp>

// probabilists Hermite polynomial

double specex::HermitePol(const int Degree, const double &x) {
  switch (Degree)
    {
    case 0: return 1; break;
    case 1: return x; break;
    case 2: return x*x-1; break;
    case 3: return x*(x*x-3); break;
    case 4: return x*x*(x*x-6)+3; break;
    case 5: return x*x*x*(x*x-10)+15*x; break;
      /*
	H0(x)=1
	H1(x)=x
	H2(x)=x^2-1
	H3(x)=x^3-3x
	H4(x)=x^4-6x^2+3
	H5(x)=x^5-10x^3+15x
	H6(x)=x^6-15x^4+45x^2-15
	H7(x)=x^7-21x^5+105x^3-105x
	H8(x)=x^8-28x^6+210x^4-420x^2+105
	H9(x)=x^9-36x^7+378x^5-1260x^3+945x
	H10(x)=x^{10}-45x^8+630x^6-3150x^4+4725x^2-945 
      */    
    default : 
      HARP_THROW("ERROR consider implementing LegendreHermite pol of abitrary degree");
      break;
    }
  
}



double specex::HermitePolDerivative(const int Degree, const double &x) {
  switch (Degree)
    {
    case 0: return 0; break;
    case 1: return 1; break;
    case 2: return 2*x; break;
    case 3: return 3*x*x-3; break;
    case 4: return 4*x*(x*x-3); break;
    case 5: return 5*x*x*(x*x-6)+15; break;
      /*
	H0(x)=1
	H1(x)=x
	H2(x)=x^2-1
	H3(x)=x^3-3x
	H4(x)=x^4-6x^2+3
	H5(x)=x^5-10x^3+15x
	H6(x)=x^6-15x^4+45x^2-15
	H7(x)=x^7-21x^5+105x^3-105x
	H8(x)=x^8-28x^6+210x^4-420x^2+105
	H9(x)=x^9-36x^7+378x^5-1260x^3+945x
	H10(x)=x^{10}-45x^8+630x^6-3150x^4+4725x^2-945 
      */    
    default : 
      HARP_THROW("ERROR consider implementing LegendreHermite pol of abitrary degree");
      break;
    }
  
}
