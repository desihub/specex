
#include <unhrp.h>

#include <specex_hermite.h>

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
    case 6: return -15+x*x*(45+x*x*(-15+x*x)); break;
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
      return x*HermitePol(Degree-1,x)-(Degree-1)*HermitePol(Degree-2,x);
    }
  
}

void specex::HermitePols(unhrp::vector_double& H, const int Degree, const double &x) {
  
  H.resize(Degree+1);
  H[0]=1; if(Degree==0) return;
  H[1]=x; if(Degree==1) return;
  for(int d=2;d<=Degree;d++) {
    H[d]=x*H[d-1]-(d-1)*H[d-2];
  }
}



void specex::HermitePolsAndDerivatives(unhrp::vector_double& H, unhrp::vector_double& dHdx, const int Degree, const double &x) {
  H.resize(Degree+1);
  dHdx.resize(Degree+1);
  H[0]=1; dHdx[0]=0; if(Degree==0) return;
  H[1]=x; dHdx[1]=1; if(Degree==1) return;
  for(int d=2;d<=Degree;d++) {
    dHdx[d]=d*H[d-1];
    H[d]=x*H[d-1]-dHdx[d-1];
  }
}

double specex::HermitePolDerivative(const int Degree, const double &x) {
  if(Degree==0) return 0;
  return Degree*HermitePol(Degree-1,x);
}

