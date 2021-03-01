#ifndef SPECEX_STAMP__H
#define SPECEX_STAMP__H

//#include <harp.hpp>
#include <harp/image.hpp>

//#include "specex_psf.h"
#include <iostream>

namespace specex {

class Stamp {
 public :
  // current bounds of psf fit stamp in image
  const harp::image* parent_image; // either one of the two
  const Stamp* parent_stamp; // either one of the two
  int begin_i;
  int begin_j;
  int end_i; 
  int end_j;
  
 
 
  void SetParent(const harp::image& image) {
    parent_image = &image;
    parent_stamp = 0;
    begin_i=0;
    end_i=parent_image->n_cols();
    begin_j=0;
    end_j=parent_image->n_rows();
  }
 
 void SetParent(const Stamp& other) {
   parent_image = 0;
   parent_stamp = &other;
   begin_i=0;
   end_i=parent_stamp->n_cols();
   begin_j=0;
   end_j=parent_stamp->n_rows();
 }
 
 
 Stamp() {
   parent_image = 0;
   parent_stamp = 0;
   begin_i=0;
   end_i=0;
   begin_j=0;
   end_j=0;
 }
 
 Stamp(const harp::image& image) {
   SetParent(image);
 }
 
 int Parent_n_cols() const {
   if(parent_image) return parent_image->n_cols();
   return parent_stamp->n_cols();
 }
 int Parent_n_rows() const {
   if(parent_image) return parent_image->n_rows();
   return parent_stamp->n_rows();
 }

 bool Contains(int i,int j) const {
   return(i>=begin_i && i<end_i && j>=begin_j && j<end_j);
 }
 
 /*
 void SetLimitsFromPSF(const SpecExPSF& PSF, const double &X, const double &Y) {
   PSF.StampLimits(X,Y,begin_i,end_i,begin_j,end_j);
   // now check image bounds
   begin_i = max(0,begin_i);
   end_i   = min(Parent_n_cols(),end_i);
   begin_j = max(0,begin_j);
   end_j   = min(Parent_n_rows(),end_j);
 }
 
 void SetLimitsFromPSF(const SpecExPSF& PSF, 
		       const double &spot_x_min, const double &spot_x_max,
		       const double &spot_y_min, const double &spot_y_max) {
   
   int k,p;
   PSF.StampLimits(spot_x_min,spot_y_min,begin_i,k,begin_j,p);
   PSF.StampLimits(spot_x_max,spot_y_max,k,end_i,p,end_j);
   
   // now check image bounds
   begin_i = max(0,begin_i);
   end_i   = min(Parent_n_cols(),end_i);
   begin_j = max(0,begin_j);
   end_j   = min(Parent_n_rows(),end_j);
 }
 */
 
  int n_cols() const {return end_i-begin_i;}
  int n_rows() const {return end_j-begin_j;}

  
  Stamp Intersection(const Stamp& other_stamp) {
    if((parent_image !=  other_stamp.parent_image) || (parent_stamp !=  other_stamp.parent_stamp) ) {
      HARP_THROW("Stamp::Intersection not same parents");
    }
    Stamp inter(*this);
    if(inter.begin_i<other_stamp.begin_i) inter.begin_i=other_stamp.begin_i;
    if(inter.end_i>other_stamp.end_i) inter.end_i=other_stamp.end_i;
    if(inter.begin_j<other_stamp.begin_j) inter.begin_j=other_stamp.begin_j;
    if(inter.end_j>other_stamp.end_j) inter.end_j=other_stamp.end_j;
    return inter;
  }

  void write(std::ostream& os) const {
    os << "[" << begin_i << ":" << end_i << "," << begin_j << ":" << end_j << "]";
  }
  friend std::ostream& operator << (std::ostream &stream, const specex::Stamp &stamp);
  
  
};
}

#endif
