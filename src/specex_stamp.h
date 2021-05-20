#ifndef SPECEX_STAMP__H
#define SPECEX_STAMP__H

#include <harp/image.hpp>
#include <iostream>
#include <unhrp_image.hpp>

namespace specex {

class Stamp {
 public :
  // current bounds of psf fit stamp in image
  const unhrp::image* parent_image; // either one of the two
  const Stamp* parent_stamp; // either one of the two
  int begin_i;
  int begin_j;
  int end_i; 
  int end_j;
  
 
 
  void SetParent(const unhrp::image& image) {
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
 
 Stamp(const unhrp::image& image) {
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
 
  int n_cols() const {return end_i-begin_i;}
  int n_rows() const {return end_j-begin_j;}

  
  Stamp Intersection(const Stamp& other_stamp) {
    if((parent_image !=  other_stamp.parent_image) || (parent_stamp !=  other_stamp.parent_stamp) ) {
      SPECEX_ERROR("Stamp::Intersection not same parents");
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
